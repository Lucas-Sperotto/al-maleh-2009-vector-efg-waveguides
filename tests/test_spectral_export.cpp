// SPDX-License-Identifier: MIT
#include "al_maleh/spectral.hpp"
#include <zlib.h>
#include <Eigen/Core>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace al_maleh;
namespace {
void require(bool condition,const char* text) {
  if(!condition)throw std::runtime_error(text);
}
void must_fail(const std::function<void()>& operation,const char* label) {
  try{operation();}catch(const std::exception&){return;}
  throw std::runtime_error(std::string(label)+" should throw");
}
std::string uncompress(const std::filesystem::path& path) {
  gzFile file=gzopen(path.string().c_str(),"rb");
  require(file!=nullptr,"cannot reopen gzip artifact");
  char buffer[8192];std::string result;int count=0;
  while((count=gzread(file,buffer,sizeof buffer))>0)result.append(buffer,static_cast<std::size_t>(count));
  const int close_status=gzclose(file);
  require(count==0&&close_status==Z_OK,"invalid gzip artifact");return result;
}
Eigen::MatrixXd read_matrix(const std::filesystem::path& path) {
  std::istringstream input(uncompress(path));std::string line;
  std::getline(input,line);
  const bool array=line=="%%MatrixMarket matrix array real general";
  require(array||line=="%%MatrixMarket matrix coordinate real general","bad Matrix Market header");
  do{std::getline(input,line);}while(!line.empty()&&line[0]=='%');
  std::istringstream shape(line);Eigen::Index rows=0,cols=0,count=0;
  shape>>rows>>cols;require(rows>0&&cols>0,"bad Matrix Market dimensions");
  Eigen::MatrixXd result=Eigen::MatrixXd::Zero(rows,cols);
  if(array) {
    for(Eigen::Index j=0;j<cols;++j)for(Eigen::Index i=0;i<rows;++i)
      require(static_cast<bool>(input>>result(i,j)),"incomplete array matrix");
  } else {
    require(static_cast<bool>(shape>>count),"missing coordinate nonzero count");
    for(Eigen::Index k=0;k<count;++k) {
      Eigen::Index i=0,j=0;double value=0;
      require(static_cast<bool>(input>>i>>j>>value)&&i>0&&i<=rows&&j>0&&j<=cols,"bad coordinate entry");
      result(i-1,j-1)+=value;
    }
  }
  std::string extra;require(!(input>>extra),"unexpected trailing matrix entries");
  require(result.allFinite(),"nonfinite exported matrix");return result;
}
std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> cells;std::istringstream row(line);std::string cell;
  while(std::getline(row,cell,','))cells.push_back(cell);return cells;
}
void round_trip(const std::filesystem::path& dir) {
  const Rectangle domain{0,2,0,1};
  const MLSEvaluator mls(make_regular_cloud(domain,5,5,1.65),1);
  const BackgroundQuadrature rule{domain,4,4,3};
  const auto assembled=assemble_cartesian(mls,rule);
  const auto spectrum=solve_generalized_symmetric(assembled.stiffness,assembled.mass);
  SpectralExportOptions options;options.nx=5;options.ny=5;options.support_factor=3.3;
  options.grid_nx=7;options.grid_ny=5;options.origin="test_\"quoted\"_origin\n";
  export_spectrum(dir,mls,rule,assembled,spectrum,options);
  const auto V=read_matrix(dir/"eigenvectors.mtx.gz");
  const auto K=read_matrix(dir/"stiffness_matrix.mtx.gz");
  const auto M=read_matrix(dir/"mass_matrix.mtx.gz");
  require(V.rows()==50&&V.cols()==50,"full eigenvector spectrum not exported");
  require((V-spectrum.eigenvectors).cwiseAbs().maxCoeff()==0,"17-digit eigenvector roundtrip changed coefficients");
  require((K-Eigen::MatrixXd(assembled.stiffness)).norm()==0,"stiffness roundtrip changed coefficients");
  require((M-Eigen::MatrixXd(assembled.mass)).norm()==0,"mass roundtrip changed coefficients");
  require((V.transpose()*M*V-Eigen::MatrixXd::Identity(50,50)).cwiseAbs().maxCoeff()<1e-6,
          "reopened vectors are not M-orthonormal");
  std::ifstream eigenvalues(dir/"eigenvalues.csv");std::string line;
  std::getline(eigenvalues,line);
  require(line=="index,lambda,kc,status,relative_residual","legacy spectrum schema changed");
  int modes=0,near_zero=0;
  while(std::getline(eigenvalues,line)) {
    const auto row=split(line);require(row.size()==5,"bad eigenvalue record");
    const int index=std::stoi(row[0]);const double lambda=std::stod(row[1]);
    require(index==modes&&index<V.cols(),"eigenvalue index does not match vector column");
    require(lambda==spectrum.eigenvalues[index],"eigenvalue roundtrip lost precision");
    const double residual=(K*V.col(index)-lambda*M*V.col(index)).norm()/
      ((K.norm()+std::abs(lambda)*M.norm())*V.col(index).norm());
    require(residual<1e-9&&std::abs(residual-std::stod(row[4]))<1e-10,"reopened eigenpair residual failed");
    if(row[3]=="near_zero")++near_zero;
    ++modes;
  }
  require(modes==50&&near_zero==spectrum.numerical_kernel_dimension,"kernel columns discarded or misclassified");
  std::ifstream dofs(dir/"dof_map.csv");std::getline(dofs,line);int dof=0;
  while(std::getline(dofs,line)) {
    const auto row=split(line);require(row.size()==6,"bad DOF map record");
    const auto& node=mls.nodes()[static_cast<std::size_t>(dof%25)];
    require(std::stoi(row[0])==dof&&std::stoi(row[1])==dof%25&&row[2]==(dof<25?"Hx":"Hy"),"DOF ordering is not Hx then Hy");
    require(std::stod(row[3])==node.position.x&&std::stod(row[4])==node.position.y&&
            std::stod(row[5])==node.support_radius,"DOF geometry mismatch");++dof;
  }
  require(dof==50,"incomplete DOF map");
  const auto curls=read_matrix(dir/"field_curl.mtx.gz");
  require(curls.rows()==35&&curls.cols()==50,"incomplete common-grid curl matrix");
  std::ifstream grid(dir/"grid.csv");std::getline(grid,line);int grid_row=0;
  while(std::getline(grid,line)) {
    const auto row=split(line);require(row.size()==4&&std::stoi(row[0])==grid_row,"bad grid row index");
    const Point2 p{std::stod(row[1]),std::stod(row[2])};
    require(std::abs(std::stod(row[3])-2.0/35)<1e-16,"incorrect common-grid weight");
    const auto shape=mls.evaluate(p);
    for(Eigen::Index column:{Eigen::Index(0),Eigen::Index(25),Eigen::Index(49)}) {
      double value=0;
      for(std::size_t k=0;k<shape.active_nodes.size();++k) {
        const auto n=static_cast<Eigen::Index>(shape.active_nodes[k]);
        value+=shape.dx[k]*V(25+n,column)-shape.dy[k]*V(n,column);
      }
      require(std::abs(value-curls(grid_row,column))<1e-9,"common-grid field column mismatch");
    }
    ++grid_row;
  }
  require(grid_row==35,"incomplete common physical grid");
  std::ifstream diagnostics(dir/"mode_diagnostics.csv");std::getline(diagnostics,line);int diagnostic_rows=0;
  while(std::getline(diagnostics,line)) {
    const auto row=split(line);require(row.size()==13,"bad full-spectrum diagnostic record");
    const int index=std::stoi(row[0]);
    require(index==diagnostic_rows&&std::abs(std::stod(row[8])-1)<1e-5,"diagnostic index or field normalization failed");
    if(spectrum.eigenvalues[index]>spectrum.zero_cutoff) {
      const auto direct=modal_diagnostics(mls,rule,V.col(index),spectrum.eigenvalues[index]);
      require(std::abs(std::stod(row[3])-direct.divergence_ratio)<1e-8&&
              std::abs(std::stod(row[4])-direct.normal_ratio)<1e-8&&
              std::abs(std::stod(row[5])-direct.curl_l2)<1e-8,"cached diagnostics differ from independent single-column evaluation");
    } else require(row[3].empty(),"near-zero divergence ratio should be undefined");
    ++diagnostic_rows;
  }
  require(diagnostic_rows==50,"not all modes have diagnostics");
  // Explicitly exercise a signed numerical zero, independently of whether the
  // eigensolver happened to produce negative roundoff on this platform.
  Spectrum signed_zero=spectrum;signed_zero.eigenvalues[0]=-1e-15;
  export_spectrum(dir/"signed_zero",mls,rule,assembled,signed_zero,options);
  std::ifstream negative(dir/"signed_zero/eigenvalues.csv");std::getline(negative,line);std::getline(negative,line);
  const auto first=split(line);
  require(std::stod(first[1])==-1e-15&&first[3]=="near_zero","negative numerical zero discarded or rounded");
  auto bad=options;bad.nx=4;
  must_fail([&]{export_spectrum(dir/"bad",mls,rule,assembled,spectrum,bad);},"inconsistent cloud dimensions");
  bad=options;bad.grid_nx=0;
  must_fail([&]{export_spectrum(dir/"bad",mls,rule,assembled,spectrum,bad);},"empty physical grid");
  auto malformed=spectrum;malformed.eigenvectors.conservativeResize(50,49);
  must_fail([&]{export_spectrum(dir/"bad",mls,rule,assembled,malformed,options);},"incomplete eigenvector matrix");
  malformed=spectrum;malformed.residuals.pop_back();
  must_fail([&]{export_spectrum(dir/"bad",mls,rule,assembled,malformed,options);},"incomplete residual vector");
}
}
int main() {
  const auto token=std::chrono::high_resolution_clock::now().time_since_epoch().count();
  const auto dir=std::filesystem::temp_directory_path()/("al_maleh_spectral_export_"+std::to_string(token));
  try {
    round_trip(dir);std::filesystem::remove_all(dir);
    std::cout<<"PASS: complete gzip/Matrix Market roundtrip, indexed eigenpairs, M-orthogonality, signed kernel, fields and diagnostics\n";
  }catch(const std::exception& error) {
    std::cerr<<"FAIL: "<<error.what()<<" (artifacts: "<<dir.string()<<")\n";return 1;
  }
}
