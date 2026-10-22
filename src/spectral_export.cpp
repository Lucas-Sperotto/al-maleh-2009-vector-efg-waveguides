// SPDX-License-Identifier: MIT
#include "al_maleh/spectral.hpp"
#include <zlib.h>
#include <Eigen/SparseCore>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef AL_MALEH_GIT_COMMIT
#define AL_MALEH_GIT_COMMIT "unknown"
#endif
#ifndef AL_MALEH_BUILD_TYPE
#define AL_MALEH_BUILD_TYPE "unknown"
#endif
#ifndef AL_MALEH_COMPILER
#define AL_MALEH_COMPILER "unknown"
#endif

namespace al_maleh {
namespace {
using RowSparse=Eigen::SparseMatrix<double,Eigen::RowMajor>;
using Triplet=Eigen::Triplet<double>;
std::ofstream output(const std::filesystem::path& path) {
  std::ofstream out(path);
  if(!out)throw std::runtime_error("Cannot create "+path.string());
  out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<std::setprecision(17);
  return out;
}
class GzipWriter {
 public:
  explicit GzipWriter(const std::filesystem::path& path):file_(gzopen(path.string().c_str(),"wb6")) {
    if(!file_)throw std::runtime_error("Cannot create compressed matrix "+path.string());
  }
  ~GzipWriter(){if(file_)gzclose(file_);}
  void write(const std::string& data) {
    const auto count=static_cast<unsigned int>(data.size());
    if(data.size()>std::numeric_limits<unsigned int>::max()||
       gzwrite(file_,data.data(),count)!=static_cast<int>(count))
      throw std::runtime_error("Cannot write compressed matrix");
  }
  void close() {
    gzFile file=file_;file_=nullptr;
    if(gzclose(file)!=Z_OK)throw std::runtime_error("Cannot close compressed matrix");
  }
 private:
  gzFile file_;
};
void matrix(const std::filesystem::path& file,const Eigen::MatrixXd& values) {
  GzipWriter out(file);
  out.write("%%MatrixMarket matrix array real general\n% Columns correspond to zero-based eigenvalues.csv index + 1.\n"+
    std::to_string(values.rows())+" "+std::to_string(values.cols())+"\n");
  for(Eigen::Index col=0;col<values.cols();++col) {
    std::ostringstream text;text<<std::setprecision(17);
    for(Eigen::Index row=0;row<values.rows();++row)text<<values(row,col)<<'\n';
    out.write(text.str());
  }
  out.close();
}
void matrix(const std::filesystem::path& file,const Eigen::SparseMatrix<double>& values) {
  GzipWriter out(file);
  out.write("%%MatrixMarket matrix coordinate real general\n% Full stored matrix; one-based row and column indices.\n"+
    std::to_string(values.rows())+" "+std::to_string(values.cols())+" "+
    std::to_string(values.nonZeros())+"\n");
  for(Eigen::Index col=0;col<values.outerSize();++col) {
    std::ostringstream text;text<<std::setprecision(17);
    for(Eigen::SparseMatrix<double>::InnerIterator value(values,col);value;++value)
      text<<value.row()+1<<' '<<value.col()+1<<' '<<value.value()<<'\n';
    out.write(text.str());
  }
  out.close();
}
std::string json_string(const std::string& value) {
  std::ostringstream out;out<<'"';
  for(unsigned char c:value) {
    if(c=='"'||c=='\\')out<<'\\'<<c;
    else if(c=='\n')out<<"\\n";
    else if(c=='\r')out<<"\\r";
    else if(c=='\t')out<<"\\t";
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<int>(c)<<std::dec;
    else out<<c;
  }
  out<<'"';return out.str();
}
std::string utc_now() {
  const std::time_t now=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  const std::tm date=*std::gmtime(&now);
  std::ostringstream out;out<<std::put_time(&date,"%Y-%m-%dT%H:%M:%SZ");return out.str();
}
const char* status(double lambda,double cutoff) {
  return lambda<=cutoff?"near_zero":"positive_candidate";
}
void finite_matrix(const Eigen::SparseMatrix<double>& matrix) {
  for(Eigen::Index col=0;col<matrix.outerSize();++col)
    for(Eigen::SparseMatrix<double>::InnerIterator value(matrix,col);value;++value)
      if(!std::isfinite(value.value()))throw std::invalid_argument("Nonfinite exported K/M");
}
RowSparse sparse(Eigen::Index rows,Eigen::Index cols,const std::vector<Triplet>& values) {
  RowSparse out(rows,cols);out.setFromTriplets(values.begin(),values.end());return out;
}
struct Sampling {
  RowSparse hx,hy,div,curl,normal;
  Eigen::VectorXd weights,boundary_weights;
};
Sampling sampling(const MLSEvaluator& mls,const BackgroundQuadrature& rule) {
  const auto n=static_cast<Eigen::Index>(mls.nodes().size());
  std::vector<Triplet> hx,hy,div,curl,normal;
  std::vector<double> weights,boundary_weights;
  for_each_quadrature_point(rule,[&](Point2 p,double w) {
    const auto row=static_cast<Eigen::Index>(weights.size());weights.push_back(w);
    const auto shape=mls.evaluate(p);
    for(std::size_t k=0;k<shape.active_nodes.size();++k) {
      const auto i=static_cast<Eigen::Index>(shape.active_nodes[k]);
      hx.emplace_back(row,i,shape.phi[k]);hy.emplace_back(row,n+i,shape.phi[k]);
      div.emplace_back(row,i,shape.dx[k]);div.emplace_back(row,n+i,shape.dy[k]);
      curl.emplace_back(row,i,-shape.dy[k]);curl.emplace_back(row,n+i,shape.dx[k]);
    }
  });
  const auto r=rule.domain;
  const double points[3]={-0.7745966692414834,0,0.7745966692414834};
  const double gauss_weights[3]={5.0/9.0,8.0/9.0,5.0/9.0};
  const auto boundary=[&](Point2 p,double w,Eigen::Index offset) {
    const auto row=static_cast<Eigen::Index>(boundary_weights.size());boundary_weights.push_back(w);
    const auto shape=mls.evaluate(p);
    for(std::size_t k=0;k<shape.active_nodes.size();++k)
      normal.emplace_back(row,offset+static_cast<Eigen::Index>(shape.active_nodes[k]),shape.phi[k]);
  };
  for(int cell=0;cell<rule.cells_y;++cell)for(int q=0;q<3;++q) {
    const double h=(r.ymax-r.ymin)/rule.cells_y;
    const double y=r.ymin+(cell+0.5*(1+points[q]))*h;
    boundary({r.xmin,y},0.5*h*gauss_weights[q],0);
    boundary({r.xmax,y},0.5*h*gauss_weights[q],0);
  }
  for(int cell=0;cell<rule.cells_x;++cell)for(int q=0;q<3;++q) {
    const double h=(r.xmax-r.xmin)/rule.cells_x;
    const double x=r.xmin+(cell+0.5*(1+points[q]))*h;
    boundary({x,r.ymin},0.5*h*gauss_weights[q],n);
    boundary({x,r.ymax},0.5*h*gauss_weights[q],n);
  }
  Sampling out;
  const auto rows=static_cast<Eigen::Index>(weights.size());
  const auto brows=static_cast<Eigen::Index>(boundary_weights.size());
  out.hx=sparse(rows,2*n,hx);out.hy=sparse(rows,2*n,hy);
  out.div=sparse(rows,2*n,div);out.curl=sparse(rows,2*n,curl);
  out.normal=sparse(brows,2*n,normal);
  out.weights=Eigen::Map<Eigen::VectorXd>(weights.data(),rows);
  out.boundary_weights=Eigen::Map<Eigen::VectorXd>(boundary_weights.data(),brows);
  return out;
}
Eigen::RowVectorXd energies(const RowSparse& op,const Eigen::MatrixXd& vectors,
                           const Eigen::VectorXd& weights) {
  const Eigen::MatrixXd fields=op*vectors;
  return (fields.array().square().colwise()*weights.array()).colwise().sum();
}
}

void export_spectrum(const std::filesystem::path& directory,
    const MLSEvaluator& mls,const BackgroundQuadrature& rule,
    const CartesianAssembly& matrices,const Spectrum& s,
    const SpectralExportOptions& options) {
  const auto nodes=static_cast<Eigen::Index>(mls.nodes().size()),ndof=2*nodes;
  const auto r=rule.domain;
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  if(options.nx<2||options.ny<2||static_cast<long long>(options.nx)*options.ny!=nodes||
     !std::isfinite(options.support_factor)||options.support_factor<=0||
     options.grid_nx<1||options.grid_ny<1||options.grid_nx>1000||options.grid_ny>1000||
     !std::isfinite(r.xmin)||!std::isfinite(r.xmax)||!std::isfinite(r.ymin)||
     !std::isfinite(r.ymax)||a<=0||b<=0||rule.cells_x<1||rule.cells_y<1||
     (rule.gauss_order<1||rule.gauss_order>5)||
     s.eigenvalues.size()!=ndof||s.eigenvectors.rows()!=ndof||s.eigenvectors.cols()!=ndof||
     matrices.mass.rows()!=ndof||matrices.mass.cols()!=ndof||
     matrices.stiffness.rows()!=ndof||matrices.stiffness.cols()!=ndof||
     s.residuals.size()!=static_cast<std::size_t>(ndof)||
     !s.eigenvalues.allFinite()||!s.eigenvectors.allFinite()||
     !std::isfinite(s.zero_cutoff)||s.zero_cutoff<0||
     !std::isfinite(s.mass_condition)||!std::isfinite(s.min_mass_eigenvalue))
    throw std::invalid_argument("Inconsistent full-spectrum export inputs");
  finite_matrix(matrices.mass);finite_matrix(matrices.stiffness);
  for(double residual:s.residuals)
    if(!std::isfinite(residual)||residual<0)throw std::invalid_argument("Invalid exported residual");
  // Cache MLS values and derivatives once per integration point, then evaluate
  // all columns together. Cost does not repeat MLS reconstruction per mode.
  const auto cache=sampling(mls,rule);
  const auto field2=(energies(cache.hx,s.eigenvectors,cache.weights)+
                     energies(cache.hy,s.eigenvectors,cache.weights)).eval();
  const auto curl2=energies(cache.curl,s.eigenvectors,cache.weights);
  const auto div2=energies(cache.div,s.eigenvectors,cache.weights);
  const auto normal2=energies(cache.normal,s.eigenvectors,cache.boundary_weights);
  if(!field2.allFinite()||!curl2.allFinite()||!div2.allFinite()||!normal2.allFinite()||
     field2.minCoeff()<=0)
    throw std::runtime_error("Invalid full-spectrum field diagnostics");
  std::filesystem::create_directories(directory);
  matrix(directory/"eigenvectors.mtx.gz",s.eigenvectors);
  matrix(directory/"mass_matrix.mtx.gz",matrices.mass);
  matrix(directory/"stiffness_matrix.mtx.gz",matrices.stiffness);
  {
    auto csv=output(directory/"eigenvalues.csv");
    csv<<"index,lambda,kc,status,relative_residual\n";
    for(Eigen::Index i=0;i<ndof;++i)
      csv<<i<<','<<s.eigenvalues[i]<<','<<std::sqrt(std::max(0.0,s.eigenvalues[i]))<<','
         <<status(s.eigenvalues[i],s.zero_cutoff)<<','<<s.residuals[static_cast<std::size_t>(i)]<<'\n';
    csv.close();
  }
  {
    auto csv=output(directory/"dof_map.csv");csv<<"dof,node,component,x,y,support_radius\n";
    for(Eigen::Index dof=0;dof<ndof;++dof) {
      const auto node=dof%nodes;const auto& p=mls.nodes()[static_cast<std::size_t>(node)];
      csv<<dof<<','<<node<<','<<(dof<nodes?"Hx":"Hy")<<','<<p.position.x<<','
         <<p.position.y<<','<<p.support_radius<<'\n';
    }
    csv.close();
  }
  {
    auto csv=output(directory/"mode_diagnostics.csv");
    csv<<"index,lambda,residual,divergence_ratio,normal_ratio,curl_norm_l2,field_norm_l2,"
         "div_norm_l2,field_energy,curl_energy,divergence_energy,normal_boundary_energy,status\n";
    for(Eigen::Index i=0;i<ndof;++i) {
      csv<<i<<','<<s.eigenvalues[i]<<','<<s.residuals[static_cast<std::size_t>(i)]<<',';
      if(s.eigenvalues[i]>s.zero_cutoff)csv<<std::sqrt(div2[i]/(s.eigenvalues[i]*field2[i]));
      csv<<','<<std::sqrt(std::max(0.0,normal2[i]*a*b/(2*(a+b)*field2[i])))<<','
         <<std::sqrt(curl2[i])<<','<<std::sqrt(field2[i])<<','<<std::sqrt(div2[i])<<','
         <<field2[i]<<','<<curl2[i]<<','<<div2[i]<<','<<normal2[i]<<','
         <<status(s.eigenvalues[i],s.zero_cutoff)<<'\n';
    }
    csv.close();
  }
  std::vector<Triplet> grid_curl;
  const Eigen::Index grid_rows=static_cast<Eigen::Index>(options.grid_nx)*options.grid_ny;
  {
    auto csv=output(directory/"grid.csv");csv<<"index,x,y,weight\n";
    for(int j=0;j<options.grid_ny;++j)for(int i=0;i<options.grid_nx;++i) {
      const Eigen::Index row=static_cast<Eigen::Index>(j)*options.grid_nx+i;
      const Point2 p{r.xmin+(i+0.5)*a/options.grid_nx,r.ymin+(j+0.5)*b/options.grid_ny};
      csv<<row<<','<<p.x<<','<<p.y<<','<<a*b/static_cast<double>(grid_rows)<<'\n';
      const auto shape=mls.evaluate(p);
      for(std::size_t k=0;k<shape.active_nodes.size();++k) {
        const auto node=static_cast<Eigen::Index>(shape.active_nodes[k]);
        grid_curl.emplace_back(row,node,-shape.dy[k]);
        grid_curl.emplace_back(row,nodes+node,shape.dx[k]);
      }
    }
    csv.close();
  }
  const Eigen::MatrixXd curl_grid=sparse(grid_rows,ndof,grid_curl)*s.eigenvectors;
  if(!curl_grid.allFinite())throw std::runtime_error("Nonfinite exported physical curl");
  matrix(directory/"field_curl.mtx.gz",curl_grid);
  {
    auto json=output(directory/"run_metadata.json");
    const Eigen::MatrixXd gram=s.eigenvectors.transpose()*(matrices.mass*s.eigenvectors);
    const double orthogonality=(gram-Eigen::MatrixXd::Identity(ndof,ndof)).cwiseAbs().maxCoeff();
    const double symmetry_k=(Eigen::MatrixXd(matrices.stiffness)-Eigen::MatrixXd(matrices.stiffness.transpose())).norm();
    const double symmetry_m=(Eigen::MatrixXd(matrices.mass)-Eigen::MatrixXd(matrices.mass.transpose())).norm();
    const char* threads=std::getenv("OMP_NUM_THREADS");
    json<<"{\n  \"schema_version\": 1,\n  \"origin\": "<<json_string(options.origin)
      <<",\n  \"created_utc\": "<<json_string(utc_now())
      <<",\n  \"source_commit\": "<<json_string(AL_MALEH_GIT_COMMIT)
      <<",\n  \"geometry\": {\"xmin\": "<<r.xmin<<", \"xmax\": "<<r.xmax
      <<", \"ymin\": "<<r.ymin<<", \"ymax\": "<<r.ymax
      <<", \"a\": "<<a<<", \"b\": "<<b<<", \"area\": "<<a*b<<"},\n"
      <<"  \"cloud\": {\"nx\": "<<options.nx<<", \"ny\": "<<options.ny<<", \"nodes\": "<<nodes
      <<", \"degree\": "<<mls.polynomial_degree()<<", \"support_factor\": "<<options.support_factor
      <<", \"support_radius\": "<<mls.nodes().front().support_radius
      <<", \"basis\": \"polynomial shifted to query point, scaled by max support radius\", \"weight\": "
      <<json_string(mls_weight_name(mls.weight()))<<"},\n"
      <<"  \"quadrature\": {\"cells_x\": "<<rule.cells_x<<", \"cells_y\": "<<rule.cells_y
      <<", \"gauss_order\": "<<rule.gauss_order<<", \"points\": "<<matrices.quadrature_points
      <<", \"boundary_gauss_order\": 3, \"max_moment_condition\": "<<matrices.max_moment_condition<<"},\n"
      <<"  \"spectrum\": {\"dofs\": "<<ndof<<", \"modes\": "<<s.eigenvalues.size()
      <<", \"zero_cutoff\": "<<s.zero_cutoff<<", \"numerical_kernel_dimension\": "<<s.numerical_kernel_dimension
      <<", \"mass_condition\": "<<s.mass_condition<<", \"min_mass_eigenvalue\": "<<s.min_mass_eigenvalue
      <<", \"max_relative_residual\": "<<*std::max_element(s.residuals.begin(),s.residuals.end())
      <<", \"max_abs_M_orthogonality_error\": "<<orthogonality
      <<", \"K_symmetry_frobenius\": "<<symmetry_k<<", \"M_symmetry_frobenius\": "<<symmetry_m<<"},\n"
      <<"  \"solver\": {\"library\": \"Eigen\", \"method\": \"GeneralizedSelfAdjointEigenSolver\", "
        "\"arithmetic\": \"IEEE 754 binary64\", \"normalization\": \"original columns V; V^T M V = I; no rescaling\", "
        "\"residual_definition\": \"||Kv-lambda Mv||_2 / ((||K||_F+|lambda| ||M||_F)||v||_2)\"},\n"
      <<"  \"conventions\": {\"dof_order\": \"Hx_0..Hx_(N-1),Hy_0..Hy_(N-1)\", "
        "\"csv_index_base\": 0, \"matrix_market_index_base\": 1, \"matrix_market_column\": \"eigenvalue index + 1\", "
        "\"kc\": \"sqrt(max(0,lambda)); signed lambda retained\", \"text_significant_digits\": 17, "
        "\"near_zero\": \"lambda <= zero_cutoff; all columns retained\", "
        "\"divergence_ratio\": \"||div H||/(sqrt(lambda)||H||); blank for near_zero\", "
        "\"normal_ratio\": \"sqrt(integral_boundary((H.n)^2) * area/perimeter / field_energy)\"},\n"
      <<"  \"physics\": {\"mu\": 1, \"epsilon\": 1, \"coordinate_scale\": \"domain dimensions in geometry; no automatic rescaling; lambda has inverse length squared units\", "
        "\"boundary\": \"natural curl(H)=0; gradient kernel unfiltered\", "
        "\"classification\": \"experimental candidates; eigenvalue proximity alone does not establish physical identity\"},\n"
      <<"  \"common_grid\": {\"nx\": "<<options.grid_nx<<", \"ny\": "<<options.grid_ny
      <<", \"rows\": "<<grid_rows<<", \"ordering\": \"j-major, i-fast\", \"location\": \"cell centers\", "
        "\"field\": \"curl(H)=dHy/dx-dHx/dy; every spectral column\", \"weights\": \"grid.csv weight column\"},\n"
      <<"  \"environment\": {\"compiler\": "<<json_string(AL_MALEH_COMPILER)
      <<", \"compiler_version\": "<<json_string(__VERSION__)<<", \"build_type\": "<<json_string(AL_MALEH_BUILD_TYPE)
      <<", \"cplusplus\": "<<__cplusplus<<", \"Eigen\": "
      <<json_string(std::to_string(EIGEN_WORLD_VERSION)+"."+std::to_string(EIGEN_MAJOR_VERSION)+"."+std::to_string(EIGEN_MINOR_VERSION))
      <<", \"zlib\": "<<json_string(zlibVersion())<<", \"OMP_NUM_THREADS\": "<<json_string(threads?threads:"unset")<<"}\n}\n";
    json.close();
  }
}
} // namespace al_maleh
