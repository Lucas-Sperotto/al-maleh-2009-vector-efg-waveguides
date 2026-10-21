// SPDX-License-Identifier: MIT
#include "al_maleh/tm_validation.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace al_maleh;
int main(int argc,char** argv) {
 try {
  // tm_validation [nx ny degree support_factor gauss_order output_directory]
  if(argc!=7)throw std::invalid_argument(
    "Usage: tm_validation nx ny degree support_factor gauss_order output_dir");
  const int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
  const int degree=std::stoi(argv[3]);
  const double factor=std::stod(argv[4]);
  const int gauss=std::stoi(argv[5]);
  const std::filesystem::path dir=argv[6];
  if(nx<3||ny<3||nx>17||ny>17||(degree!=1&&degree!=2)||
     !std::isfinite(factor)||factor<=0||
     (gauss!=2&&gauss!=3&&gauss!=4))
    throw std::invalid_argument("Invalid M4 run parameters");
  const Rectangle r{0,2,0,1};
  const double hx=2.0/(nx-1),hy=1.0/(ny-1);
  const double radius=factor*std::max(hx,hy);
  const MLSEvaluator mls(make_regular_cloud(r,nx,ny,radius),degree);
  const BackgroundQuadrature quadrature{r,nx-1,ny-1,gauss};
  const auto matrices=assemble_cartesian(mls,quadrature);
  const auto s=solve_generalized_symmetric(matrices.stiffness,matrices.mass);
  const auto reference=validate_tm_modes(mls,quadrature,matrices,s,2,2);

  std::filesystem::create_directories(dir);
  std::ofstream csv(dir/"tm_matches.csv");
  if(!csv)throw std::runtime_error("Cannot create matches file");
  csv<<std::setprecision(17)
     <<"m,n,lambda_exact,kc_exact,nearest_index,nearest_lambda,nearest_relative_error,"
       "nearest_curl_overlap,nearest_cluster_overlap,nearest_cluster_size,"
       "strongest_index,strongest_lambda,strongest_overlap,strongest_relative_error,"
       "divergence_ratio,normal_boundary_ratio,algebraic_residual,status\n";
  for(const auto& t:reference) {
    csv<<t.reference.m<<','<<t.reference.n<<','
       <<t.reference.lambda<<','<<t.reference.kc<<','
       <<t.nearest_index<<','<<t.nearest_lambda<<','
       <<t.nearest_relative_error<<','<<t.nearest_overlap<<','
       <<t.nearest_cluster_overlap<<','<<t.nearest_cluster_size<<','
       <<t.strongest_index<<','<<t.strongest_lambda<<','
       <<t.strongest_overlap<<','<<t.strongest_relative_error<<','
       <<t.divergence_ratio<<','<<t.normal_ratio<<','
       <<t.algebraic_residual<<','
       <<(t.supported_candidate?"shape_supported_candidate":"unconfirmed_candidate")
       <<'\n';
  }
  std::ofstream spectrum(dir/"candidate_spectrum.csv");
  if(!spectrum)throw std::runtime_error("Cannot create spectrum file");
  spectrum<<std::setprecision(17)
          <<"index,lambda,residual,class,divergence_ratio,normal_boundary_ratio\n";
  int candidates=0;
  for(Eigen::Index i=0;i<s.eigenvalues.size();++i) {
    if(s.eigenvalues[i]<=s.zero_cutoff)continue;
    if(candidates>=16)break;
    const auto d=modal_diagnostics(
      mls,quadrature,s.eigenvectors.col(i),s.eigenvalues[i]);
    spectrum<<i<<','<<s.eigenvalues[i]<<','
      <<s.residuals[static_cast<std::size_t>(i)]
      <<",positive_candidate,"<<d.divergence_ratio<<','<<d.normal_ratio<<'\n';
    ++candidates;
  }
  std::ofstream meta(dir/"run_metadata.txt");
  if(!meta)throw std::runtime_error("Cannot create metadata");
  meta<<std::setprecision(17)
      <<"model=M4_cartesian_EFG_validation\n"
      <<"domain=rectangle_2x1\n"
      <<"status=experimental_unverified_physics\n"
      <<"nx="<<nx<<"\nny="<<ny<<"\ndegree="<<degree
      <<"\nsupport_factor="<<factor<<"\nsupport_radius="<<radius
      <<"\ngauss_order="<<gauss<<"\n"
      <<"dofs="<<s.eigenvalues.size()
      <<"\nquadrature_points="<<matrices.quadrature_points
      <<"\nnumerical_kernel="<<s.numerical_kernel_dimension
      <<"\nzero_cutoff="<<s.zero_cutoff
      <<"\nmass_condition="<<s.mass_condition
      <<"\nmass_min_eigenvalue="<<s.min_mass_eigenvalue
      <<"\nmls_max_condition="<<matrices.max_moment_condition
      <<"\ncandidates_checked="<<candidates<<'\n';
  const auto& tm11=reference.front();
  std::cout<<std::setprecision(10)
      <<"M4 nx="<<nx<<" ny="<<ny<<" degree="<<degree
      <<" support="<<factor<<" gauss="<<gauss
      <<" nullity="<<s.numerical_kernel_dimension
      <<" mass_condition="<<s.mass_condition
      <<" TM11_nearest_relative_error="<<tm11.nearest_relative_error
      <<" TM11_shape_overlap="<<tm11.nearest_cluster_overlap
      <<" TM11_status="<<(tm11.supported_candidate?"shape_supported":"unconfirmed")
      <<" output="<<dir.string()<<'\n';
 }catch(const std::exception& e) {
  std::cerr<<"M4 ERROR: "<<e.what()<<'\n';
  return 1;
 }
}
