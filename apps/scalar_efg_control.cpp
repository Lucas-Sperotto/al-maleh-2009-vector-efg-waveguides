// SPDX-License-Identifier: MIT
#include "al_maleh/scalar_control.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {
void save_mtx(const std::filesystem::path& path,const Eigen::MatrixXd& A){
  std::ofstream out(path);
  if(!out)throw std::runtime_error("Cannot open Matrix Market export "+path.string());
  out<<"%%MatrixMarket matrix array real general\n"
     <<A.rows()<<" "<<A.cols()<<"\n"<<std::setprecision(17);
  for(Eigen::Index col=0;col<A.cols();++col)
    for(Eigen::Index row=0;row<A.rows();++row)
      out<<A(row,col)<<"\n";
  if(!out)throw std::runtime_error("Matrix Market export write failed");
}
}
int main(int argc,char** argv){
  try{
    if(argc!=4)throw std::invalid_argument("Usage: scalar_efg_control nx ny output_dir");
    const int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
    const std::filesystem::path out=argv[3];
    if(nx<4||ny<4||nx>17||ny>17||nx*ny>300)
      throw std::invalid_argument("Scalar control grid must be 4..17 and <=300 nodes");
    const Rectangle domain{0,2,0,1};
    const double support_factor=3.3;
    const double radius=support_factor*std::max(2.0/(nx-1),1.0/(ny-1));
    const MLSEvaluator mls(make_regular_cloud(domain,nx,ny,radius),1);
    const BackgroundQuadrature q{domain,nx-1,ny-1,3};
    const auto matrices=assemble_scalar_control(mls,q);
    const auto solution=solve_scalar_control(matrices);
    std::filesystem::create_directories(out);
    save_mtx(out/"stiffness_matrix.mtx",matrices.stiffness);
    save_mtx(out/"mass_matrix.mtx",matrices.mass);
    save_mtx(out/"boundary_collocation.mtx",matrices.boundary);
    save_mtx(out/"boundary_nullspace.mtx",matrices.nullspace);
    save_mtx(out/"eigenvectors.mtx",solution.coefficients);
    std::ofstream vals(out/"eigenvalues.csv");
    if(!vals)throw std::runtime_error("Cannot open scalar eigenvalue CSV");
    vals<<std::setprecision(17)
        <<"index,lambda,kc,projected_relative_residual,pec_collocation_relative_residual\n";
    for(Eigen::Index i=0;i<solution.eigenvalues.size();++i)
      vals<<i<<','<<solution.eigenvalues[i]<<','
          <<std::sqrt(solution.eigenvalues[i])<<','
          <<solution.projected_residuals[static_cast<std::size_t>(i)]<<','
          <<solution.collocation_residuals[static_cast<std::size_t>(i)]<<'\n';

    const double pi=std::acos(-1.0);
    const std::vector<std::pair<int,int>> modes{
        {1,1},{2,1},{3,1},{1,2},{2,2},{4,1}};
    std::ofstream comparison(out/"tm_reference.csv");
    if(!comparison)throw std::runtime_error("Cannot open scalar TM reference");
    comparison<<std::setprecision(17)
       <<"m,n,lambda_analytic,index_nearest,lambda_nearest,relative_lambda_error,"
         "scalar_overlap,boundary_trace_ratio,normal_flux_ratio,"
         "formal_divergence_l2,interpretation\n";
    for(const auto& mode:modes){
      const double exact=pi*pi*(mode.first*mode.first/4.0+
                                mode.second*mode.second);
      Eigen::Index best=-1;
      double error=std::numeric_limits<double>::infinity();
      for(Eigen::Index i=0;i<solution.eigenvalues.size();++i){
        const double e=std::abs(solution.eigenvalues[i]-exact)/exact;
        if(e<error){error=e;best=i;}
      }
      const auto f=evaluate_scalar_control_mode(
          mls,q,solution.coefficients.col(best),mode.first,mode.second);
      comparison<<mode.first<<','<<mode.second<<','<<exact<<','
                <<best<<','<<solution.eigenvalues[best]<<','<<error<<','
                <<f.relative_shape_overlap<<','<<f.boundary_trace_ratio<<','
                <<f.normal_flux_ratio<<','<<f.divergence_l2<<','
                <<(exact==5*pi*pi?"degenerate_needs_joint_match":
                   "exploratory_nearest_eigenvalue")<<'\n';
    }
    std::ofstream meta(out/"run_metadata.txt");
    if(!meta)throw std::runtime_error("Cannot open run metadata");
    meta<<std::setprecision(17)
        <<"model=independent_scalar_EFG_TM_control_not_vector_EFG\n"
        <<"boundary_condition=PEC_Ez_dirichlet_collocation_only\n"
        <<"boundary_between_collocation_points_not_exactly_zero\n"
        <<"operator=-Laplacian_scalar_weak_form\n"
        <<"normalization=full_coefficients_mass_orthonormal\n"
        <<"reduced_equation=Qt_S_Q_y_lambda_Qt_M_Q_y\n"
        <<"residual=projected_onto_constraint_nullspace\n"
        <<"domain_a=2\ndomain_b=1\n"
        <<"nx="<<nx<<"\nny="<<ny<<"\n"
        <<"nodes="<<nx*ny<<"\n"
        <<"boundary_points="<<matrices.boundary_points<<"\n"
        <<"boundary_rank="<<matrices.boundary_rank<<"\n"
        <<"reduced_dofs="<<solution.eigenvalues.size()<<"\n"
        <<"boundary_constraint_norm="<<matrices.boundary_constraint_residual<<"\n"
        <<"mass_condition="<<solution.mass_condition<<"\n"
        <<"mass_orthogonality_max_error="<<solution.mass_orthogonality_error<<"\n"
        <<"support_factor="<<support_factor<<"\n"
        <<"support_radius="<<radius<<"\n"
        <<"quadrature_order="<<q.gauss_order<<"\n"
        <<"quadrature_points="<<matrices.quadrature_points<<"\n"
        <<"highest_mls_moment_condition="<<matrices.max_moment_condition<<"\n";
    std::cout<<std::setprecision(12)
       <<"M6c scalar EFG "<<nx<<'x'<<ny<<" nodes="<<nx*ny
       <<" boundary_rank="<<matrices.boundary_rank
       <<" reduced_dofs="<<solution.eigenvalues.size()
       <<" lambda_1="<<solution.eigenvalues[0]
       <<" mass_condition="<<solution.mass_condition
       <<" output="<<out<<"\n";
  }catch(const std::exception& e){
    std::cerr<<"M6c ERROR: "<<e.what()<<"\n";return 1;
  }
}
