// SPDX-License-Identifier: MIT
// No change to the baseline curl-curl eigenproblem: diagnostic operators only.
#include "al_maleh/constraints.hpp"
#include "al_maleh/spectral.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace al_maleh;
int main(int argc,char** argv){
  try{
    if(argc!=6)
      throw std::invalid_argument("Usage: constraint_subspace nx ny index_1 index_2 output_dir");
    const int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
    const int i=std::stoi(argv[3]),j=std::stoi(argv[4]);
    const std::filesystem::path out=argv[5];
    if(nx<3||ny<3||nx>17||ny>17||i<0||j<=i||2*nx*ny>600||j>=2*nx*ny)
      throw std::invalid_argument("Invalid EFG nodal shape or distinct eigen-indices");
    const Rectangle r{0,2,0,1};
    const auto radius=3.3*std::max(2.0/(nx-1),1.0/(ny-1));
    const MLSEvaluator mls(make_regular_cloud(r,nx,ny,radius),1);
    const BackgroundQuadrature quadrature{r,nx-1,ny-1,3};
    const auto assembly=assemble_cartesian(mls,quadrature);
    const auto constraint=assemble_constraints(mls,quadrature,3);
    const auto spectrum=solve_generalized_symmetric(assembly.stiffness,assembly.mass);
    Eigen::MatrixXd V(2*nx*ny,2);
    V.col(0)=spectrum.eigenvectors.col(i);
    V.col(1)=spectrum.eigenvectors.col(j);
    const double pi=std::acos(-1.0);
    const double lambda_ref=5*pi*pi;
    const auto metrics=evaluate_subspace_constraints(
       constraint,assembly.mass,V,lambda_ref,r);

    // Verify basis invariance on the actual computed doublet, not merely
    // in a synthetic patch-test: rotations may mix nonidentical eigenvalues.
    Eigen::Matrix2d rotation;
    const double theta=0.713;
    rotation<<std::cos(theta),-std::sin(theta),
              std::sin(theta), std::cos(theta);
    const auto rotated=evaluate_subspace_constraints(
        constraint,assembly.mass,(V*rotation).eval(),lambda_ref,r);
    auto agrees=[](double x,double y){return std::abs(x-y)<1e-8;};
    if(!agrees(metrics.divergence_ratio.minimum,rotated.divergence_ratio.minimum)||
       !agrees(metrics.divergence_ratio.maximum,rotated.divergence_ratio.maximum)||
       !agrees(metrics.normal_ratio.minimum,rotated.normal_ratio.minimum)||
       !agrees(metrics.normal_ratio.maximum,rotated.normal_ratio.maximum))
      throw std::runtime_error("Measured doublet constraint extrema change under basis rotation");
    // A known physical direction is a rigorous lower bound on the
    // supremum of the constraint over the entire subspace.
    for(int k=0;k<2;++k){
      const Eigen::VectorXd v=V.col(k);
      const double norm=v.dot(assembly.mass*v);
      const double d=std::sqrt(std::max(0.0,v.dot(constraint.divergence*v))
                               /(lambda_ref*norm));
      if(d>metrics.divergence_ratio.maximum+1e-10)
        throw std::runtime_error("Invariant subspace maximum violates Rayleigh lower bound");
    }
    std::filesystem::create_directories(out);
    std::ofstream csv(out/"subspace_constraints.csv");
    if(!csv)throw std::runtime_error("Cannot write subspace_constraints.csv");
    csv<<std::setprecision(17)
       <<"nx,ny,dofs,index_1,index_2,lambda_1,lambda_2,lambda_reference,"
         "divergence_ratio_min,divergence_ratio_max,normal_ratio_min,normal_ratio_max,"
         "divergence_worst_coeff_1,divergence_worst_coeff_2,"
         "normal_worst_coeff_1,normal_worst_coeff_2,passes_divergence_threshold\n"
       <<nx<<','<<ny<<','<<2*nx*ny<<','<<i<<','<<j<<','
       <<spectrum.eigenvalues[i]<<','<<spectrum.eigenvalues[j]<<','
       <<lambda_ref<<','
       <<metrics.divergence_ratio.minimum<<','<<metrics.divergence_ratio.maximum<<','
       <<metrics.normal_ratio.minimum<<','<<metrics.normal_ratio.maximum<<','
       <<metrics.divergence_ratio.maximizing_coefficients[0]<<','
       <<metrics.divergence_ratio.maximizing_coefficients[1]<<','
       <<metrics.normal_ratio.maximizing_coefficients[0]<<','
       <<metrics.normal_ratio.maximizing_coefficients[1]<<','
       <<(metrics.divergence_ratio.maximum<=0.25?"true":"false")<<'\n';

    std::ofstream low(out/"positive_below_tm11.csv");
    if(!low)throw std::runtime_error("Cannot write positive_below_tm11.csv");
    const double tm11=1.25*pi*pi;
    low<<std::setprecision(17)
       <<"index,lambda,divergence_ratio,normal_ratio,relative_residual\n";
    int low_count=0;
    for(Eigen::Index n=0;n<spectrum.eigenvalues.size();++n){
      const double lambda=spectrum.eigenvalues[n];
      if(lambda<=spectrum.zero_cutoff||lambda>=tm11)continue;
      const auto v=spectrum.eigenvectors.col(n).eval();
      const double norm=v.dot(assembly.mass*v);
      const double div2=v.dot(constraint.divergence*v);
      const double normal2=v.dot(constraint.normal_boundary*v);
      const double div=std::sqrt(std::max(0.0,div2)/(lambda*norm));
      const double normal=std::sqrt(std::max(0.0,normal2)/(3.0*norm));
      low<<n<<','<<lambda<<','<<div<<','<<normal<<','
         <<spectrum.residuals[static_cast<std::size_t>(n)]<<'\n';
      ++low_count;
    }
    std::ofstream meta(out/"run_metadata.txt");
    meta<<std::setprecision(17)
        <<"status=diagnostic_only_no_penalty\n"
        <<"boundary_constraint=normal_flux_squared\n"
        <<"divergence_constraint=L2_divergence_squared\n"
        <<"gauss_volume=3\ngauss_segmented_boundary=3\n"
        <<"support_factor=3.3\n"
        <<"nx="<<nx<<"\nny="<<ny<<"\n"
        <<"index_1="<<i<<"\nindex_2="<<j<<"\n"
        <<"dofs="<<2*nx*ny<<"\n"
        <<"numerical_kernel="<<spectrum.numerical_kernel_dimension<<"\n"
        <<"positive_below_tm11="<<low_count<<"\n"
        <<"mass_condition="<<spectrum.mass_condition<<"\n";
    std::cout<<std::setprecision(9)
       <<"M6b "<<nx<<'x'<<ny<<" DOFs="<<2*nx*ny
       <<" invariant_div_range=["<<metrics.divergence_ratio.minimum
       <<','<<metrics.divergence_ratio.maximum<<']'
       <<" invariant_normal_range=["<<metrics.normal_ratio.minimum
       <<','<<metrics.normal_ratio.maximum<<']'
       <<" low_positive="<<low_count<<" output="<<out<<"\n";
  }catch(const std::exception& e){
    std::cerr<<"M6b ERROR: "<<e.what()<<"\n";return 1;
  }
}
