// SPDX-License-Identifier: MIT
// M6d: diagnostic only; K/M and the eigenvalues are unchanged.
#include "al_maleh/gradient_audit.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {
void write_matrix(const std::filesystem::path& path,const Eigen::MatrixXd& mat){
  std::ofstream out(path);
  if(!out)throw std::runtime_error("Cannot create Matrix Market file");
  out<<"%%MatrixMarket matrix array real general\n"
     <<mat.rows()<<" "<<mat.cols()<<"\n"<<std::setprecision(17);
  for(Eigen::Index j=0;j<mat.cols();++j)
    for(Eigen::Index i=0;i<mat.rows();++i)
      out<<mat(i,j)<<"\n";
  if(!out)throw std::runtime_error("Failed writing Matrix Market file");
}
std::vector<std::pair<std::string,int>> named_modes(int nx,int ny){
  if(nx==9&&ny==7)return {{"TM11",54},{"TM21",60},{"TM31",63},{"TM12",66},
                          {"doublet_A",70},{"doublet_B",71}};
  if(nx==13&&ny==11)return {{"TM11",100},{"TM21",108},{"TM31",117},{"TM12",123},
                            {"doublet_A",127},{"doublet_B",128}};
  if(nx==17&&ny==15)return {{"TM11",153},{"TM21",163},{"TM31",178},{"TM12",185},
                            {"doublet_A",192},{"doublet_B",193}};
  return {};
}
}
int main(int argc,char** argv) {
  try{
    if(argc!=4&&argc!=5)
      throw std::invalid_argument(
          "Usage: gradient_kernel_audit nx ny output_dir [rank_relative_cutoff]");
    const int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
    const std::filesystem::path out=argv[3];
    if(nx<4||ny<4||nx>17||ny>17||2*nx*ny>600)
      throw std::invalid_argument("Invalid M6d grid size");
    GradientAuditOptions options;
    if(argc==5)options.relative_rank_cutoff=std::stod(argv[4]);
    const Rectangle r{0,2,0,1};
    const double radius=3.3*std::max(2.0/(nx-1),1.0/(ny-1));
    const MLSEvaluator mls(make_regular_cloud(r,nx,ny,radius),1);
    const BackgroundQuadrature q{r,nx-1,ny-1,3};
    const auto assembly=assemble_cartesian(mls,q);
    const auto spectrum=solve_generalized_symmetric(assembly.stiffness,assembly.mass);
    const auto pairing=assemble_gradient_pairing(mls,q);
    const auto projection=build_gradient_projection(assembly.mass,pairing,options);
    const double pi=std::acos(-1.0),first_tm=1.25*pi*pi;
    std::filesystem::create_directories(out);
    const auto original_modes=named_modes(nx,ny);

    std::ofstream modes(out/"eigenmode_gradient.csv");
    modes<<std::setprecision(17)
         <<"index,lambda,original_class,relative_evp_residual,"
           "gradient_energy_fraction,complement_energy_fraction,"
           "projected_gradient_curl_rayleigh,projected_gradient_norm\n";
    int low=0;
    double low_grad_sum=0,low_min=1,low_max=0;
    for(Eigen::Index j=0;j<spectrum.eigenvalues.size();++j){
      const auto m=evaluate_gradient_mode(projection,assembly.mass,
                            assembly.stiffness,spectrum.eigenvectors.col(j));
      const double lambda=spectrum.eigenvalues[j];
      const bool below=(lambda>spectrum.zero_cutoff&&lambda<first_tm);
      if(below) {
        ++low;
        low_grad_sum+=m.gradient_energy_fraction;
        low_min=std::min(low_min,m.gradient_energy_fraction);
        low_max=std::max(low_max,m.gradient_energy_fraction);
      }
      modes<<j<<','<<lambda<<','
          <<(lambda<=spectrum.zero_cutoff?"near_zero":
             below?"positive_below_TM11":"positive_other")<<','
          <<spectrum.residuals[static_cast<std::size_t>(j)]<<','
          <<m.gradient_energy_fraction<<','
          <<m.complement_energy_fraction<<','
          <<m.projected_gradient_curl_rayleigh<<','
          <<m.projected_vector_norm<<"\n";
    }
    if(!modes)throw std::runtime_error("Cannot write complete M6d spectrum audit");

    std::ofstream selected(out/"selected_modes.csv");
    selected<<std::setprecision(17)
            <<"label,index,lambda,gradient_energy_fraction,"
              "complement_energy_fraction,gradient_curl_rayleigh\n";
    for(const auto& [name,index]:original_modes){
      if(index>=spectrum.eigenvalues.size())
        throw std::runtime_error("Archived selected index outside spectrum");
      const auto m=evaluate_gradient_mode(projection,assembly.mass,
                  assembly.stiffness,spectrum.eigenvectors.col(index));
      selected<<name<<','<<index<<','<<spectrum.eigenvalues[index]<<','
              <<m.gradient_energy_fraction<<','
              <<m.complement_energy_fraction<<','
              <<m.projected_gradient_curl_rayleigh<<"\n";
    }

    const auto n=static_cast<Eigen::Index>(mls.nodes().size());
    Eigen::VectorXd x(n),y(n),quadratic(n),reference(n);
    for(Eigen::Index j=0;j<n;++j) {
      const Point2 p=mls.nodes()[static_cast<std::size_t>(j)].position;
      x[j]=p.x;y[j]=p.y;
      quadratic[j]=p.x*p.x+p.y*p.y;
      reference[j]=std::sin(pi*p.x/2)*std::sin(pi*p.y);
    }
    std::ofstream probe(out/"gradient_fit_probes.csv");
    probe<<std::setprecision(17)<<"scalar_probe,relative_l2_best_vector_fit_error\n";
    for(const auto& [name,v]:std::vector<std::pair<std::string,Eigen::VectorXd>>{
       {"x",x},{"y",y},{"x2_plus_y2",quadratic},{"sin_pi_x_over_2_sin_pi_y",reference}}){
      probe<<name<<','
           <<scalar_gradient_fit_relative_error(pairing,projection,v)<<"\n";
    }

    write_matrix(out/"gradient_coupling.mtx",pairing.coupling);
    write_matrix(out/"scalar_gradient_energy.mtx",pairing.scalar_energy);
    write_matrix(out/"projected_scalar_gram.mtx",projection.projected_scalar_gram);
    write_matrix(out/"gradient_projector_basis.mtx",projection.basis);
    // Full eigenpairs were already archived in M6a; persist the current vector
    // mass and eigenbasis here as well for standalone projection verification.
    write_matrix(out/"vector_mass.mtx",Eigen::MatrixXd(assembly.mass));
    write_matrix(out/"eigenvectors.mtx",spectrum.eigenvectors);
    std::ofstream singular(out/"gradient_gram_eigenvalues.csv");
    singular<<std::setprecision(17)<<"index,sigma,kept\n";
    for(Eigen::Index i=0;i<projection.projected_scalar_eigenvalues.size();++i)
      singular<<i<<','<<projection.projected_scalar_eigenvalues[i]<<','
              <<(projection.projected_scalar_eigenvalues[i]>projection.rank_cutoff
                 ?"true":"false")<<"\n";

    std::ofstream metadata(out/"metadata.txt");
    metadata<<std::setprecision(17)
            <<"model=M6d_diagnostic_only_unmodified_cartesian_EFG\n"
            <<"gradient_pairing=int_Ni_dot_grad_phi_j\n"
            <<"scalar_gradient_energy=int_grad_phi_i_dot_grad_phi_j\n"
            <<"operator_projection=M_inverse_H\n"
            <<"projector=Q_Qtrans_M\n"
            <<"not_verified_de_rham_compatible=true\n"
            <<"nx="<<nx<<"\nny="<<ny<<"\nnodes="<<n
            <<"\nvector_dofs="<<2*n
            <<"\nzero_cutoff="<<spectrum.zero_cutoff
            <<"\nnumerical_kernel="<<spectrum.numerical_kernel_dimension
            <<"\npositive_below_TM11="<<low
            <<"\nprojected_gradient_rank="<<projection.numerical_rank
            <<"\ngradient_rank_relative_cutoff="<<options.relative_rank_cutoff
            <<"\ngradient_rank_absolute_cutoff="<<projection.rank_cutoff
            <<"\nvector_mass_condition="<<projection.mass_condition
            <<"\nmax_M_basis_orthogonality_error="
                <<projection.maximum_m_orthogonality_error
            <<"\ngradient_energy_matrix_relative_discrepancy="
                <<projection.relative_gradient_reconstruction_frobenius
            <<"\nlow_positive_gradient_fraction_mean="
                <<(low?low_grad_sum/low:0)
            <<"\nlow_positive_gradient_fraction_min="<<(low?low_min:0)
            <<"\nlow_positive_gradient_fraction_max="<<low_max
            <<"\nvolume_quadrature_points="<<pairing.quadrature_points
            <<"\nsupport_factor=3.3\nquadrature_order=3\n";
    if(!metadata||!probe||!selected||!singular)
      throw std::runtime_error("M6d output stream failed");
    std::cout<<std::setprecision(12)
             <<"M6d "<<nx<<"x"<<ny<<" rank="<<projection.numerical_rank
             <<" low_positive="<<low<<" low_gradient_mean="
             <<(low?low_grad_sum/low:0)
             <<" mass_orthogonality="
             <<projection.maximum_m_orthogonality_error
             <<" output="<<out<<"\n";
  }catch(const std::exception& e) {
    std::cerr<<"M6d ERROR: "<<e.what()<<"\n";return 1;
  }
}
