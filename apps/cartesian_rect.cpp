// SPDX-License-Identifier: MIT
#include "al_maleh/cartesian.hpp"

#include <Eigen/Core>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    // Usage: cartesian_rect [nx ny [degree [support_factor]]]
    const int nx=argc>1?std::stoi(argv[1]):7;
    const int ny=argc>2?std::stoi(argv[2]):5;
    const int degree=argc>3?std::stoi(argv[3]):1;
    const double factor=argc>4?std::stod(argv[4]):4.0;
    if (argc>5 || nx<2 || ny<2 || !(factor>0.0))
      throw std::invalid_argument("Usage: cartesian_rect [nx ny [degree=1|2 [support_factor]]]");
    const al_maleh::Rectangle rect{0,2,0,1};
    const double dx=2.0/(nx-1),dy=1.0/(ny-1);
    const double radius=factor*std::max(dx,dy);
    al_maleh::MLSEvaluator mls(al_maleh::make_regular_cloud(rect,nx,ny,radius),degree);
    const al_maleh::BackgroundQuadrature quad{rect,nx-1,ny-1,3};
    const auto matrices=al_maleh::assemble_cartesian(mls,quad);
    const Eigen::Index n=static_cast<Eigen::Index>(mls.nodes().size());
    Eigen::VectorXd ex=Eigen::VectorXd::Zero(2*n);
    ex.head(n).setOnes();
    std::cout<<"model=cartesian_efg_m2_prototype\n"
             <<"domain=rectangular_2x1\n"
             <<"boundary=natural_curlH_zero_unconstrained_gradient_kernel\n"
             <<"nodes="<<n<<" dofs="<<2*n<<" degree="<<degree<<" support="<<radius<<"\n"
             <<"quadrature_points="<<matrices.quadrature_points
             <<" max_active_nodes="<<matrices.max_active_nodes<<"\n"
             <<"max_moment_condition="<<matrices.max_moment_condition<<"\n"
             <<"K_nonzeros="<<matrices.stiffness.nonZeros()
             <<" M_nonzeros="<<matrices.mass.nonZeros()<<"\n"
             <<"constant_field_mass="<<ex.dot(matrices.mass*ex)
             <<" constant_field_curl_energy="<<ex.dot(matrices.stiffness*ex)<<"\n";
    std::cout<<"note=no_generalized_eigen_solver_yet\n";
  } catch(const std::exception& e) {
    std::cerr<<"ERROR: "<<e.what()<<"\n";
    return 1;
  }
}
