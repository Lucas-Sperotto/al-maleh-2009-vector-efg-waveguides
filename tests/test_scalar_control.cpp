// SPDX-License-Identifier: MIT
#include "al_maleh/scalar_control.hpp"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace al_maleh;
namespace {
void require(bool ok,const std::string& reason){
  if(!ok)throw std::runtime_error(reason);
}
void near(double got,double expected,double tol,const char* name){
  require(std::isfinite(got)&&std::abs(got-expected)<=tol,
          std::string(name)+" got "+std::to_string(got)+
          " instead of "+std::to_string(expected));
}
void test_patch(int degree){
  const Rectangle r{0,2,0,1};
  const auto nodes=make_regular_cloud(r,5,5,1.65);
  const MLSEvaluator mls(nodes,degree);
  const BackgroundQuadrature q{r,4,4,3};
  const auto sys=assemble_scalar_control(mls,q);
  const Eigen::Index n=nodes.size();
  require(sys.boundary_points==16,"Expected 16 boundary collocation nodes");
  require(sys.boundary.rows()==16&&sys.boundary.cols()==n,"PEC matrix size");
  require(sys.nullspace.rows()==n&&sys.nullspace.cols()>0,"Nullspace size");
  require(sys.boundary_rank>0&&sys.boundary_rank<=16,"Boundary rank");
  require(sys.nullspace.cols()==n-sys.boundary_rank,"Reduced dimension");
  require(sys.quadrature_points==144,"Quadrature size");
  near((sys.stiffness-sys.stiffness.transpose()).norm(),0,1e-11,"Scalar stiffness symmetry");
  near((sys.mass-sys.mass.transpose()).norm(),0,1e-11,"Scalar mass symmetry");
  near((sys.boundary*sys.nullspace).norm(),0,1e-8,"Collocation nullspace");
  Eigen::VectorXd constant=Eigen::VectorXd::Ones(n);
  Eigen::VectorXd x(n);
  for(Eigen::Index i=0;i<n;++i)x[i]=nodes[static_cast<std::size_t>(i)].position.x;
  near(constant.dot(sys.stiffness*constant),0,1e-8,"Constant scalar gradient energy");
  near(constant.dot(sys.mass*constant),2,1e-8,"Scalar constant mass");
  near(x.dot(sys.stiffness*x),2,1e-8,"Linear scalar gradient energy");
  near(x.dot(sys.mass*x),8.0/3.0,1e-8,"Linear scalar x mass");
  const auto sol=solve_scalar_control(sys);
  require(sol.eigenvalues.size()==sys.nullspace.cols(),"Complete reduced spectrum");
  require(sol.coefficients.rows()==n,"Complete full coefficients");
  require(sol.eigenvalues.minCoeff()>0,"Scalar PEC positive spectrum");
  require(sol.mass_condition>0&&std::isfinite(sol.mass_condition),"Scalar mass condition");
  require(sol.mass_orthogonality_error<1e-7,"Scalar M orthogonality");
  require(sol.projected_residuals.size()==static_cast<std::size_t>(sol.eigenvalues.size()),"Missing projected residuals");
  for(std::size_t i=0;i<sol.projected_residuals.size();++i){
    require(sol.projected_residuals[i]<1e-8,"Projected residual too large");
    require(sol.collocation_residuals[i]<1e-8,"Exact boundary-node constraint lost");
  }
  auto metrics=evaluate_scalar_control_mode(mls,q,sol.coefficients.col(0),1,1);
  require(metrics.relative_shape_overlap>=0&&metrics.relative_shape_overlap<=1,"Overlap range");
  require(metrics.normal_flux_ratio>=0&&std::isfinite(metrics.normal_flux_ratio),"Normal trace range");
  require(metrics.boundary_trace_ratio>=0&&std::isfinite(metrics.boundary_trace_ratio),"PEC trace range");
  near(metrics.divergence_l2,0,1e-12,"Rotated gradient formal divergence");
}
void test_refinement(){
  const Rectangle r{0,2,0,1};
  const MLSEvaluator mls(make_regular_cloud(r,9,7,0.825),1);
  const BackgroundQuadrature q{r,8,6,3};
  const auto matrices=assemble_scalar_control(mls,q);
  const auto sol=solve_scalar_control(matrices);
  const double exact=1.25*std::acos(-1.)*std::acos(-1.);
  const auto m=evaluate_scalar_control_mode(mls,q,sol.coefficients.col(0),1,1);
  require(sol.eigenvalues.size()==matrices.nullspace.cols(),"Refined EVP dimensions");
  require(std::isfinite(sol.eigenvalues[0])&&sol.eigenvalues[0]>0,
          "Refined scalar first eigenvalue");
  require((matrices.boundary*sol.coefficients.col(0)).norm()<1e-7,
          "Refined scalar boundary collocation");
  // Record, but do not force, physical validation when studying
  // scalar collocation versus the unconstrained vector formulation.
  std::cout<<"M6c 9x7 lambda1="<<sol.eigenvalues[0]
           <<" relative_TM11="<<std::abs(sol.eigenvalues[0]-exact)/exact
           <<" TM11_overlap="<<m.relative_shape_overlap
           <<" PEC_trace="<<m.boundary_trace_ratio
           <<" normal_flux="<<m.normal_flux_ratio<<"\n";
}
} // namespace
int main(){
  try{
    test_patch(1);
    test_patch(2);
    test_refinement();
    std::cout<<"PASS: scalar MLS weak forms, PEC collocation rank and nullspace, full spectrum, algebraic residuals and field reconstruction\n";
  }catch(const std::exception& e){
    std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;
  }
}
