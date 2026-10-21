// SPDX-License-Identifier: MIT
#include "al_maleh/constraints.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
using namespace al_maleh;
namespace {
void check_near(double x,double y,double tol,const char* label) {
  if(!std::isfinite(x)||std::abs(x-y)>tol)
    throw std::runtime_error(std::string(label)+" unexpected value "+
       std::to_string(x)+" expected "+std::to_string(y));
}
void expect_invalid(const std::function<void()>& f) {
  try{ f(); }catch(const std::invalid_argument&){return;}
  throw std::runtime_error("Expected invalid_argument");
}
double energy(const Eigen::SparseMatrix<double>& a,const Eigen::VectorXd& v){
  return v.dot(a*v);
}
void run(int degree,int boundary_order) {
  const Rectangle r{0,2,0,1};
  const auto nodes=make_regular_cloud(r,5,5,3.0);
  const MLSEvaluator mls(nodes,degree);
  const BackgroundQuadrature q{r,4,4,3};
  const auto a=assemble_cartesian(mls,q);
  const auto c=assemble_constraints(mls,q,boundary_order);
  const Eigen::Index N=nodes.size();
  if(c.boundary_points!=static_cast<std::size_t>(2*(q.cells_x+q.cells_y)*boundary_order))
    throw std::runtime_error("Wrong boundary quadrature count");
  const Eigen::MatrixXd D=Eigen::MatrixXd(c.divergence);
  const Eigen::MatrixXd B=Eigen::MatrixXd(c.normal_boundary);
  check_near((D-D.transpose()).norm(),0,1e-10,"divergence symmetry");
  check_near((B-B.transpose()).norm(),0,1e-10,"boundary symmetry");
  Eigen::VectorXd cx=Eigen::VectorXd::Zero(2*N),
                  cy=Eigen::VectorXd::Zero(2*N),
                  rotation=Eigen::VectorXd::Zero(2*N),
                  grad=Eigen::VectorXd::Zero(2*N);
  for(Eigen::Index i=0;i<N;++i){
    const Point2 p=nodes[static_cast<std::size_t>(i)].position;
    cx[i]=1;
    cy[N+i]=1;
    rotation[i]=-p.y;rotation[N+i]=p.x;
    grad[i]=p.x;grad[N+i]=p.y;
  }
  // Divergence of constant and rotation fields is zero.
  // Divergence of grad(.5(x²+y²))=2 -> integral div²=2*4=8.
  check_near(energy(c.divergence,cx),0,1e-9,"constant div");
  check_near(energy(c.divergence,cy),0,1e-9,"constant y div");
  check_near(energy(c.divergence,rotation),0,1e-8,"rotation div");
  check_near(energy(c.divergence,grad),8,1e-8,"linear gradient div");
  // Integrated squared boundary normal components for rectangle 2x1.
  check_near(energy(c.normal_boundary,cx),2,1e-8,"x constant boundary");
  check_near(energy(c.normal_boundary,cy),4,1e-8,"y constant boundary");
  check_near(energy(c.normal_boundary,rotation),6,1e-8,"rotation normal");
  check_near(energy(c.normal_boundary,grad),6,1e-8,"gradient normal");
  Eigen::MatrixXd V(2*N,2);
  V.col(0)=cx;V.col(1)=grad;
  const auto first=evaluate_subspace_constraints(c,a.mass,V,5.0,r);
  const double theta=0.713;
  Eigen::Matrix2d R;
  R<<std::cos(theta),-std::sin(theta),std::sin(theta),std::cos(theta);
  const auto rotated=evaluate_subspace_constraints(c,a.mass,(V*R).eval(),5.0,r);
  for(const auto* m:{&first,&rotated}){
    if(m->dimension!=2)throw std::runtime_error("Bad subspace dimension");
    check_near(m->divergence_ratio.minimum,0,1e-7,"min divergence");
    if(m->divergence_ratio.maximum<0.4)throw std::runtime_error("Max divergence not physical");
    if(m->normal_ratio.minimum<0||m->normal_ratio.maximum<0)
      throw std::runtime_error("Negative constraint metric");
  }
  // Near a true zero singular constraint, sqrt(roundoff) is not smooth.\n  // Compare the squared extremal ratios (generalized eigenvalues) instead.\n  check_near(first.divergence_ratio.minimum*first.divergence_ratio.minimum,\n             rotated.divergence_ratio.minimum*rotated.divergence_ratio.minimum,\n             1e-10,"squared div minimum invariance");
  check_near(first.divergence_ratio.maximum,rotated.divergence_ratio.maximum,1e-10,"div maximum invariance");
  check_near(first.normal_ratio.minimum,rotated.normal_ratio.minimum,1e-10,"wall minimum invariance");
  check_near(first.normal_ratio.maximum,rotated.normal_ratio.maximum,1e-10,"wall maximum invariance");
  // Both extrema must reproduce the quotient of their own maximizing vector.
  const Eigen::VectorXd x=V*first.divergence_ratio.maximizing_coefficients;
  const double quotient=std::sqrt(energy(c.divergence,x)/(5.0*energy(a.mass,x)));
  check_near(quotient,first.divergence_ratio.maximum,1e-9,"maximizer quotient");
  expect_invalid([&]{ (void)evaluate_subspace_constraints(c,a.mass,
                 (V*Eigen::Matrix2d::Zero()).eval(),5,r); });
  expect_invalid([&]{ (void)evaluate_subspace_constraints(c,a.mass,V,0,r); });
  expect_invalid([&]{ (void)assemble_constraints(mls,q,5); });
}
} // namespace
int main() {
  try {
    for(int degree:{1,2}) for(int order:{2,3,4})run(degree,order);
    std::cout<<"PASS: D/B patch energies, segmented boundary quadrature, SPD Gram, basis-invariant extrema\n";
  }catch(const std::exception& e){
    std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;
  }
}
