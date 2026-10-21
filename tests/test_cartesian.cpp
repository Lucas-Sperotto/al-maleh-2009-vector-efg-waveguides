// SPDX-License-Identifier: MIT
#include "al_maleh/cartesian.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {

void near(double actual, double expected, double tol, const std::string& what) {
  if (!std::isfinite(actual) || std::abs(actual-expected) > tol)
    throw std::runtime_error(what + " expected=" + std::to_string(expected) +
                             " got=" + std::to_string(actual));
}

void reject(const std::function<void()>& fn, const std::string& label) {
  try { fn(); } catch (const std::exception&) { return; }
  throw std::runtime_error(label+" should have failed");
}

double quadratic(const Eigen::SparseMatrix<double>& mat,
                 const Eigen::VectorXd& x) {
  return x.dot(mat*x);
}

void test_quadrature() {
  BackgroundQuadrature rule{{0.0, 2.0, 0.0, 1.0}, 3, 2, 3};
  double volume=0, x=0, y=0, x2=0, xy=0, y2=0;
  int evaluations=0;
  for_each_quadrature_point(rule, [&](Point2 p, double w) {
    ++evaluations; volume+=w; x+=w*p.x; y+=w*p.y;
    x2+=w*p.x*p.x; xy+=w*p.x*p.y; y2+=w*p.y*p.y;
  });
  near(volume,2.0,2e-14,"integration area");
  near(x,2.0,2e-14,"integration x");
  near(y,1.0,2e-14,"integration y");
  near(x2,8.0/3.0,2e-14,"integration x²");
  near(y2,2.0/3.0,2e-14,"integration y²");
  near(xy,1.0,2e-14,"integration xy");
  if(evaluations!=54)throw std::runtime_error("Wrong quadrature size");

  reject([&]{ BackgroundQuadrature bad=rule; bad.gauss_order=5;
    for_each_quadrature_point(bad, [](Point2,double){}); },"gauss order");
  reject([&]{ BackgroundQuadrature bad=rule; bad.cells_y=0;
    for_each_quadrature_point(bad, [](Point2,double){}); },"zero cells");
  reject([&]{ BackgroundQuadrature bad=rule; bad.domain.xmax=0;
    for_each_quadrature_point(bad, [](Point2,double){}); },"invalid domain");
}

// Original Chapter 5 of Al Maleh (2009) uses 1x1 through 5x5 Gauss
// points and four background cells for each nodal rectangle (Fig. 9).
// Verify physical weights and Gauss exactness independently of MLS.
void test_author_ch5_quadrature() {
  const Rectangle rect{0,2,0,1};
  for(int order=1;order<=5;++order) {
    BackgroundQuadrature rule{rect,2*(5-1),2*(3-1),order};
    int count=0;
    double area=0,integral_x=0,integral_y=0;
    for_each_quadrature_point(rule,[&](Point2 p,double w){
      ++count;
      area+=w;integral_x+=w*p.x;integral_y+=w*p.y;
      if(!(p.x>=0&&p.x<=2&&p.y>=0&&p.y<=1&&w>0))
        throw std::runtime_error("Outside quadrature cell");
    });
    if(count!=8*4*order*order)
      throw std::runtime_error("Wrong Ch5 quadrature points per integration cell");
    near(area,2,1e-13,"Ch5 2x1 area");
    near(integral_x,2,1e-13,"Ch5 x first moment");
    near(integral_y,1,1e-13,"Ch5 y first moment");
  }
  // Five-point Gauss-Legendre in each integration axis is exact up to
  // degree nine for polynomials on each individual background cell.
  double high_order=0;
  for_each_quadrature_point({rect,8,4,5},[&](Point2 p,double w){
    high_order+=w*std::pow(p.x,8)*std::pow(p.y,8);
  });
  near(high_order,(std::pow(2.0,9)/9.0)/9.0,3e-12,
       "Ch5 Gauss5 polynomial degree 8 integral");
  // The 1x1 rule intentionally cannot integrate even a quadratic
  // on a cell exactly, so this MUST be tested as a different rule.
  double midpoint_x2=0;
  for_each_quadrature_point({rect,8,4,1},[&](Point2 p,double w){
    midpoint_x2+=w*p.x*p.x;
  });
  if(std::abs(midpoint_x2-8.0/3.0)<1e-4)
    throw std::runtime_error("One-point Gauss unexpectedly exact on x²");
}

void test_cartesian_assembly(int degree) {
  const Rectangle rect{0,2,0,1};
  const auto nodes=make_regular_cloud(rect,5,5,3.0);
  const MLSEvaluator mls(nodes,degree);
  const BackgroundQuadrature quadrature{rect,4,4,3};
  const CartesianAssembly a=assemble_cartesian(mls,quadrature);
  const auto N=static_cast<Eigen::Index>(nodes.size());
  if(a.stiffness.rows()!=2*N || a.mass.rows()!=2*N)
    throw std::runtime_error("Unexpected matrix DOFs");
  if(a.quadrature_points!=144 || a.max_active_nodes!=nodes.size())
    throw std::runtime_error("Unexpected integration diagnostics");
  if(!(a.max_moment_condition>0.0) || !std::isfinite(a.max_moment_condition))
    throw std::runtime_error("Invalid moment condition");

  const Eigen::MatrixXd K=Eigen::MatrixXd(a.stiffness);
  const Eigen::MatrixXd M=Eigen::MatrixXd(a.mass);
  if(!K.allFinite() || !M.allFinite())
    throw std::runtime_error("Non-finite assembled matrices");
  near((K-K.transpose()).cwiseAbs().maxCoeff(),0.0,2e-12,"K symmetry");
  near((M-M.transpose()).cwiseAbs().maxCoeff(),0.0,2e-12,"M symmetry");
  near(M.block(0,N,N,N).norm(),0.0,2e-13,"Mxy zero");
  near(M.block(N,0,N,N).norm(),0.0,2e-13,"Myx zero");

  Eigen::VectorXd constant_x=Eigen::VectorXd::Zero(2*N),
                  constant_y=Eigen::VectorXd::Zero(2*N),
                  rotation=Eigen::VectorXd::Zero(2*N),
                  gradient=Eigen::VectorXd::Zero(2*N);
  for(Eigen::Index i=0;i<N;++i){
    constant_x[i]=1.0;
    constant_y[N+i]=1.0;
    const Point2 p=nodes[static_cast<std::size_t>(i)].position;
    rotation[i]=-p.y; rotation[N+i]=p.x;  // curl=2, div=0
    gradient[i]=p.x; gradient[N+i]=p.y;  // curl=0
  }

  near(quadratic(a.mass,constant_x),2.0,1e-9,"constant x L2");
  near(quadratic(a.mass,constant_y),2.0,1e-9,"constant y L2");
  near(quadratic(a.stiffness,constant_x),0.0,1e-9,"constant x curl");
  near(quadratic(a.stiffness,constant_y),0.0,1e-9,"constant y curl");
  near(quadratic(a.stiffness,rotation),8.0,1e-8,"curl=2 energy");
  near(quadratic(a.mass,rotation),10.0/3.0,1e-8,"rotation L2");
  near(quadratic(a.stiffness,gradient),0.0,1e-8,"gradient kernel");

  Eigen::VectorXd trial=Eigen::VectorXd::Zero(2*N);
  for(Eigen::Index i=0;i<trial.size();++i)
    trial[i]=std::sin(0.17*static_cast<double>(i+1));
  if(quadratic(a.mass,trial)<=0.0)
    throw std::runtime_error("Negative or zero mass for nonzero trial");
  if(quadratic(a.stiffness,trial)<-1e-8)
    throw std::runtime_error("K is not positive semidefinite");
}

} // namespace
int main() {
  try {
    test_quadrature();
    test_author_ch5_quadrature();
    test_cartesian_assembly(1);
    test_cartesian_assembly(2);
    std::cout<<"PASS: rectangle quadrature, Cartesian K/M symmetry, constant and rotation patches, curl-grad kernel (degree 1/2), Ch5 1..5 Gauss and four subcells\n";
  }catch(const std::exception& e) {
    std::cerr<<"FAIL: "<<e.what()<<"\n";
    return 1;
  }
  return 0;
}
