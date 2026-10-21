// SPDX-License-Identifier: MIT
#include "al_maleh/gradient_audit.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace al_maleh;
namespace {
void demand(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
void near(double value,double expected,double tolerance,const char* label) {
  if(!std::isfinite(value)||std::abs(value-expected)>tolerance)
    throw std::runtime_error(std::string(label)+" got "+std::to_string(value)+
                             " expected "+std::to_string(expected));
}
void invalid(const std::function<void()>& f){
  try{f();}catch(const std::invalid_argument&){return;}
  throw std::runtime_error("Expected invalid_argument");
}
void patch(int degree) {
  const Rectangle rect{0,2,0,1};
  const auto nodes=make_regular_cloud(rect,5,5,3.0);
  const MLSEvaluator mls(nodes,degree);
  const BackgroundQuadrature q{rect,4,4,3};
  const auto a=assemble_cartesian(mls,q);
  const auto pair=assemble_gradient_pairing(mls,q);
  const auto projection=build_gradient_projection(a.mass,pair);
  const Eigen::Index n=nodes.size();
  demand(pair.quadrature_points==a.quadrature_points,"M6d volume quadrature mismatch");
  demand(pair.coupling.rows()==2*n&&pair.coupling.cols()==n,
         "Bad vector/scalar gradient pairing");
  demand(pair.scalar_energy.rows()==n,"Bad scalar gradient energy");
  near((pair.scalar_energy-pair.scalar_energy.transpose()).norm(),0,1e-10,
       "Scalar gradient Gram symmetry");
  demand(projection.numerical_rank>=2&&projection.numerical_rank<n,
         "Expected at least constant scalar null direction");
  demand(projection.maximum_m_orthogonality_error<1e-6,
         "M-orthogonal basis wrong");
  demand(std::isfinite(projection.relative_gradient_reconstruction_frobenius),
         "Nonfinite projection discrepancy");

  Eigen::VectorXd x(n),y(n);
  Eigen::VectorXd cx=Eigen::VectorXd::Zero(2*n);
  for(Eigen::Index i=0;i<n;++i){
    const auto p=nodes[static_cast<std::size_t>(i)].position;
    x[i]=p.x;y[i]=p.y;
    cx[i]=1;
  }
  near(x.dot(pair.scalar_energy*x),2,1e-7,"Gradient x energy");
  near(y.dot(pair.scalar_energy*y),2,1e-7,"Gradient y energy");
  demand(scalar_gradient_fit_relative_error(pair,projection,x)<1e-5,
         "Affine x gradient must fit Cartesian basis");
  demand(scalar_gradient_fit_relative_error(pair,projection,y)<1e-5,
         "Affine y gradient must fit Cartesian basis");
  const auto m=evaluate_gradient_mode(projection,a.mass,a.stiffness,cx);
  demand(m.gradient_energy_fraction>1-1e-6,
         "Constant vector x must be representable by scalar gradient");
  near(m.gradient_energy_fraction+m.complement_energy_fraction,1,1e-12,
       "Energy projection sum");
  demand(m.projected_gradient_curl_rayleigh<1e-4,
         "Affine-gradient projected vector curl leakage too large");

  const Eigen::VectorXd rhs=a.mass*cx;
  const Eigen::VectorXd projected=
      projection.basis*(projection.basis.transpose()*rhs);
  const Eigen::VectorXd projected_again=
      projection.basis*(projection.basis.transpose()*(a.mass*projected));
  near((projected-projected_again).norm()/projected.norm(),0,1e-6,
       "L2 gradient projector must be idempotent");

  // A change of scalar basis that is orthogonal must not change the
  // represented gradient range or the projected energy of any vector.
  GradientPairing transformed=pair;
  Eigen::MatrixXd R=Eigen::MatrixXd::Identity(n,n);
  const double t=0.713;
  R(0,0)=std::cos(t);R(0,1)=-std::sin(t);
  R(1,0)=std::sin(t);R(1,1)=std::cos(t);
  transformed.coupling=(pair.coupling*R).eval();
  transformed.scalar_energy=(R.transpose()*pair.scalar_energy*R).eval();
  const auto changed=build_gradient_projection(a.mass,transformed);
  demand(changed.numerical_rank==projection.numerical_rank,
         "Scalar basis rotation changed gradient numerical rank");
  const auto rotated=evaluate_gradient_mode(changed,a.mass,a.stiffness,cx);
  near(rotated.gradient_energy_fraction,m.gradient_energy_fraction,2e-6,
       "Gradient projection energy not invariant under scalar basis change");
  near(rotated.projected_gradient_curl_rayleigh,
       m.projected_gradient_curl_rayleigh,1e-5,
       "Projected-gradient curl not basis invariant");
  invalid([&] {
    GradientAuditOptions o; o.relative_rank_cutoff=0;
    (void)build_gradient_projection(a.mass,pair,o);
  });
  invalid([&] {
    (void)scalar_gradient_fit_relative_error(
        pair,projection,Eigen::VectorXd::Ones(n-1));
  });
}
void spectral() {
  const Rectangle rect{0,2,0,1};
  const MLSEvaluator mls(make_regular_cloud(rect,9,7,0.825),1);
  const BackgroundQuadrature q{rect,8,6,3};
  const auto matrices=assemble_cartesian(mls,q);
  const auto pairing=assemble_gradient_pairing(mls,q);
  const auto projection=build_gradient_projection(matrices.mass,pairing);
  const auto spectrum=solve_generalized_symmetric(matrices.stiffness,matrices.mass);
  int count=0;
  for(Eigen::Index i=0;i<spectrum.eigenvalues.size();++i){
    const auto m=evaluate_gradient_mode(
        projection,matrices.mass,matrices.stiffness,
        spectrum.eigenvectors.col(i));
    demand(m.gradient_energy_fraction>=0&&m.gradient_energy_fraction<=1,
           "Every eigenmode must have valid projected energy fraction");
    demand(std::isfinite(m.projected_gradient_curl_rayleigh),
           "Projected gradient curl energy must be finite");
    if(spectrum.eigenvalues[i]>spectrum.zero_cutoff&&
       spectrum.eigenvalues[i]<1.25*std::acos(-1.)*std::acos(-1.))
      ++count;
  }
  demand(count==49,"Historical low-frequency count changed");
  std::cout<<"M6d 9x7 projected gradient rank="<<projection.numerical_rank
           <<" low="<<count<<"\n";
}
}
int main(){
  try{
    patch(1);
    patch(2);
    spectral();
    std::cout<<"PASS: gradient L2 coupling, polynomial fit, basis invariance,"
                " projector idempotence, all modal fractions, unmodified EVP\n";
  }catch(const std::exception& e){
    std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;
  }
}
