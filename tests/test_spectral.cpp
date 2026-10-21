// SPDX-License-Identifier: MIT
#include "al_maleh/spectral.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace al_maleh;
namespace {
void near(double x,double y,double tol,const char* label){
  if(!std::isfinite(x)||std::abs(x-y)>tol)
    throw std::runtime_error(std::string(label)+
      " expected="+std::to_string(y)+" got="+std::to_string(x));
}
Eigen::SparseMatrix<double> diagonal(double a,double b,double c){
  Eigen::SparseMatrix<double> out(3,3);
  out.insert(0,0)=a;out.insert(1,1)=b;out.insert(2,2)=c;
  out.makeCompressed();return out;
}
void must_fail(const std::function<void()>& f,const char* label){
  try{f();}catch(const std::exception&){return;}
  throw std::runtime_error(std::string(label)+" should have raised exception");
}
void known_generalized_problem(){
  const auto K=diagonal(0,4,10),M=diagonal(1,2,2);
  const Spectrum s=solve_generalized_symmetric(K,M);
  if(s.numerical_kernel_dimension!=1)
    throw std::runtime_error("Incorrect numerical kernel");
  near(s.eigenvalues[0],0,1e-12,"lambda0");
  near(s.eigenvalues[1],2,1e-12,"lambda1");
  near(s.eigenvalues[2],5,1e-12,"lambda2");
  near(s.mass_condition,2,1e-12,"mass condition");
  for(double r:s.residuals)near(r,0,1e-12,"eigenpair residual");
  const Eigen::MatrixXd gram=s.eigenvectors.transpose()*
     Eigen::MatrixXd(M)*s.eigenvectors;
  near((gram-Eigen::MatrixXd::Identity(3,3)).norm(),0,1e-12,
       "M orthogonality");
  must_fail([&]{(void)solve_generalized_symmetric(K,diagonal(1,0,2));},
            "singular mass");
  must_fail([&]{(void)solve_generalized_symmetric(K,diagonal(1,-2,2));},
            "indefinite mass");
  Eigen::SparseMatrix<double> badK=K;badK.coeffRef(0,1)=1;
  must_fail([&]{(void)solve_generalized_symmetric(badK,M);},
            "asymmetric stiffness");
  SpectralOptions opts;opts.max_dense_dofs=2;
  must_fail([&]{(void)solve_generalized_symmetric(K,M,opts);},
            "dense problem dimension limit");
}
void analytical_tm(){
  const auto m=rectangular_tm_reference({0,2,0,1},4,4);
  const double pi=std::acos(-1.0);
  if(m.size()!=16||m[0].m!=1||m[0].n!=1)
    throw std::runtime_error("TM11 order");
  near(m[0].lambda,1.25*pi*pi,1e-12,"TM11 analytical");
  near(m[1].lambda,2*pi*pi,1e-12,"TM21 analytical");
}
void assembled_spectrum(){
  const Rectangle r{0,2,0,1};
  const MLSEvaluator mls(make_regular_cloud(r,5,5,1.65),1);
  const BackgroundQuadrature quad{r,4,4,3};
  const auto assembled=assemble_cartesian(mls,quad);
  const auto s=solve_generalized_symmetric(
       assembled.stiffness,assembled.mass);
  if(s.eigenvalues.size()!=50||s.numerical_kernel_dimension<2||
     s.numerical_kernel_dimension>=50)
    throw std::runtime_error("Unexpected assembled spectrum/kernel");
  const int idx=s.numerical_kernel_dimension;
  if(s.eigenvalues[idx]<=s.zero_cutoff)
    throw std::runtime_error("Missing positive eigenvalue candidate");
  near(s.residuals[static_cast<std::size_t>(idx)],0,1e-6,
       "assembled residual");
  const auto diag=modal_diagnostics(
       mls,quad,s.eigenvectors.col(idx),s.eigenvalues[idx]);
  near(diag.field_l2,1,1e-5,"field L2 normalization");
  if(!std::isfinite(diag.divergence_ratio)||
     !std::isfinite(diag.normal_ratio)||diag.curl_l2<=0)
    throw std::runtime_error("Bad field metrics");
}
}
int main(){
 try{
  known_generalized_problem();
  analytical_tm();
  assembled_spectrum();
  std::cout<<"PASS: known eigenvalues, rank checks, M-orthogonality, TM references and eigenfield metrics\n";
 }catch(const std::exception& e){
  std::cerr<<"FAIL: "<<e.what()<<'\n';
  return 1;
 }
}
