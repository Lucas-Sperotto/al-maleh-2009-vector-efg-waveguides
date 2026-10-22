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
// Al Maleh (2009) Table 4 (Cartesian, 3x5 grid) and Tables 2/3: quartic
// spline Eq. (4-2), linear basis, four cells per nodal rectangle. r=1.5h is
// an inference of this audit (Table 4 r/h scan), NOT a value printed in the
// thesis. 5 near-zero values, 14 small spurious values and TM11 (row 20).
void author_profile_table4_cartesian(){
  const Rectangle r{0,2,0,1};
  const double pi=std::acos(-1.0);
  const double table4[3][15]={
    {1.16e-03,9.87e-03,1.17e-02,3.94e-02,5.73e-02,2.71e-01,6.85e-01,1.05,2.53,2.66,4.68,10.9,11.3,11.9,14.5},
    {1.08e-03,1.07e-02,1.13e-02,3.96e-02,6.62e-02,2.84e-01,7.24e-01,1.11,2.53,2.82,4.76,10.9,11.6,11.9,14.5},
    {1.03e-03,9.91e-03,1.08e-02,3.81e-02,6.60e-02,2.84e-01,7.24e-01,1.11,2.52,2.82,4.76,10.8,11.6,11.9,14.5}};
  const double table2[3]={17.199,17.210,17.206};  // TM11 3x5, Gauss 2..4
  const MLSEvaluator mls(make_regular_cloud(r,5,3,0.75),1,1e12,MLSWeight::QuarticSpline);
  for(int g=2;g<=4;++g){
    const auto a=assemble_cartesian(mls,{r,8,4,g});
    const auto s=solve_generalized_symmetric(a.stiffness,a.mass);
    if(s.numerical_kernel_dimension!=5)
      throw std::runtime_error("Table 4: five near-zero Cartesian eigenvalues expected");
    for(int k=0;k<15;++k){
      const double v=s.eigenvalues[5+k];
      if(std::abs(v/table4[g-2][k]-1)>0.005)
        throw std::runtime_error("Table 4 Cartesian row "+std::to_string(k+6)+
                                 " Gauss "+std::to_string(g)+" got "+std::to_string(v));
    }
    near(100*(s.eigenvalues[19]/(1.25*pi*pi)-1),table2[g-2],0.002,"Table 2 TM11 3x5");
  }
}
}
int main(){
 try{
  known_generalized_problem();
  analytical_tm();
  assembled_spectrum();
  author_profile_table4_cartesian();
  std::cout<<"PASS: known eigenvalues, rank checks, M-orthogonality, TM references, eigenfield metrics "
             "and Al Maleh Table 4/Table 2 Cartesian spectrum with the identified Ch.4 profile\n";
 }catch(const std::exception& e){
  std::cerr<<"FAIL: "<<e.what()<<'\n';
  return 1;
 }
}
