// SPDX-License-Identifier: MIT
#include "al_maleh/rot_irrot.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {
void demand(bool ok,const char* label){
  if(!ok)throw std::runtime_error(label);
}
void near(double x,double target,double eps,const char* label){
  if(!std::isfinite(x)||std::abs(x-target)>eps)
    throw std::runtime_error(std::string(label)+": "+std::to_string(x));
}
void patch(int degree) {
  const Rectangle domain{0,2,0,1};
  const auto nodes=make_regular_cloud(domain,5,3,1.65);
  const MLSEvaluator mls(nodes,degree);
  const BackgroundQuadrature q{domain,8,4,3};
  const auto a=assemble_rot_irrot(mls,q);
  const auto n=static_cast<Eigen::Index>(nodes.size());
  demand(a.quadrature_points==288,"Wrong four-cell Ch6 quadrature count");
  near((a.stiffness-a.stiffness.transpose()).norm(),0,1e-8,
       "M5 K symmetry");
  near((a.mass-a.mass.transpose()).norm(),0,1e-8,"M5 M symmetry");
  near(a.stiffness.topRows(n).norm(),0,1e-12,
       "curl grad must be zero (full K I block)");
  near(a.stiffness.leftCols(n).norm(),0,1e-12,
       "K coupling curl grad must vanish");
  near((a.mass.topLeftCorner(n,n)-a.mass.bottomRightCorner(n,n)).norm(),
       0,1e-10,"M_II and M_RR should match");
  near((a.mass.topRightCorner(n,n)+a.mass.topRightCorner(n,n).transpose()).norm(),
       0,1e-10,"IR mass cross block must be skew-symmetric");
  Eigen::VectorXd random=Eigen::VectorXd::LinSpaced(2*n,-0.2,0.8);
  demand(random.dot(a.mass*random)>=-1e-11,"Mass must be semidefinite");
  demand(random.dot(a.stiffness*random)>=-1e-11,"K must be semidefinite");
  // Discrete partition of unity: sum_i ∇Ni = 0, and J∇ constant = 0.
  Eigen::VectorXd constant=Eigen::VectorXd::Zero(2*n);
  constant.head(n).setOnes();
  near((a.mass*constant).norm(),0,1e-8,
       "Gradient of scalar partition of unity must be mass-null");
  constant.head(n).setZero();constant.tail(n).setOnes();
  near((a.mass*constant).norm(),0,1e-8,
       "Rotation of scalar partition of unity must be mass-null");
  // A vector generated from affine scalar field x has curl=0 and
  // normalized mass equal to the domain area.
  Eigen::VectorXd affine=Eigen::VectorXd::Zero(2*n);
  for(Eigen::Index i=0;i<n;++i)
    affine[i]=nodes[static_cast<std::size_t>(i)].position.x;
  near(affine.dot(a.mass*affine),2,2e-7,"Mass norm of grad x");
  near(affine.dot(a.stiffness*affine),0,1e-11,
       "Exact irrotational patch stiffness");
  const auto f=evaluate_rot_irrot(mls.evaluate({0.63,0.42}),affine,n);
  near(f.hx,1,3e-8,"grad x first component");
  near(f.hy,0,3e-8,"grad x second component");
  near(f.curl,0,2e-10,"curl grad x");
  near(f.divergence,0,2e-7,"lap grad x");
  // J∇(x²+y²)/2 = (-y,x), curl = 2, div = 0;
  // exact for quadratic polynomial reproduction (degree=2).
  if(degree==2){
    Eigen::VectorXd whirl=Eigen::VectorXd::Zero(2*n);
    for(Eigen::Index i=0;i<n;++i) {
      const auto p=nodes[static_cast<std::size_t>(i)].position;
      whirl[n+i]=0.5*(p.x*p.x+p.y*p.y);
    }
    const auto f2=evaluate_rot_irrot(mls.evaluate({0.67,0.39}),whirl,n);
    near(f2.rot_x,-0.39,1e-7,"Jgrad quadratic x");
    near(f2.rot_y,0.67,1e-7,"Jgrad quadratic y");
    near(f2.curl,2,1e-6,"curl Jgrad equals Laplacian 2");
    near(f2.divergence,0,1e-12,"div Jgrad exactly zero");
    near(whirl.dot(a.stiffness*whirl),8,1e-6,
         "Original Ch6 rotational patch energy");
  }
  const auto s=solve_rot_irrot_rank_reduced(a);
  demand(s.removed_mass_nullity>=2,
         "M5 exact partition-of-unity mass nullity not detected");
  if(degree==1){
    // Baseline of source Fig. 9/Ch6 Table 4 with our declared MLS:
    // four exactly dependent directions plus fourteen curl-null
    // eigenvectors in the retained 26-dimensional mass space.
    demand(s.removed_mass_nullity==4,
           "M5 affine MLS should expose four mass dependencies");
    demand(s.retained_rank==26,
           "M5 3x5 original Ch6 rank changed");
    demand(s.numerical_curl_kernel==14,
           "M5 3x5 fourteen curl-zero eigenmodes changed");
  }
  demand(s.retained_rank+s.removed_mass_nullity==2*n,
         "Missing mass ranks");
  demand(s.eigenvalues.size()==s.retained_rank,
         "Reduced spectrum size");
  demand(s.M_orthogonality_error<1e-6,
         "Reduced eigenvectors must be mass orthonormal");
  demand(s.eigenvalues.allFinite(),
         "M5 eigenvalues must be finite");
  near(s.eigenvalues.maxCoeff(),s.eigenvalues.maxCoeff(),1e-12,"EVP finite");
  for(auto x:s.reduced_residuals)demand(std::isfinite(x)&&x<1e-7,
         "Reduced EVP residual must be tiny");
  std::cout<<"Ch6 degree="<<degree<<" full_dofs="<<2*n
           <<" retained_mass_rank="<<s.retained_rank
           <<" mass_nullity="<<s.removed_mass_nullity
           <<" zero="<<s.numerical_curl_kernel<<"\n";
}

// Theorem (exact, any cloud/quadrature with polynomial reproduction): if
// f=phi+i psi is a holomorphic polynomial reproduced by the MLS space, then
// H=grad(phi)+z x grad(psi)=0 identically (Cauchy-Riemann) and psi is
// harmonic, so the coefficient vector lies in ker(M) AND ker(K): the pencil
// (K,M) is singular. Degree 1: f=1,i,z,iz (4 vectors); degree 2 adds z^2, iz^2.
void singular_pencil_theorem(int degree) {
  const Rectangle domain{0,2,0,1};
  const auto nodes=make_regular_cloud(domain,5,3,1.65);
  const MLSEvaluator mls(nodes,degree);
  const auto a=assemble_rot_irrot(mls,{domain,8,4,3});
  const auto n=static_cast<Eigen::Index>(nodes.size());
  std::vector<Eigen::VectorXd> kernel;
  auto add=[&](auto phi,auto psi){
    Eigen::VectorXd z(2*n);
    for(Eigen::Index i=0;i<n;++i){
      const auto p=nodes[static_cast<std::size_t>(i)].position;
      z[i]=phi(p.x,p.y); z[n+i]=psi(p.x,p.y);
    }
    kernel.push_back(z);
  };
  add([](double,double){return 1.0;},[](double,double){return 0.0;});
  add([](double,double){return 0.0;},[](double,double){return 1.0;});
  add([](double x,double){return x;},[](double,double y){return y;});        // z
  add([](double,double y){return -y;},[](double x,double){return x;});       // i z
  if(degree==2){
    add([](double x,double y){return x*x-y*y;},[](double x,double y){return 2*x*y;});  // z^2
    add([](double x,double y){return -2*x*y;},[](double x,double y){return x*x-y*y;}); // i z^2
  }
  for(const auto& z:kernel){
    near((a.mass*z).norm()/(a.mass.norm()*z.norm()),0,1e-13,"holomorphic vector not mass-null");
    near((a.stiffness*z).norm()/(a.stiffness.norm()*z.norm()),0,1e-13,
         "holomorphic vector not stiffness-null (pencil would not be singular)");
  }
  const auto s=solve_rot_irrot_rank_reduced(a);
  demand(s.removed_mass_nullity==static_cast<int>(kernel.size()),
         "mass nullity differs from the holomorphic count for this cloud");
  // ker(M) inside ker(K) => projected and original finite spectra coincide.
  demand(s.mass_kernel_stiffness_coupling<1e-12,"||KZ||/||K|| certificate failed");
  demand(s.max_full_residual<1e-12,"full residual K h - lambda M h not certified");
  demand(s.max_reduced_residual<1e-12,"projected residual not small");
}

// The projected residual Q^T(Kh-lambda Mh) is ~0 BY CONSTRUCTION. Discarding
// genuine mass directions (cutoff 1e-2) changes the problem: only the full
// residual and the coupling ||KZ|| reveal it.
void projected_residual_is_not_a_certificate() {
  const Rectangle domain{0,2,0,1};
  const MLSEvaluator mls(make_regular_cloud(domain,5,3,1.65),1);
  const auto a=assemble_rot_irrot(mls,{domain,8,4,2});
  const auto bad=solve_rot_irrot_rank_reduced(a,1e-2,0.0);
  demand(bad.max_reduced_residual<1e-12,"projected residual should stay tiny");
  demand(bad.max_full_residual>1e-4,"full residual should expose the over-reduction");
  demand(bad.mass_kernel_stiffness_coupling>1e-3,"coupling should expose the over-reduction");
  // The default cutoff sits on the plateau [1e-15, mu_min/mu_max ~ 1e-4].
  const auto good=solve_rot_irrot_rank_reduced(a);
  const auto tight=solve_rot_irrot_rank_reduced(a,1e-14,0.0);
  demand(good.retained_rank==tight.retained_rank,"rank changes inside the plateau");
  near((good.eigenvalues-tight.eigenvalues).cwiseAbs().maxCoeff(),0,1e-10,
       "spectrum changes inside the rank plateau");
}

// Independent closed-form pencils: these checks do not use MLS assembly or
// compare two calls to the same eigensolver as their mathematical reference.
void exact_pencil_certificates() {
  // A dense change of orthonormal coordinates must preserve the quotient
  // spectrum {3,5} and a common two-dimensional nullspace.
  Eigen::Vector4d u(1,2,3,4);
  u.normalize();
  const Eigen::Matrix4d U=Eigen::Matrix4d::Identity()-2*u*u.transpose();
  RotIrrotAssembly compatible;
  compatible.mass=U*Eigen::Vector4d(1,2,0,0).asDiagonal()*U.transpose();
  compatible.stiffness=U*Eigen::Vector4d(3,10,0,0).asDiagonal()*U.transpose();
  const auto good=solve_rot_irrot_rank_reduced(compatible);
  demand(good.retained_rank==2&&good.removed_mass_nullity==2,
         "Closed-form common-kernel rank is wrong");
  near(good.eigenvalues[0],3,2e-13,"Common-kernel first finite eigenvalue");
  near(good.eigenvalues[1],5,2e-13,"Common-kernel second finite eigenvalue");
  demand(good.mass_kernel_stiffness_coupling<1e-14&&good.max_full_residual<1e-14,
         "Common-kernel pencil certificates failed");

  // M=diag(1,1,0), K positive definite with K_13=1. Eliminating the third
  // row gives v3=-v1/2, hence finite eigenvalues {3/2,3}. Orthogonal mass
  // projection instead gives {2,3}: it cannot enforce that third equation.
  RotIrrotAssembly coupled;
  coupled.mass=Eigen::Vector3d(1,1,0).asDiagonal();
  coupled.stiffness.resize(3,3);
  coupled.stiffness<<2,0,1, 0,3,0, 1,0,2;
  const auto bad=solve_rot_irrot_rank_reduced(coupled);
  near(bad.eigenvalues[0],2,1e-14,"Coupled projection eigenvalue");
  near(bad.eigenvalues[1],3,1e-14,"Uncoupled projection eigenvalue");
  demand(bad.max_reduced_residual<1e-14,"Coupled projected residual");
  near(bad.max_full_residual,1/(std::sqrt(19.0)+2*std::sqrt(2.0)),1e-14,
       "Full residual must retain discarded-row equation");
  near(bad.mass_kernel_stiffness_coupling,std::sqrt(5.0/19.0),1e-14,
       "Closed-form stiffness/kernel coupling");
  const Eigen::Vector3d exact_vector(1,0,-0.5);
  near((coupled.stiffness*exact_vector-1.5*coupled.mass*exact_vector).norm(),
       0,1e-14,"Schur-complement eigenpair of original pencil");

  // Even zero coupling AND zero full residual cannot establish completeness:
  // an excessive mass cutoff may silently remove a genuine finite zero mode.
  RotIrrotAssembly truncated;
  truncated.mass=Eigen::Vector3d(1,1,1e-4).asDiagonal();
  truncated.stiffness=Eigen::Vector3d(2,3,0).asDiagonal();
  const auto full=solve_rot_irrot_rank_reduced(truncated);
  const auto cut=solve_rot_irrot_rank_reduced(truncated,1e-3,0);
  demand(full.retained_rank==3&&full.numerical_curl_kernel==1,
         "Positive-mass zero mode must be retained at default cutoff");
  demand(cut.retained_rank==2&&cut.numerical_curl_kernel==0,
         "Excessive cutoff should expose lost finite zero mode");
  near(cut.mass_kernel_stiffness_coupling,0,1e-15,"Truncated zero mode has KZ=0");
  near(cut.max_full_residual,0,1e-15,"Remaining eigenpairs still solve full pencil");
  near(cut.mass_eigen_max_removed,1e-4,1e-15,
       "Mass diagnostic must reveal discarded positive direction");
}

// Al Maleh (2009) Table 4, 3x5 grid, rotational/irrotational, 2x2 Gauss per
// cell, four cells per nodal rectangle: "14 zeros", 5.91, 12.9, 15.0, 15.5,
// 16.8 (TM11, bold). The thesis gives Eq. (4-2) but no radius: r=1.5h is an
// inference of this audit (unique minimum of an r/h scan against Table 4).
void author_profile_table4() {
  const Rectangle domain{0,2,0,1};
  const MLSEvaluator mls(make_regular_cloud(domain,5,3,0.75),1,1e12,
                         MLSWeight::QuarticSpline);
  const auto s=solve_rot_irrot_rank_reduced(assemble_rot_irrot(mls,{domain,8,4,2}));
  demand(s.numerical_curl_kernel==14,"author profile should give 14 zeros");
  const double table4[]={5.91,12.9,15.0,15.5,16.8};
  for(int k=0;k<5;++k){
    const double v=s.eigenvalues[14+k];
    if(std::abs(v/table4[k]-1)>0.005)
      throw std::runtime_error("Table 4 rot/irrot value "+std::to_string(k+2)+
                               " not reproduced: "+std::to_string(v));
  }
  const double pi=std::acos(-1.0);
  near(100*(s.eigenvalues[18]/(1.25*pi*pi)-1),36.183,0.002,"Table 5 TM11 3x5 2x2 error");
}
}
int main(){
  try{patch(1);patch(2);
    singular_pencil_theorem(1);singular_pencil_theorem(2);
    projected_residual_is_not_a_certificate();
    exact_pencil_certificates();
    author_profile_table4();
    std::cout<<"PASS: original Al Maleh Ch6 gradient, Jgrad, Laplacian, "
                "mass singularity, rank reduction, exact vector patches; singular "
                "pencil theorem (deg 1/2), full-residual certificate, Table 4 profile\n";
  }catch(const std::exception& e){
    std::cerr<<"M5 FAIL: "<<e.what()<<"\n";return 1;
  }
}
