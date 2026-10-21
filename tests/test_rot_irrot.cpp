// SPDX-License-Identifier: MIT
#include "al_maleh/rot_irrot.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

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
}
int main(){
  try{patch(1);patch(2);
    std::cout<<"PASS: original Al Maleh Ch6 gradient, Jgrad, Laplacian, "
                "mass singularity, rank reduction, exact vector patches\n";
  }catch(const std::exception& e){
    std::cerr<<"M5 FAIL: "<<e.what()<<"\n";return 1;
  }
}
