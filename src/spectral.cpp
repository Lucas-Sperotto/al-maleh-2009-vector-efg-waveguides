// SPDX-License-Identifier: MIT
#include "al_maleh/spectral.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace al_maleh {
namespace {
void field(const MLSResult& s,const Eigen::VectorXd& c,Eigen::Index n,
           double& hx,double& hy,double& div,double& curl){
  hx=hy=div=curl=0;
  for(std::size_t a=0;a<s.active_nodes.size();++a){
    const auto i=static_cast<Eigen::Index>(s.active_nodes[a]);
    hx+=s.phi[a]*c[i];hy+=s.phi[a]*c[n+i];
    div+=s.dx[a]*c[i]+s.dy[a]*c[n+i];
    curl+=s.dx[a]*c[n+i]-s.dy[a]*c[i];
  }
}
}
Spectrum solve_generalized_symmetric(const Eigen::SparseMatrix<double>& ks,
    const Eigen::SparseMatrix<double>& ms,const SpectralOptions& options){
  if(!std::isfinite(options.zero_absolute)||options.zero_absolute<0||
     !std::isfinite(options.zero_relative)||options.zero_relative<0||
     !std::isfinite(options.max_mass_condition)||
     options.max_mass_condition<1||!std::isfinite(options.symmetry_tolerance)||
     options.symmetry_tolerance<=0||options.max_dense_dofs<1)
    throw std::invalid_argument("Invalid spectrum options");
  if(ks.rows()<=0||ks.rows()!=ks.cols()||ms.rows()!=ks.rows()||
     ms.cols()!=ks.cols())
    throw std::invalid_argument("K/M must be matching square matrices");
  if(ks.rows()>options.max_dense_dofs)
    throw std::invalid_argument("Dense reference solver dimension exceeded");
  Eigen::MatrixXd K(ks),M(ms);
  if(!K.allFinite()||!M.allFinite())
    throw std::invalid_argument("Nonfinite stiffness or mass matrix");
  if((K-K.transpose()).norm()>
       options.symmetry_tolerance*std::max(1.0,K.norm())||
     (M-M.transpose()).norm()>
       options.symmetry_tolerance*std::max(1.0,M.norm()))
    throw std::invalid_argument("Nonsymmetric input K or M");
  const Eigen::MatrixXd A=(0.5*(K+K.transpose())).eval();
  const Eigen::MatrixXd B=(0.5*(M+M.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mass(B,Eigen::EigenvaluesOnly);
  if(mass.info()!=Eigen::Success)
    throw std::runtime_error("Cannot compute mass spectrum");
  Spectrum result;
  result.min_mass_eigenvalue=mass.eigenvalues().minCoeff();
  if(!(result.min_mass_eigenvalue>0))
    throw std::runtime_error("Non-positive-definite mass; inspect rank");
  result.mass_condition=mass.eigenvalues().maxCoeff()/result.min_mass_eigenvalue;
  if(!std::isfinite(result.mass_condition)||
     result.mass_condition>options.max_mass_condition)
    throw std::runtime_error("Mass condition number beyond configured limit");
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen(A,B);
  if(eigen.info()!=Eigen::Success)
    throw std::runtime_error("Generalized symmetric eigensolver failed");
  result.eigenvalues=eigen.eigenvalues();
  result.eigenvectors=eigen.eigenvectors();
  if(!result.eigenvalues.allFinite()||!result.eigenvectors.allFinite())
    throw std::runtime_error("Nonfinite eigensolution");
  result.zero_cutoff=std::max(options.zero_absolute,
      options.zero_relative*std::max(0.0,result.eigenvalues.maxCoeff()));
  if(result.eigenvalues.minCoeff() < -result.zero_cutoff)
    throw std::runtime_error("Materially negative eigenvalue of curl-curl");
  const double kn=K.norm(),mn=M.norm();
  for(Eigen::Index j=0;j<result.eigenvalues.size();++j){
    const Eigen::VectorXd v=result.eigenvectors.col(j);
    const double lambda=result.eigenvalues[j];
    const double denom=std::max(1e-300,(kn+std::abs(lambda)*mn)*v.norm());
    const double residual=(A*v-lambda*B*v).norm()/denom;
    if(!std::isfinite(residual)||std::abs(v.dot(B*v)-1)>1e-5)
      throw std::runtime_error("Eigenpair residual or M normalization invalid");
    result.residuals.push_back(residual);
    if(lambda<=result.zero_cutoff)++result.numerical_kernel_dimension;
  }
  return result;
}
std::vector<TMMode> rectangular_tm_reference(
    Rectangle r,int mmax,int nmax){
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  if(!std::isfinite(a)||!std::isfinite(b)||a<=0||b<=0||
     mmax<1||nmax<1||mmax>100||nmax>100)
    throw std::invalid_argument("Invalid rectangular TM benchmark");
  std::vector<TMMode> out;
  const double pi=std::acos(-1.0);
  for(int m=1;m<=mmax;++m)for(int n=1;n<=nmax;++n){
    const double lambda=std::pow(m*pi/a,2)+std::pow(n*pi/b,2);
    out.push_back({m,n,lambda,std::sqrt(lambda)});
  }
  std::sort(out.begin(),out.end(),
    [](const TMMode& x,const TMMode& y){
      if(x.lambda!=y.lambda)return x.lambda<y.lambda;
      return x.m<y.m||(x.m==y.m&&x.n<y.n);
    });
  return out;
}
FieldMetrics modal_diagnostics(const MLSEvaluator& mls,
    const BackgroundQuadrature& rule,
    const Eigen::VectorXd& v,double lambda){
  const auto n=static_cast<Eigen::Index>(mls.nodes().size());
  if(v.size()!=2*n||!v.allFinite()||!std::isfinite(lambda)||lambda<=0)
    throw std::invalid_argument("Invalid eigenfield inputs");
  double field2=0,div2=0,curl2=0;
  for_each_quadrature_point(rule,[&](Point2 p,double w){
    const auto s=mls.evaluate(p);
    double hx,hy,div,curl;field(s,v,n,hx,hy,div,curl);
    field2+=w*(hx*hx+hy*hy);
    div2+=w*div*div;curl2+=w*curl*curl;
  });
  if(!(field2>0)||!std::isfinite(field2))
    throw std::runtime_error("Nonpositive modal field norm");
  const auto r=rule.domain;
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  const double points[3]={-0.7745966692414834,0,0.7745966692414834};
  const double weights[3]={5.0/9.0,8.0/9.0,5.0/9.0};
  double boundary2=0;
  for(int q=0;q<3;++q){
    const double t=0.5*(1+points[q]);
    for(double x:{r.xmin,r.xmax}){
      const auto s=mls.evaluate({x,r.ymin+t*b});
      double hx,hy,div,curl;field(s,v,n,hx,hy,div,curl);
      boundary2+=0.5*b*weights[q]*hx*hx;
    }
    for(double y:{r.ymin,r.ymax}){
      const auto s=mls.evaluate({r.xmin+t*a,y});
      double hx,hy,div,curl;field(s,v,n,hx,hy,div,curl);
      boundary2+=0.5*a*weights[q]*hy*hy;
    }
  }
  FieldMetrics d;
  d.field_l2=std::sqrt(field2);
  d.curl_l2=std::sqrt(curl2);
  d.div_l2=std::sqrt(div2);
  d.divergence_ratio=std::sqrt(div2/(lambda*field2));
  d.normal_ratio=std::sqrt(std::max(0.0,boundary2*(a*b/(2*(a+b)))/field2));
  return d;
}
} // namespace al_maleh
