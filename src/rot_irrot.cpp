// SPDX-License-Identifier: MIT
#include "al_maleh/rot_irrot.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace al_maleh {
RotIrrotAssembly assemble_rot_irrot(const MLSEvaluator& mls,
                                    const BackgroundQuadrature& rule) {
  const auto n=static_cast<Eigen::Index>(mls.nodes().size());
  if(n<3 || n>120)
    throw std::invalid_argument("M5 dense prototype requires 3..120 nodes");
  RotIrrotAssembly a;
  a.stiffness=Eigen::MatrixXd::Zero(2*n,2*n);
  a.mass=Eigen::MatrixXd::Zero(2*n,2*n);
  for_each_quadrature_point(rule,[&](Point2 point,double w){
    const auto shape=mls.evaluate(point);
    ++a.quadrature_points;
    a.max_mls_condition=std::max(a.max_mls_condition,shape.condition_number);
    const auto count=shape.active_nodes.size();
    if(shape.dx.size()!=count||shape.dy.size()!=count||
       shape.dxx.size()!=count||shape.dxy.size()!=count||
       shape.dyy.size()!=count)
      throw std::runtime_error("Incomplete MLS gradients/Hessians for M5");
    for(std::size_t ii=0;ii<count;++ii) {
      const auto i=static_cast<Eigen::Index>(shape.active_nodes[ii]);
      const double ix=shape.dx[ii],iy=shape.dy[ii];
      const double ilap=shape.dxx[ii]+shape.dyy[ii];
      for(std::size_t jj=0;jj<count;++jj) {
        const auto j=static_cast<Eigen::Index>(shape.active_nodes[jj]);
        const double jx=shape.dx[jj],jy=shape.dy[jj];
        const double jlap=shape.dxx[jj]+shape.dyy[jj];
        // N_i^I = (dx Ni,dy Ni); N_i^R = (-dy Ni,dx Ni)
        // curl N_i^I=0; div N_i^R=0; curl N_i^R=laplacian(Ni).
        const double gg=w*(ix*jx+iy*jy);
        a.mass(i,j)+=gg;
        a.mass(n+i,n+j)+=gg;
        a.mass(i,n+j)+=w*(-ix*jy+iy*jx);
        a.mass(n+i,j)+=w*(-iy*jx+ix*jy);
        a.stiffness(n+i,n+j)+=w*ilap*jlap;
      }
    }
  });
  if(!a.mass.allFinite()||!a.stiffness.allFinite())
    throw std::runtime_error("M5 assembled matrices contain nonfinite data");
  const double mnorm=std::max(1.0,a.mass.norm());
  const double knorm=std::max(1.0,a.stiffness.norm());
  if((a.mass-a.mass.transpose()).norm()>1e-11*mnorm ||
     (a.stiffness-a.stiffness.transpose()).norm()>1e-11*knorm)
    throw std::runtime_error("M5 variational forms are not symmetric");
  return a;
}

RotIrrotField evaluate_rot_irrot(const MLSResult& shape,
    const Eigen::VectorXd& coefficient,Eigen::Index n) {
  if(n<1||coefficient.size()!=2*n||!coefficient.allFinite())
    throw std::invalid_argument("M5 modal coefficient dimension mismatch");
  RotIrrotField v;
  for(std::size_t k=0;k<shape.active_nodes.size();++k) {
    const auto i=static_cast<Eigen::Index>(shape.active_nodes[k]);
    if(i<0||i>=n)
      throw std::invalid_argument("MLS active node out of bounds");
    const double a=coefficient[i],b=coefficient[n+i];
    const double dx=shape.dx[k],dy=shape.dy[k];
    const double lap=shape.dxx[k]+shape.dyy[k];
    v.irrot_x+=dx*a;v.irrot_y+=dy*a;
    v.rot_x+=-dy*b;v.rot_y+=dx*b;
    v.divergence+=lap*a;
    v.curl+=lap*b;
  }
  v.hx=v.irrot_x+v.rot_x;v.hy=v.irrot_y+v.rot_y;
  return v;
}

RankReducedSpectrum solve_rot_irrot_rank_reduced(
    const RotIrrotAssembly& forms,double relative,double absolute) {
  if(!(relative>0&&relative<1)||!(absolute>=0)||
     !std::isfinite(relative)||!std::isfinite(absolute))
    throw std::invalid_argument("Invalid singular-mass rank cutoff");
  const Eigen::MatrixXd& K=forms.stiffness;
  const Eigen::MatrixXd& M=forms.mass;
  if(K.rows()==0||K.rows()!=K.cols()||K.rows()!=M.rows()||
     M.rows()!=M.cols()||!K.allFinite()||!M.allFinite())
    throw std::invalid_argument("Bad M5 EVP dimensions");
  if((K-K.transpose()).norm()>1e-9*std::max(1.0,K.norm())||
     (M-M.transpose()).norm()>1e-9*std::max(1.0,M.norm()))
    throw std::invalid_argument("M5 matrices not symmetric");
  const Eigen::MatrixXd A=(0.5*(K+K.transpose())).eval();
  const Eigen::MatrixXd B=(0.5*(M+M.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> se(B);
  if(se.info()!=Eigen::Success)
    throw std::runtime_error("M5 mass eigensolver failed");
  const double mu_max=se.eigenvalues().maxCoeff();
  if(!(mu_max>0)||!std::isfinite(mu_max))
    throw std::runtime_error("M5 mass matrix has no positive spectrum");
  RankReducedSpectrum out;
  out.mass_eigen_max=mu_max;
  out.mass_rank_cutoff=std::max(absolute,relative*mu_max);
  if(se.eigenvalues().minCoeff()<-1e-9*mu_max)
    throw std::runtime_error("M5 mass matrix materially indefinite");
  std::vector<Eigen::Index> indices;
  for(Eigen::Index i=0;i<se.eigenvalues().size();++i)
    if(se.eigenvalues()[i]>out.mass_rank_cutoff)
      indices.push_back(i);
  out.retained_rank=static_cast<int>(indices.size());
  out.removed_mass_nullity=static_cast<int>(M.rows())-out.retained_rank;
  if(out.retained_rank<2)
    throw std::runtime_error("M5 mass reduction leaves no usable field space");
  out.mass_eigen_min_retained=se.eigenvalues()[indices.front()];
  Eigen::MatrixXd Q(M.rows(),out.retained_rank);
  for(int i=0;i<out.retained_rank;++i)
    Q.col(i)=se.eigenvectors().col(indices[static_cast<std::size_t>(i)])/
      std::sqrt(se.eigenvalues()[indices[static_cast<std::size_t>(i)]]);
  const Eigen::MatrixXd reduced=(Q.transpose()*A*Q).eval();
  const Eigen::MatrixXd sym=(0.5*(reduced+reduced.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(sym);
  if(eigen.info()!=Eigen::Success)
    throw std::runtime_error("M5 reduced symmetric eigensolver failed");
  out.eigenvalues=eigen.eigenvalues();
  out.eigenvectors=Q*eigen.eigenvectors();
  if(!out.eigenvalues.allFinite()||!out.eigenvectors.allFinite())
    throw std::runtime_error("Nonfinite M5 eigensolution");
  const auto gram=(out.eigenvectors.transpose()*B*out.eigenvectors).eval();
  out.M_orthogonality_error=(gram-Eigen::MatrixXd::Identity(out.retained_rank,
                                            out.retained_rank)).cwiseAbs().maxCoeff();
  if(out.M_orthogonality_error>1e-5)
    throw std::runtime_error("M5 M-orthonormality failed");
  const double maxabs=out.eigenvalues.cwiseAbs().maxCoeff();
  out.eigenvalue_zero_cutoff=std::max(1e-9,1e-9*maxabs);
  const double Kn=A.norm(),Mn=B.norm();
  for(Eigen::Index j=0;j<out.eigenvalues.size();++j) {
    const auto v=out.eigenvectors.col(j).eval();
    const double lambda=out.eigenvalues[j];
    const Eigen::VectorXd residue=(A*v-lambda*B*v).eval();
    const double denominator=std::max(1e-300,
        (Kn+std::abs(lambda)*Mn)*v.norm());
    out.full_residuals.push_back(residue.norm()/denominator);
    out.reduced_residuals.push_back((Q.transpose()*residue).norm()/denominator);
    if(std::abs(lambda)<=out.eigenvalue_zero_cutoff)
      ++out.numerical_curl_kernel;
  }
  return out;
}
} // namespace al_maleh
