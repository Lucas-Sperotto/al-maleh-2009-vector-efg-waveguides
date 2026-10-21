// SPDX-License-Identifier: MIT
#include "al_maleh/constraints.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace al_maleh {
namespace {
struct Gauss1D { std::vector<double> x, w; };
Gauss1D gauss(int order) {
  switch(order) {
    case 2: return {{-0.57735026918962576451,0.57735026918962576451},{1,1}};
    case 3: return {{-0.77459666924148337704,0,0.77459666924148337704},
                    {5.0/9.0,8.0/9.0,5.0/9.0}};
    case 4: return {{-0.86113631159405257522,-0.3399810435848562648,
                      0.3399810435848562648,0.86113631159405257522},
                    {0.34785484513745385737,0.65214515486254614263,
                     0.65214515486254614263,0.34785484513745385737}};
    default: throw std::invalid_argument("Boundary quadrature order must be 2, 3, or 4");
  }
}
InvariantRange range(const Eigen::MatrixXd& gram,
                     const Eigen::MatrixXd& projected,
                     double scaling) {
  const Eigen::MatrixXd symmetric = (0.5*(projected+projected.transpose())).eval();
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ev(symmetric, gram);
  if(ev.info()!=Eigen::Success)
    throw std::runtime_error("Constraint subspace generalized solve failed");
  const double min_e=ev.eigenvalues().minCoeff(),max_e=ev.eigenvalues().maxCoeff();
  const double tolerance=1e-10*std::max(1.0,max_e);
  if(min_e < -tolerance)
    throw std::runtime_error("Constraint quadratic form not positive semidefinite");
  InvariantRange output;
  output.minimum=std::sqrt(std::max(0.0,min_e)*scaling);
  output.maximum=std::sqrt(std::max(0.0,max_e)*scaling);
  output.maximizing_coefficients=ev.eigenvectors().col(ev.eigenvalues().size()-1);
  return output;
}
} // namespace

ConstraintAssembly assemble_constraints(
    const MLSEvaluator& mls,const BackgroundQuadrature& q,
    int boundary_gauss_order) {
  const std::size_t n=mls.nodes().size();
  if(n==0||n>600)throw std::invalid_argument("Constraint assembly supports 1..600 nodes");
  const auto nd=static_cast<Eigen::Index>(2*n);
  Eigen::MatrixXd D=Eigen::MatrixXd::Zero(nd,nd);
  Eigen::MatrixXd B=Eigen::MatrixXd::Zero(nd,nd);
  // Match the domain and integration discretization used by K/M.
  for_each_quadrature_point(q,[&](Point2 p,double w) {
    const MLSResult s=mls.evaluate(p);
    for(std::size_t a=0;a<s.active_nodes.size();++a) {
      const Eigen::Index i=static_cast<Eigen::Index>(s.active_nodes[a]);
      for(std::size_t b=0;b<s.active_nodes.size();++b) {
        const Eigen::Index j=static_cast<Eigen::Index>(s.active_nodes[b]);
        D(i,j)+=w*s.dx[a]*s.dx[b];
        D(i,static_cast<Eigen::Index>(n)+j)+=w*s.dx[a]*s.dy[b];
        D(static_cast<Eigen::Index>(n)+i,j)+=w*s.dy[a]*s.dx[b];
        D(static_cast<Eigen::Index>(n)+i,static_cast<Eigen::Index>(n)+j)
          +=w*s.dy[a]*s.dy[b];
      }
    }
  });
  if(q.cells_x<1||q.cells_y<1)
    throw std::invalid_argument("Invalid background grid");
  // Keep boundary quadrature segmented by integration cells: one global
  // three-point rule can severely alias higher-frequency field traces.
  const auto g=gauss(boundary_gauss_order);
  const Rectangle r=q.domain;
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  ConstraintAssembly result;
  auto integrate_edge=[&](Point2 p,double nx,double ny,double weight){
    const MLSResult s=mls.evaluate(p);
    ++result.boundary_points;
    for(std::size_t i=0;i<s.active_nodes.size();++i){
      const Eigen::Index ii=static_cast<Eigen::Index>(s.active_nodes[i]);
      for(std::size_t j=0;j<s.active_nodes.size();++j){
        const Eigen::Index jj=static_cast<Eigen::Index>(s.active_nodes[j]);
        const double basis=weight*s.phi[i]*s.phi[j];
        B(ii,jj)+=nx*nx*basis;
        B(ii,static_cast<Eigen::Index>(n)+jj)+=nx*ny*basis;
        B(static_cast<Eigen::Index>(n)+ii,jj)+=ny*nx*basis;
        B(static_cast<Eigen::Index>(n)+ii,static_cast<Eigen::Index>(n)+jj)
           +=ny*ny*basis;
      }
    }
  };
  for(int cell=0;cell<q.cells_y;++cell) {
    const double step=b/q.cells_y,origin=r.ymin+cell*step;
    for(std::size_t j=0;j<g.x.size();++j) {
      const double y=origin+0.5*step*(1+g.x[j]);
      const double w=0.5*step*g.w[j];
      integrate_edge({r.xmin,y},-1,0,w);
      integrate_edge({r.xmax,y},1,0,w);
    }
  }
  for(int cell=0;cell<q.cells_x;++cell) {
    const double step=a/q.cells_x,origin=r.xmin+cell*step;
    for(std::size_t j=0;j<g.x.size();++j) {
      const double x=origin+0.5*step*(1+g.x[j]);
      const double w=0.5*step*g.w[j];
      integrate_edge({x,r.ymin},0,-1,w);
      integrate_edge({x,r.ymax},0,1,w);
    }
  }
  result.divergence=D.sparseView(0,0);
  result.normal_boundary=B.sparseView(0,0);
  result.divergence.makeCompressed();
  result.normal_boundary.makeCompressed();
  return result;
}

SubspaceConstraintMetrics evaluate_subspace_constraints(
    const ConstraintAssembly& constraints,
    const Eigen::SparseMatrix<double>& mass,
    const Eigen::MatrixXd& V,
    double lambda_reference, Rectangle r) {
  const Eigen::Index n=mass.rows(),k=V.cols();
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  if(n<1||mass.cols()!=n||k<1||V.rows()!=n||
     constraints.divergence.rows()!=n||constraints.divergence.cols()!=n||
     constraints.normal_boundary.rows()!=n||constraints.normal_boundary.cols()!=n||
     !V.allFinite()||!std::isfinite(lambda_reference)||lambda_reference<=0||
     !std::isfinite(a)||!std::isfinite(b)||a<=0||b<=0)
    throw std::invalid_argument("Invalid constraint subspace inputs");
  Eigen::MatrixXd G=(V.transpose()*(mass*V)).eval();
  G=(0.5*(G+G.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gram_eig(G,Eigen::EigenvaluesOnly);
  if(gram_eig.info()!=Eigen::Success||gram_eig.eigenvalues().minCoeff()<=
     1e-12*std::max(1.0,gram_eig.eigenvalues().maxCoeff()))
    throw std::invalid_argument("Singular or nearly dependent M-subspace basis");
  SubspaceConstraintMetrics result;
  result.dimension=k;
  result.divergence_ratio=range(G,V.transpose()*(constraints.divergence*V),
                               1.0/lambda_reference);
  result.normal_ratio=range(G,V.transpose()*(constraints.normal_boundary*V),
                            (a*b)/(2*(a+b)));
  return result;
}
} // namespace al_maleh
