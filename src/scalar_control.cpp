// SPDX-License-Identifier: MIT
#include "al_maleh/scalar_control.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace al_maleh {
namespace {
struct Gauss1D { std::vector<double> x, w; };
Gauss1D boundary_gauss(int order) {
  switch(order) {
    case 2:return {{-0.57735026918962576451,0.57735026918962576451},{1,1}};
    case 3:return {{-0.77459666924148337704,0,0.77459666924148337704},
                   {5.0/9.0,8.0/9.0,5.0/9.0}};
    case 4:return {{-0.86113631159405257522,-0.33998104358485626480,
                     0.33998104358485626480,0.86113631159405257522},
                   {0.34785484513745385737,0.65214515486254614263,
                    0.65214515486254614263,0.34785484513745385737}};
    default:throw std::invalid_argument("Boundary Gauss order must be 2..4");
  }
}
void check_domain(Rectangle r) {
  if(!std::isfinite(r.xmin)||!std::isfinite(r.xmax)||
     !std::isfinite(r.ymin)||!std::isfinite(r.ymax)||
     r.xmin>=r.xmax||r.ymin>=r.ymax)
    throw std::invalid_argument("Invalid scalar EFG rectangular domain");
}
struct Sample {double phi=0,dx=0,dy=0,dxy=0;};
Sample reconstruct(const MLSResult& s,const Eigen::VectorXd& coefficients){
  Sample result;
  for(std::size_t i=0;i<s.active_nodes.size();++i){
    const auto id=static_cast<Eigen::Index>(s.active_nodes[i]);
    const double c=coefficients[id];
    result.phi+=c*s.phi[i];
    result.dx+=c*s.dx[i];
    result.dy+=c*s.dy[i];
    result.dxy+=c*s.dxy[i];
  }
  return result;
}
}

ScalarControlSystem assemble_scalar_control(
    const MLSEvaluator& mls,const BackgroundQuadrature& q){
  check_domain(q.domain);
  const auto n=static_cast<Eigen::Index>(mls.nodes().size());
  if(n<4||n>600)throw std::invalid_argument("Scalar reference supports 4..600 nodes");
  ScalarControlSystem result;
  result.stiffness=Eigen::MatrixXd::Zero(n,n);
  result.mass=Eigen::MatrixXd::Zero(n,n);
  for_each_quadrature_point(q,[&](Point2 p,double weight){
    const auto s=mls.evaluate(p);
    ++result.quadrature_points;
    result.max_moment_condition=std::max(result.max_moment_condition,s.condition_number);
    for(std::size_t i=0;i<s.active_nodes.size();++i){
      const auto a=static_cast<Eigen::Index>(s.active_nodes[i]);
      for(std::size_t j=0;j<s.active_nodes.size();++j){
        const auto b=static_cast<Eigen::Index>(s.active_nodes[j]);
        result.stiffness(a,b)+=weight*(s.dx[i]*s.dx[j]+s.dy[i]*s.dy[j]);
        result.mass(a,b)+=weight*s.phi[i]*s.phi[j];
      }
    }
  });
  const double a=q.domain.xmax-q.domain.xmin,b=q.domain.ymax-q.domain.ymin;
  const double tol=1e-10*std::max({1.0,a,b});
  std::vector<Point2> points;
  for(const auto& node:mls.nodes()){
    const Point2 p=node.position;
    if(p.x<q.domain.xmin-tol||p.x>q.domain.xmax+tol||
       p.y<q.domain.ymin-tol||p.y>q.domain.ymax+tol)
      throw std::invalid_argument("Cloud node lies outside integration rectangle");
    if(std::abs(p.x-q.domain.xmin)<tol||
       std::abs(p.x-q.domain.xmax)<tol||
       std::abs(p.y-q.domain.ymin)<tol||
       std::abs(p.y-q.domain.ymax)<tol)
      points.push_back(p);
  }
  if(points.size()<4||points.size()>=static_cast<std::size_t>(n))
    throw std::invalid_argument("Expected interior and at least four PEC collocation nodes");
  result.boundary_points=points.size();
  result.boundary=Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(points.size()),n);
  for(Eigen::Index row=0;row<result.boundary.rows();++row){
    const auto s=mls.evaluate(points[static_cast<std::size_t>(row)]);
    for(std::size_t i=0;i<s.active_nodes.size();++i)
      result.boundary(row,static_cast<Eigen::Index>(s.active_nodes[i]))=s.phi[i];
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(result.boundary,Eigen::ComputeFullV);
  if(svd.info()!=Eigen::Success||svd.singularValues().size()==0||
     !(svd.singularValues()[0]>0))
    throw std::runtime_error("PEC boundary constraint SVD failed");
  const double cutoff=1e-11*svd.singularValues()[0];
  for(Eigen::Index i=0;i<svd.singularValues().size();++i)
    if(svd.singularValues()[i]>cutoff)++result.boundary_rank;
  if(result.boundary_rank>=n)
    throw std::runtime_error("No unconstrained EFG scalar DOFs");
  result.nullspace=svd.matrixV().rightCols(n-result.boundary_rank);
  result.boundary_constraint_residual=(result.boundary*result.nullspace).norm();
  if(result.boundary_constraint_residual>1e-8*
       std::max(1.0,result.boundary.norm()*result.nullspace.norm()))
    throw std::runtime_error("PEC constraint nullspace inaccurate");
  result.reduced_stiffness=
      (result.nullspace.transpose()*result.stiffness*result.nullspace).eval();
  result.reduced_mass=
      (result.nullspace.transpose()*result.mass*result.nullspace).eval();
  result.reduced_stiffness=
      (0.5*(result.reduced_stiffness+result.reduced_stiffness.transpose())).eval();
  result.reduced_mass=
      (0.5*(result.reduced_mass+result.reduced_mass.transpose())).eval();
  return result;
}

ScalarControlSpectrum solve_scalar_control(const ScalarControlSystem& system){
  const auto n=system.mass.rows(),k=system.nullspace.cols();
  if(n<1||system.mass.cols()!=n||system.stiffness.rows()!=n||
     system.stiffness.cols()!=n||system.boundary.cols()!=n||
     system.nullspace.rows()!=n||k<1||
     system.reduced_mass.rows()!=k||system.reduced_mass.cols()!=k||
     system.reduced_stiffness.rows()!=k||system.reduced_stiffness.cols()!=k)
    throw std::invalid_argument("Inconsistent constrained scalar EFG matrices");
  if(!system.mass.allFinite()||!system.stiffness.allFinite()||
     !system.boundary.allFinite()||!system.nullspace.allFinite()||
     !system.reduced_stiffness.allFinite()||!system.reduced_mass.allFinite())
    throw std::invalid_argument("Nonfinite constrained scalar EFG matrices");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> masscheck(
      system.reduced_mass,Eigen::EigenvaluesOnly);
  if(masscheck.info()!=Eigen::Success||masscheck.eigenvalues().minCoeff()<=0)
    throw std::runtime_error("Constrained scalar EFG mass not SPD");
  ScalarControlSpectrum output;
  output.mass_condition=masscheck.eigenvalues().maxCoeff()/
                        masscheck.eigenvalues().minCoeff();
  if(!std::isfinite(output.mass_condition)||output.mass_condition>1e13)
    throw std::runtime_error("Constrained scalar mass too ill-conditioned");
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      system.reduced_stiffness,system.reduced_mass);
  if(solver.info()!=Eigen::Success)
    throw std::runtime_error("Constrained scalar generalized eigensolver failed");
  output.eigenvalues=solver.eigenvalues();
  output.coefficients=system.nullspace*solver.eigenvectors();
  if(!output.eigenvalues.allFinite()||!output.coefficients.allFinite())
    throw std::runtime_error("Nonfinite constrained scalar eigensolution");
  if(output.eigenvalues.minCoeff()<=0)
    throw std::runtime_error("Non-positive scalar Dirichlet eigenvalue");
  const Eigen::MatrixXd gram=output.coefficients.transpose()*
      system.mass*output.coefficients;
  output.mass_orthogonality_error=
      (gram-Eigen::MatrixXd::Identity(k,k)).cwiseAbs().maxCoeff();
  if(output.mass_orthogonality_error>1e-6)
    throw std::runtime_error("Scalar EFG eigenvectors not M-orthonormal");
  const double Sn=system.reduced_stiffness.norm(),Mn=system.reduced_mass.norm();
  for(Eigen::Index i=0;i<k;++i){
    const double lambda=output.eigenvalues[i];
    const Eigen::VectorXd y=solver.eigenvectors().col(i);
    const double residual=(system.reduced_stiffness*y-lambda*system.reduced_mass*y).norm()/
        std::max(1e-300,(Sn+lambda*Mn)*y.norm());
    const double boundary=(system.boundary*output.coefficients.col(i)).norm()/
        output.coefficients.col(i).norm();
    if(!std::isfinite(residual)||!std::isfinite(boundary)||
       residual>1e-7||boundary>1e-7)
      throw std::runtime_error("Scalar constrained eigenpair failed algebraic or boundary test");
    output.projected_residuals.push_back(residual);
    output.collocation_residuals.push_back(boundary);
  }
  return output;
}

ScalarControlFieldMetrics evaluate_scalar_control_mode(
    const MLSEvaluator& mls,const BackgroundQuadrature& q,
    const Eigen::VectorXd& v,int m,int n,int boundary_gauss_order) {
  check_domain(q.domain);
  if(v.size()!=static_cast<Eigen::Index>(mls.nodes().size())||!v.allFinite()||
     m<1||n<1)throw std::invalid_argument("Bad scalar mode query");
  const auto g=boundary_gauss(boundary_gauss_order);
  const Rectangle r=q.domain;
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  const double pi=std::acos(-1.);
  double psi2=0,ref2=0,cross=0,h2=0,div2=0,trace2=0,normal2=0;
  for_each_quadrature_point(q,[&](Point2 p,double w){
    const auto s=reconstruct(mls.evaluate(p),v);
    const double ref=std::sin(m*pi*(p.x-r.xmin)/a)*
                     std::sin(n*pi*(p.y-r.ymin)/b);
    psi2+=w*s.phi*s.phi;
    ref2+=w*ref*ref;
    cross+=w*ref*s.phi;
    h2+=w*(s.dx*s.dx+s.dy*s.dy);
    // Analytically, div(-partial_y psi,partial_x psi)=0.
    const double div=-s.dxy+s.dxy;
    div2+=w*div*div;
  });
  if(!(psi2>0)||!(ref2>0)||!(h2>0))
    throw std::runtime_error("Zero scalar control modal energy");
  auto edge=[&](Point2 p,double nx,double ny,double w){
    const auto s=reconstruct(mls.evaluate(p),v);
    trace2+=w*s.phi*s.phi;
    const double hnormal=nx*(-s.dy)+ny*s.dx;
    normal2+=w*hnormal*hnormal;
  };
  for(int cy=0;cy<q.cells_y;++cy)for(std::size_t j=0;j<g.x.size();++j){
    const double step=b/q.cells_y;
    const double y=r.ymin+(cy+0.5*(1+g.x[j]))*step;
    edge({r.xmin,y},-1,0,0.5*step*g.w[j]);
    edge({r.xmax,y},1,0,0.5*step*g.w[j]);
  }
  for(int cx=0;cx<q.cells_x;++cx)for(std::size_t j=0;j<g.x.size();++j){
    const double step=a/q.cells_x;
    const double x=r.xmin+(cx+0.5*(1+g.x[j]))*step;
    edge({x,r.ymin},0,-1,0.5*step*g.w[j]);
    edge({x,r.ymax},0,1,0.5*step*g.w[j]);
  }
  ScalarControlFieldMetrics result;
  result.relative_shape_overlap=std::min(1.0,std::abs(cross)/std::sqrt(psi2*ref2));
  result.boundary_trace_ratio=std::sqrt((a*b/(2*(a+b)))*trace2/psi2);
  result.normal_flux_ratio=std::sqrt((a*b/(2*(a+b)))*normal2/h2);
  result.divergence_l2=std::sqrt(div2);
  return result;
}
} // namespace al_maleh
