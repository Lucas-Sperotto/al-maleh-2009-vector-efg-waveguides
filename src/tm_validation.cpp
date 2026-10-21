// SPDX-License-Identifier: MIT
#include "al_maleh/tm_validation.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace al_maleh {
namespace {
void validate_options(const MatchOptions& o) {
  const double values[]={o.degeneracy_relative_tolerance,o.maximum_relative_error,
    o.minimum_shape_overlap,o.maximum_divergence_ratio,o.maximum_normal_ratio,
    o.maximum_algebraic_residual};
  for(double d:values)if(!std::isfinite(d)||d<0)
    throw std::invalid_argument("Invalid negative/nonfinite matching threshold");
  if(o.minimum_shape_overlap>1)
    throw std::invalid_argument("Shape overlap must be in [0,1]");
}
double relative_error(double numerical,double exact) {
  return std::abs(numerical-exact)/exact;
}
}
std::vector<TMShapeMatch> validate_tm_modes(
  const MLSEvaluator& mls,const BackgroundQuadrature& integration,
  const CartesianAssembly& assembly,const Spectrum& spectrum,
  int max_m,int max_n,const MatchOptions& options) {
  validate_options(options);
  const auto n=static_cast<Eigen::Index>(mls.nodes().size());
  const auto ndof=2*n;
  if(n==0||assembly.stiffness.rows()!=ndof||assembly.mass.rows()!=ndof||
     spectrum.eigenvalues.size()!=ndof||
     spectrum.eigenvectors.rows()!=ndof||
     spectrum.eigenvectors.cols()!=ndof||
     spectrum.residuals.size()!=static_cast<std::size_t>(ndof))
    throw std::invalid_argument("Mismatched geometry/assembled matrices/spectrum");
  const Rectangle r=integration.domain;
  const double a=r.xmax-r.xmin,b=r.ymax-r.ymin;
  const double pi=std::acos(-1.0);
  const auto reference=rectangular_tm_reference(r,max_m,max_n);
  std::vector<TMShapeMatch> matches;
  matches.reserve(reference.size());
  Eigen::VectorXd curl_norms=Eigen::VectorXd::Zero(ndof);
  for(Eigen::Index j=0;j<ndof;++j) {
    const auto v=spectrum.eigenvectors.col(j);
    const double value=v.dot(assembly.stiffness*v);
    // Roundoff may give a tiny negative value in the exact kernel.
    if(value < -1e-8)throw std::runtime_error("Materially negative curl energy");
    curl_norms[j]=std::sqrt(std::max(0.0,value));
  }
  for(const auto& tm:reference) {
    Eigen::VectorXd target=Eigen::VectorXd::Zero(ndof);
    double target_norm2=0;
    for_each_quadrature_point(integration,[&](Point2 p,double weight){
      const double sx=std::sin(tm.m*pi*(p.x-r.xmin)/a);
      const double sy=std::sin(tm.n*pi*(p.y-r.ymin)/b);
      const double f=sx*sy;
      target_norm2+=weight*f*f;
      const MLSResult v=mls.evaluate(p);
      for(std::size_t k=0;k<v.active_nodes.size();++k) {
        const auto idx=static_cast<Eigen::Index>(v.active_nodes[k]);
        target[idx] -= weight*f*v.dy[k];
        target[n+idx] += weight*f*v.dx[k];
      }
    });
    if(!(target_norm2>0))throw std::runtime_error("Invalid reference-mode norm");
    const double scalar_norm=std::sqrt(target_norm2);
    int nearest=-1,best=-1;
    double nearest_error=std::numeric_limits<double>::infinity();
    double strongest=-1;
    Eigen::VectorXd overlaps=Eigen::VectorXd::Zero(ndof);
    for(Eigen::Index i=0;i<ndof;++i) {
      if(spectrum.eigenvalues[i]<=spectrum.zero_cutoff)continue;
      const double err=relative_error(spectrum.eigenvalues[i],tm.lambda);
      const double ov=curl_norms[i] > 0 ?
           std::abs(target.dot(spectrum.eigenvectors.col(i)))/
            (curl_norms[i]*scalar_norm):0;
      if(!std::isfinite(ov))throw std::runtime_error("Nonfinite shape overlap");
      // Cauchy-Schwarz: guard against inconsistent quadrature.
      if(ov>1.000001)throw std::runtime_error("Shape overlap exceeds unity");
      overlaps[i]=std::min(1.0,ov);
      if(err<nearest_error) {
        nearest_error=err;nearest=static_cast<int>(i);
      }
      if(ov>strongest) {
        strongest=ov;best=static_cast<int>(i);
      }
    }
    TMShapeMatch out;
    out.reference=tm;
    if(nearest>=0) {
      out.nearest_index=nearest;
      out.nearest_lambda=spectrum.eigenvalues[nearest];
      out.nearest_relative_error=nearest_error;
      out.nearest_overlap=overlaps[nearest];
      double sumsq=0;
      for(Eigen::Index i=0;i<ndof;++i) {
        if(spectrum.eigenvalues[i]<=spectrum.zero_cutoff)continue;
        const double distance=std::abs(spectrum.eigenvalues[i]-out.nearest_lambda);
        const double scale=std::max({1.0,std::abs(spectrum.eigenvalues[i]),
                                    std::abs(out.nearest_lambda)});
        if(distance<=options.degeneracy_relative_tolerance*scale) {
          ++out.nearest_cluster_size;
          sumsq+=overlaps[i]*overlaps[i];
        }
      }
      out.nearest_cluster_overlap=std::min(1.0,std::sqrt(sumsq));
      out.algebraic_residual=spectrum.residuals[static_cast<std::size_t>(nearest)];
      const auto diagnostics=modal_diagnostics(
           mls,integration,spectrum.eigenvectors.col(nearest),out.nearest_lambda);
      out.divergence_ratio=diagnostics.divergence_ratio;
      out.normal_ratio=diagnostics.normal_ratio;
      out.supported_candidate=
        nearest_error<=options.maximum_relative_error &&
        out.nearest_cluster_overlap>=options.minimum_shape_overlap &&
        out.divergence_ratio<=options.maximum_divergence_ratio &&
        out.normal_ratio<=options.maximum_normal_ratio &&
        out.algebraic_residual<=options.maximum_algebraic_residual;
    }
    if(best>=0) {
      out.strongest_index=best;
      out.strongest_lambda=spectrum.eigenvalues[best];
      out.strongest_overlap=overlaps[best];
      out.strongest_relative_error=
        relative_error(out.strongest_lambda,tm.lambda);
    }
    matches.push_back(out);
  }
  return matches;
}
} // namespace al_maleh
