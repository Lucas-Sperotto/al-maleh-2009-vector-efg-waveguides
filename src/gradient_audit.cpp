// SPDX-License-Identifier: MIT
#include "al_maleh/gradient_audit.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace al_maleh {
GradientPairing assemble_gradient_pairing(
    const MLSEvaluator& mls, const BackgroundQuadrature& quadrature) {
  const Eigen::Index n = static_cast<Eigen::Index>(mls.nodes().size());
  if(n<3 || n>600) throw std::invalid_argument("M6d supports 3..600 scalar nodes");
  GradientPairing out;
  out.coupling = Eigen::MatrixXd::Zero(2*n,n);
  out.scalar_energy = Eigen::MatrixXd::Zero(n,n);

  for_each_quadrature_point(quadrature,[&](Point2 p,double weight) {
    const MLSResult s = mls.evaluate(p);
    const auto count=s.active_nodes.size();
    if(s.phi.size()!=count||s.dx.size()!=count||s.dy.size()!=count)
      throw std::runtime_error("Inconsistent MLS basis vectors");
    ++out.quadrature_points;
    for(std::size_t a=0;a<count;++a) {
      const Eigen::Index i=static_cast<Eigen::Index>(s.active_nodes[a]);
      for(std::size_t b=0;b<count;++b) {
        const Eigen::Index j=static_cast<Eigen::Index>(s.active_nodes[b]);
        out.coupling(i,j)   += weight*s.phi[a]*s.dx[b];
        out.coupling(n+i,j) += weight*s.phi[a]*s.dy[b];
        out.scalar_energy(i,j) +=
            weight*(s.dx[a]*s.dx[b]+s.dy[a]*s.dy[b]);
      }
    }
  });
  if(!out.coupling.allFinite()||!out.scalar_energy.allFinite())
    throw std::runtime_error("Nonfinite gradient/MLS volume pairing");
  return out;
}

GradientProjection build_gradient_projection(
    const Eigen::SparseMatrix<double>& vector_mass,
    const GradientPairing& pairing,
    const GradientAuditOptions& options) {
  const Eigen::Index n=pairing.scalar_energy.rows(),d=2*n;
  if(n<1||pairing.scalar_energy.cols()!=n||
     pairing.coupling.rows()!=d||pairing.coupling.cols()!=n||
     vector_mass.rows()!=d||vector_mass.cols()!=d||
     !pairing.scalar_energy.allFinite()||!pairing.coupling.allFinite()||
     !std::isfinite(options.relative_rank_cutoff)||
     options.relative_rank_cutoff<=0||options.relative_rank_cutoff>=1||
     !std::isfinite(options.maximum_mass_condition)||
     options.maximum_mass_condition<=1)
    throw std::invalid_argument("Invalid M6d gradient projection input");
  Eigen::MatrixXd M(vector_mass);
  if(!M.allFinite()||
     (M-M.transpose()).norm()>1e-10*std::max(1.0,M.norm()))
    throw std::invalid_argument("M6d mass must be finite and symmetric");
  M=(0.5*(M+M.transpose())).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ms(M,Eigen::EigenvaluesOnly);
  if(ms.info()!=Eigen::Success||!(ms.eigenvalues().minCoeff()>0))
    throw std::runtime_error("M6d vector mass must be positive definite");
  GradientProjection result;
  result.mass_condition=ms.eigenvalues().maxCoeff()/ms.eigenvalues().minCoeff();
  if(!std::isfinite(result.mass_condition)||
     result.mass_condition>options.maximum_mass_condition)
    throw std::runtime_error("M6d vector mass too ill-conditioned");
  Eigen::LLT<Eigen::MatrixXd> chol(M);
  if(chol.info()!=Eigen::Success)
    throw std::runtime_error("M6d mass Cholesky failed");
  result.coefficients=chol.solve(pairing.coupling);
  Eigen::MatrixXd projected=pairing.coupling.transpose()*result.coefficients;
  projected=(0.5*(projected+projected.transpose())).eval();
  if(!projected.allFinite())throw std::runtime_error("Nonfinite projected scalar Gram");
  result.projected_scalar_gram=projected;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eg(projected);
  if(eg.info()!=Eigen::Success)
    throw std::runtime_error("M6d gradient Gram eigensolver failed");
  result.projected_scalar_eigenvalues=eg.eigenvalues();
  const double largest=eg.eigenvalues().maxCoeff();
  if(!(largest>0))throw std::runtime_error("Projected gradient Gram has no positive directions");
  if(eg.eigenvalues().minCoeff() < -1e-9*largest)
    throw std::runtime_error("Projected gradient Gram is not positive semidefinite");
  result.rank_cutoff=options.relative_rank_cutoff*largest;
  const Eigen::Index first=std::distance(
    eg.eigenvalues().data(),
    std::lower_bound(eg.eigenvalues().data(),
                     eg.eigenvalues().data()+n,result.rank_cutoff));
  result.numerical_rank=n-first;
  if(result.numerical_rank<1||result.numerical_rank>=d)
    throw std::runtime_error("Invalid numerical rank for projected MLS gradients");

  // q_i = M^-1 H u_i / sqrt(sigma_i), for sigma_i > rank_cutoff.
  // Thus q_i^T M q_j = delta_ij and P_g=Q Q^T M.
  result.basis=result.coefficients*
      eg.eigenvectors().rightCols(result.numerical_rank);
  for(Eigen::Index i=0;i<result.numerical_rank;++i)
    result.basis.col(i)/=std::sqrt(eg.eigenvalues()[first+i]);
  const Eigen::MatrixXd gram=result.basis.transpose()*M*result.basis;
  result.maximum_m_orthogonality_error=
      (gram-Eigen::MatrixXd::Identity(result.numerical_rank,result.numerical_rank))
      .cwiseAbs().maxCoeff();
  if(!result.basis.allFinite()||
     result.maximum_m_orthogonality_error>5e-6)
    throw std::runtime_error("M6d projected basis is not M-orthonormal");

  // The exact energy ||grad phi(c)||² = c^T S c, whereas the best
  // Cartesian vector representation has energy c^T G c.
  const double exact=pairing.scalar_energy.norm();
  result.relative_gradient_reconstruction_frobenius=
      (pairing.scalar_energy-projected).norm()/
      std::max(1e-300,exact);
  if(!std::isfinite(result.relative_gradient_reconstruction_frobenius))
    throw std::runtime_error("Invalid gradient approximation discrepancy");
  return result;
}

GradientModeMetrics evaluate_gradient_mode(
    const GradientProjection& projection,
    const Eigen::SparseMatrix<double>& vector_mass,
    const Eigen::SparseMatrix<double>& vector_stiffness,
    const Eigen::VectorXd& v) {
  const Eigen::Index d=vector_mass.rows();
  if(d<1||vector_mass.cols()!=d||
     vector_stiffness.rows()!=d||vector_stiffness.cols()!=d||
     v.size()!=d||!v.allFinite()||projection.basis.rows()!=d||
     projection.basis.cols()!=projection.numerical_rank||
     projection.numerical_rank<1)
    throw std::invalid_argument("Invalid M6d vector mode or projected basis");
  const double norm2=v.dot(vector_mass*v);
  if(!(norm2>0)||!std::isfinite(norm2))
    throw std::runtime_error("Invalid vector mass energy");
  const Eigen::VectorXd alpha=projection.basis.transpose()*(vector_mass*v);
  const Eigen::VectorXd grad=projection.basis*alpha;
  const double projected2=grad.dot(vector_mass*grad);
  const double frac=projected2/norm2;
  if(!std::isfinite(frac)||frac < -1e-8||frac > 1+1e-6)
    throw std::runtime_error("Gradient projection violates energy bound");
  GradientModeMetrics result;
  result.gradient_energy_fraction=std::clamp(frac,0.0,1.0);
  result.complement_energy_fraction=1-result.gradient_energy_fraction;
  result.projected_vector_norm=std::sqrt(std::max(0.0,projected2));
  if(projected2>1e-20*norm2) {
    const double k2=grad.dot(vector_stiffness*grad);
    if(!std::isfinite(k2)||k2 < -1e-8*std::max(1.,std::abs(k2)))
      throw std::runtime_error("Negative projected-gradient curl energy");
    result.projected_gradient_curl_rayleigh=std::max(0.0,k2)/projected2;
  }
  return result;
}

double scalar_gradient_fit_relative_error(
    const GradientPairing& pairing,
    const GradientProjection& projection,
    const Eigen::VectorXd& c) {
  if(c.size()!=pairing.scalar_energy.rows()||!c.allFinite()||
     projection.projected_scalar_gram.rows()!=c.size())
    throw std::invalid_argument("Invalid scalar gradient coefficients");
  const double exact=c.dot(pairing.scalar_energy*c);
  if(!(exact>1e-20)||!std::isfinite(exact))
    throw std::invalid_argument("Scalar gradient has zero or invalid energy");
  const double captured=c.dot(projection.projected_scalar_gram*c);
  const double frac=captured/exact;
  if(!std::isfinite(frac)||frac < -1e-8||frac>1+1e-6)
    throw std::runtime_error("L2 gradient fit energy inequality violated");
  return std::sqrt(std::max(0.0,1-std::clamp(frac,0.0,1.0)));
}

} // namespace al_maleh
