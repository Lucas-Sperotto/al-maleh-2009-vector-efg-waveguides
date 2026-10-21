// SPDX-License-Identifier: MIT
#pragma once

#include "al_maleh/cartesian.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace al_maleh {

// Independent quadratic forms: no penalty is applied to the eigenproblem.
// D = integral_Omega div(N_i) div(N_j);
// B = integral_boundary (n.N_i)(n.N_j).
struct ConstraintAssembly {
  Eigen::SparseMatrix<double> divergence;
  Eigen::SparseMatrix<double> normal_boundary;
  std::size_t boundary_points = 0;
};

ConstraintAssembly assemble_constraints(
    const MLSEvaluator& mls, const BackgroundQuadrature& quadrature,
    int boundary_gauss_order = 3);

struct InvariantRange {
  double minimum = 0.0;
  double maximum = 0.0;
  // Coefficients c in the input basis V for the worst-case field Vc.
  Eigen::VectorXd maximizing_coefficients;
};

struct SubspaceConstraintMetrics {
  InvariantRange divergence_ratio;
  InvariantRange normal_ratio;
  Eigen::Index dimension = 0;
};

// For an input basis V (not necessarily M-orthonormal), compute extrema
// over all nonzero Vc of:
// sqrt((Vc)^T D (Vc) / (lambda_reference * (Vc)^T M (Vc)))
// sqrt((area/perimeter) (Vc)^T B (Vc) / ((Vc)^T M (Vc))).
// This is invariant to any nonsingular change of basis V -> V R.
// Individual eigenvectors may mix without altering these ranges.
SubspaceConstraintMetrics evaluate_subspace_constraints(
    const ConstraintAssembly& constraints,
    const Eigen::SparseMatrix<double>& mass,
    const Eigen::MatrixXd& subspace,
    double lambda_reference, Rectangle domain);

} // namespace al_maleh
