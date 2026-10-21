// SPDX-License-Identifier: MIT
#pragma once

#include "al_maleh/cartesian.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <vector>

namespace al_maleh {

// Independent SCALAR MLS/EFG control for the TM longitudinal potential Ez.
// This is NOT the Cartesian vector EFG formulation and does NOT fix its
// spectral pollution. Boundary values are enforced only at cloud boundary
// collocation points (not strongly along the entire continuum boundary).
struct ScalarControlSystem {
  Eigen::MatrixXd stiffness;     // int grad(phi_i).grad(phi_j)
  Eigen::MatrixXd mass;          // int phi_i phi_j
  Eigen::MatrixXd boundary;      // C(i,j) = phi_j(boundary_point_i)
  Eigen::MatrixXd nullspace;     // columns Q spanning ker(C)
  Eigen::MatrixXd reduced_stiffness, reduced_mass;
  std::size_t quadrature_points = 0;
  std::size_t boundary_points = 0;
  Eigen::Index boundary_rank = 0;
  double max_moment_condition = 0.0;
  double boundary_constraint_residual = 0.0;
};

ScalarControlSystem assemble_scalar_control(
    const MLSEvaluator& mls, const BackgroundQuadrature& quadrature);

struct ScalarControlSpectrum {
  Eigen::VectorXd eigenvalues;   // lambda = kc^2, ascending
  Eigen::MatrixXd coefficients; // full scalar MLS coefficient vectors
  std::vector<double> projected_residuals;
  std::vector<double> collocation_residuals;
  double mass_condition = 0.0;
  double mass_orthogonality_error = 0.0;
};

ScalarControlSpectrum solve_scalar_control(const ScalarControlSystem& system);

// Exact scalar field reconstruction from the same MLS coefficients.
// H = (-partial_y Ez, partial_x Ez) is analytically solenoidal.
// The collocation constraint Cc=0 does not imply zero trace between nodes.
struct ScalarControlFieldMetrics {
  double relative_shape_overlap = 0.0;
  double boundary_trace_ratio = 0.0;
  double normal_flux_ratio = 0.0;
  double divergence_l2 = 0.0;
};

ScalarControlFieldMetrics evaluate_scalar_control_mode(
    const MLSEvaluator& mls, const BackgroundQuadrature& quadrature,
    const Eigen::VectorXd& coefficients, int m, int n,
    int boundary_gauss_order = 3);

} // namespace al_maleh
