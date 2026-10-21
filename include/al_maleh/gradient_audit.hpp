// SPDX-License-Identifier: MIT
#pragma once
#include "al_maleh/cartesian.hpp"
#include "al_maleh/spectral.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace al_maleh {

// Discrete L2 pairing between the Cartesian vector MLS basis N_i and
// gradients of the *scalar* MLS basis phi_j. These finite-dimensional
// spaces are NOT a verified compatible discrete de Rham sequence.
struct GradientPairing {
  Eigen::MatrixXd coupling;     // H_{ij} = integral N_i . grad(phi_j)
  Eigen::MatrixXd scalar_energy;// S_{ij} = integral grad(phi_i).grad(phi_j)
  std::size_t quadrature_points = 0;
};

// Uses the same MLS evaluator and volume quadrature as the original EFG K/M.
// No constrained eigenproblem, penalty or modification of K/M.
GradientPairing assemble_gradient_pairing(
    const MLSEvaluator& mls, const BackgroundQuadrature& quadrature);

struct GradientAuditOptions {
  // Rank of the M-projected gradient space is numerically defined, not
  // assumed equal to scalar_nodes-1.
  double relative_rank_cutoff = 1e-8;
  double maximum_mass_condition = 1e13;
};

struct GradientProjection {
  Eigen::MatrixXd basis;       // Q: M-orthonormal Cartesian vector columns
  Eigen::MatrixXd coefficients;// M^-1 H: best L2 vector fit to every grad phi_j
  Eigen::MatrixXd projected_scalar_gram; // H^T M^-1 H
  Eigen::VectorXd projected_scalar_eigenvalues; // ascending, entire spectrum
  Eigen::Index numerical_rank = 0;
  double rank_cutoff = 0.0;
  double mass_condition = 0.0;
  double maximum_m_orthogonality_error = 0.0;
  double relative_gradient_reconstruction_frobenius = 0.0;
};

GradientProjection build_gradient_projection(
    const Eigen::SparseMatrix<double>& vector_mass,
    const GradientPairing& pairing,
    const GradientAuditOptions& options = {});

struct GradientModeMetrics {
  double gradient_energy_fraction = 0.0; // ||P_g v||_M² / ||v||_M²
  double complement_energy_fraction = 0.0;
  double projected_gradient_curl_rayleigh = 0.0;
  double projected_vector_norm = 0.0; // ||P_g v||_M
};

// For a vector v, P_g v = Q Q^T M v and Q^T M Q = I. The split
// is a DIAGNOSTIC: it is not an invariant physical Helmholtz
// decomposition unless the discrete spaces are compatible.
GradientModeMetrics evaluate_gradient_mode(
    const GradientProjection& projection,
    const Eigen::SparseMatrix<double>& vector_mass,
    const Eigen::SparseMatrix<double>& vector_stiffness,
    const Eigen::VectorXd& mode);

// Relative best-fit error in L2 for grad(sum_j scalar_coeff_j phi_j).
// Compares exact integrated scalar gradient energy c^T S c to the
// energy captured by the vector MLS subspace c^T H^T M^-1 H c.
double scalar_gradient_fit_relative_error(
    const GradientPairing& pairing,
    const GradientProjection& projection,
    const Eigen::VectorXd& scalar_coefficients);

} // namespace al_maleh
