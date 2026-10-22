// SPDX-License-Identifier: MIT
#pragma once
#include "al_maleh/cartesian.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <vector>

namespace al_maleh {
// Original Al Maleh (2009), Chapter 6, Eqs. 6-1 .. 6-6.
// Internal coefficient ordering: [all gradient coefficients; all J-gradient
// coefficients]. Author's interleaved listing differs by permutation only.
struct RotIrrotAssembly {
  Eigen::MatrixXd stiffness, mass;
  std::size_t quadrature_points = 0;
  double max_mls_condition = 0;
};
RotIrrotAssembly assemble_rot_irrot(const MLSEvaluator& mls,
                                    const BackgroundQuadrature& rule);

struct RankReducedSpectrum {
  Eigen::VectorXd eigenvalues;      // all retained EVP values incl. negative
  Eigen::MatrixXd eigenvectors;     // 2*N x rank, M-orthonormal
  std::vector<double> reduced_residuals, full_residuals;
  int retained_rank=0, removed_mass_nullity=0, numerical_curl_kernel=0;
  double mass_eigen_max=0, mass_eigen_min_retained=0, mass_rank_cutoff=0;
  double mass_eigen_max_removed=0;  // largest discarded mass eigenvalue (signed)
  double eigenvalue_zero_cutoff=0, M_orthogonality_error=0;
  // ||K Z||_F/||K||_F for the orthonormal discarded mass directions Z.
  // ~eps certifies ker(M_removed) in ker(K): the projected pencil then has
  // exactly the finite spectrum of the original singular pencil (K,M).
  double mass_kernel_stiffness_coupling=0;
  double max_reduced_residual=0, max_full_residual=0;
};
// A singular M is intrinsic to gradient/J-gradient MLS bases. Project onto
// the positive-mass subspace with an explicit spectral cutoff; report removed
// dimensions, never synthesize mass eigenvalues or drop negative λ silently.
// The projected residual Q^T(Kh-λMh) is ~0 by construction and does NOT
// certify the reduction: use max_full_residual and the coupling above.
RankReducedSpectrum solve_rot_irrot_rank_reduced(
    const RotIrrotAssembly& forms,
    double relative_mass_rank_cutoff=1e-10,
    double absolute_mass_rank_cutoff=1e-12);

struct RotIrrotField {
  double hx=0,hy=0;
  double irrot_x=0,irrot_y=0;
  double rot_x=0,rot_y=0;
  double divergence=0,curl=0;
};
RotIrrotField evaluate_rot_irrot(const MLSResult& shape,
    const Eigen::VectorXd& coefficients,Eigen::Index node_count);
} // namespace al_maleh
