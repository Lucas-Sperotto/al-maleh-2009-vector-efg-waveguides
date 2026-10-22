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
  // removed_mass_nullity counts discarded directions at the chosen cutoff;
  // it need not equal the exact nullity if positive mass modes are discarded.
  int retained_rank=0, removed_mass_nullity=0, numerical_curl_kernel=0;
  double mass_eigen_max=0, mass_eigen_min_retained=0, mass_rank_cutoff=0;
  double mass_eigen_max_removed=0;  // largest discarded mass eigenvalue (signed)
  double eigenvalue_zero_cutoff=0, M_orthogonality_error=0;
  // ||K Z||_F/||K||_F for the orthonormal discarded mass directions Z.
  // Small coupling and full residuals certify the returned eigenpairs, but
  // completeness also requires the discarded space to be exactly ker(M).
  // Otherwise even KZ=0 can hide a discarded finite zero eigenvalue. When
  // Z spans ker(M) and KZ=0, projection gives the spectrum on the quotient
  // by this common kernel. Check the mass gap and cutoff stability as well.
  double mass_kernel_stiffness_coupling=0;
  double max_reduced_residual=0, max_full_residual=0;
};
// A singular M is intrinsic to gradient/J-gradient MLS bases. Project onto
// the positive-mass subspace with an explicit spectral cutoff; report removed
// dimensions, never synthesize mass eigenvalues or drop negative λ silently.
// The projected residual Q^T(Kh-λMh) is ~0 by construction and does NOT
// certify the reduction: use max_full_residual, the coupling above and the
// discarded/retained mass gap. reduced_residuals use Q^T times the full
// residual, with the same full-pencil normalization as full_residuals.
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
