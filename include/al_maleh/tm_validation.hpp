// SPDX-License-Identifier: MIT
#pragma once
#include "al_maleh/spectral.hpp"
#include <vector>

namespace al_maleh {
struct TMShapeMatch {
  TMMode reference;
  int nearest_index=-1;
  double nearest_lambda=0;
  double nearest_relative_error=0;
  double nearest_overlap=0;
  double nearest_cluster_overlap=0; // projection to near-degenerate eigenspace
  int nearest_cluster_size=0;
  int strongest_index=-1;
  double strongest_lambda=0;
  double strongest_overlap=0;
  double strongest_relative_error=0;
  double divergence_ratio=0;
  double normal_ratio=0;
  double algebraic_residual=0;
  bool supported_candidate=false; // not a physical-mode certification
};
struct MatchOptions {
  double degeneracy_relative_tolerance=1e-6;
  double maximum_relative_error=0.15;
  double minimum_shape_overlap=0.75;
  double maximum_divergence_ratio=0.25;
  double maximum_normal_ratio=0.25;
  double maximum_algebraic_residual=1e-7;
};
// Fits analytic TM longitudinal Ez using curl(H_t) in L2.
// A subspace projection handles rotations within (near-)degenerate modes.
std::vector<TMShapeMatch> validate_tm_modes(
  const MLSEvaluator& mls,
  const BackgroundQuadrature& integration,
  const CartesianAssembly& assembly,
  const Spectrum& spectrum,
  int max_m,int max_n,
  const MatchOptions& options={});
} // namespace al_maleh
