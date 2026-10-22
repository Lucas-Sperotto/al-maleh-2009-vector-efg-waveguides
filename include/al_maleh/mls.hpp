// SPDX-License-Identifier: MIT
// Core MLS shapes for an independent reconstruction of Al Maleh (2009).
// Mathematical approach adapted from Lucas-Sperotto/efgmi-hcurl-diffuse-whitney.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace al_maleh {

struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

struct Node {
  Point2 position;
  double support_radius = 1.0;
};

struct MLSResult {
  std::vector<std::size_t> active_nodes;
  std::vector<double> phi, dx, dy, dxx, dxy, dyy;
  double condition_number = 0.0;
};

// Radial weight w(q), q=|x-x_i|/r_i, zero for q>=1. Both profiles are C2.
// WendlandC2:    (1-q)^4 (4q+1)            -- engineering choice of this
//                project (historical default, NOT taken from the thesis).
// QuarticSpline: 1-6q^2+8q^3-3q^4          -- Al Maleh (2009), Eq. (4-2).
// The thesis does not state the numerical value of the radius r.
enum class MLSWeight { WendlandC2, QuarticSpline };
const char* mls_weight_name(MLSWeight weight) noexcept;
MLSWeight parse_mls_weight(const std::string& name);

// Polynomial basis of degree 1 or 2, shifted to each query point and scaled
// by the largest active support radius (affine-invariant; condition_number is then
// a geometric measure of the local node distribution, not of cloud size).
// Compact radial supports with per-node radii; no nodal interpolation claim.
class MLSEvaluator {
 public:
  explicit MLSEvaluator(std::vector<Node> nodes,
                        int polynomial_degree = 1,
                        double max_condition = 1.0e12,
                        MLSWeight weight = MLSWeight::WendlandC2);

  MLSResult evaluate(Point2 point) const;
  int polynomial_degree() const noexcept { return degree_; }
  MLSWeight weight() const noexcept { return weight_; }
  const std::vector<Node>& nodes() const noexcept { return nodes_; }

 private:
  std::vector<Node> nodes_;
  int degree_ = 1;
  double max_condition_ = 1.0e12;
  MLSWeight weight_ = MLSWeight::WendlandC2;
};

}  // namespace al_maleh
