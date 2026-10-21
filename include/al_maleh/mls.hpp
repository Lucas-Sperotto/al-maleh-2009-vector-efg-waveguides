// SPDX-License-Identifier: MIT
// Core MLS shapes for an independent reconstruction of Al Maleh (2009).
// Mathematical approach adapted from Lucas-Sperotto/efgmi-hcurl-diffuse-whitney.
#pragma once

#include <cstddef>
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

// A fixed, centered/scaled polynomial basis is used: degree 1 or degree 2.
// Wendland C2 compact supports with per-node radii; no nodal interpolation claim.
class MLSEvaluator {
 public:
  explicit MLSEvaluator(std::vector<Node> nodes,
                        int polynomial_degree = 1,
                        double max_condition = 1.0e12);

  MLSResult evaluate(Point2 point) const;
  int polynomial_degree() const noexcept { return degree_; }
  const std::vector<Node>& nodes() const noexcept { return nodes_; }

 private:
  std::vector<Node> nodes_;
  int degree_ = 1;
  double origin_x_ = 0.0;
  double origin_y_ = 0.0;
  double scale_ = 1.0;
  double max_condition_ = 1.0e12;
};

}  // namespace al_maleh
