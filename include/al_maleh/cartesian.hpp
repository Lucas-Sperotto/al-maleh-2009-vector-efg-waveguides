// SPDX-License-Identifier: MIT
#pragma once

#include "al_maleh/mls.hpp"

#include <Eigen/SparseCore>
#include <cstddef>
#include <functional>
#include <vector>

namespace al_maleh {

struct Rectangle {
  double xmin = 0.0, xmax = 2.0, ymin = 0.0, ymax = 1.0;
};

struct BackgroundQuadrature {
  Rectangle domain;
  int cells_x = 8;
  int cells_y = 4;
  int gauss_order = 3; // 1..5 Gauss points per dimension (Al Maleh Ch. 5)
};

using QuadratureVisitor = std::function<void(Point2, double)>;

// A background mesh performs *integration only*: EFG shape functions
// are evaluated from the cloud; no finite-element basis is constructed.
void for_each_quadrature_point(const BackgroundQuadrature& rule,
                               const QuadratureVisitor& visitor);

std::vector<Node> make_regular_cloud(Rectangle rect, int nx, int ny,
                                     double support_radius);

struct CartesianAssembly {
  Eigen::SparseMatrix<double> stiffness;  // K: curl-curl
  Eigen::SparseMatrix<double> mass;       // M: vector L2
  std::size_t quadrature_points = 0;
  std::size_t max_active_nodes = 0;
  double max_moment_condition = 0.0;
};

// Vector field H_h=(sum phi_i Hx_i, sum phi_i Hy_i).
// DOF ordering: [Hx_0 ... Hx_(N-1), Hy_0 ... Hy_(N-1)].
// Material is uniform, dimensionless mu=epsilon=1 for this phase.
// PEC TM uses the *natural* curl(H)=0 boundary condition in this
// variational prototype; the huge gradient nullspace is NOT filtered.
CartesianAssembly assemble_cartesian(const MLSEvaluator& mls,
                                      const BackgroundQuadrature& rule);

} // namespace al_maleh
