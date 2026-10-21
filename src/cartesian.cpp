// SPDX-License-Identifier: MIT
#include "al_maleh/cartesian.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace al_maleh {
namespace {

void validate(Rectangle r) {
  if (!std::isfinite(r.xmin) || !std::isfinite(r.xmax) ||
      !std::isfinite(r.ymin) || !std::isfinite(r.ymax) ||
      !(r.xmin < r.xmax) || !(r.ymin < r.ymax)) {
    throw std::invalid_argument("Rectangle requires finite ordered limits");
  }
}

struct GaussRule1D {
  std::vector<double> points, weights;
};

GaussRule1D gauss_legendre(int order) {
  switch (order) {
    case 2:
      return {{-0.57735026918962576451, +0.57735026918962576451},
              {1.0, 1.0}};
    case 3:
      return {{-0.77459666924148337704, 0.0, +0.77459666924148337704},
              {5.0/9.0, 8.0/9.0, 5.0/9.0}};
    case 4:
      return {{-0.86113631159405257522, -0.33998104358485626480,
                0.33998104358485626480, 0.86113631159405257522},
              {0.34785484513745385737, 0.65214515486254614263,
               0.65214515486254614263, 0.34785484513745385737}};
    default:
      throw std::invalid_argument("Quadrature order must be 2, 3, or 4");
  }
}

} // namespace

void for_each_quadrature_point(const BackgroundQuadrature& rule,
                               const QuadratureVisitor& visitor) {
  validate(rule.domain);
  if (rule.cells_x < 1 || rule.cells_y < 1)
    throw std::invalid_argument("Quadrature needs positive cell counts");
  if (!visitor)
    throw std::invalid_argument("Quadrature visitor is empty");
  const GaussRule1D g = gauss_legendre(rule.gauss_order);
  const double hx = (rule.domain.xmax - rule.domain.xmin) / rule.cells_x;
  const double hy = (rule.domain.ymax - rule.domain.ymin) / rule.cells_y;
  const double jacobian = hx * hy / 4.0;

  for (int cy = 0; cy < rule.cells_y; ++cy) {
    const double y0 = rule.domain.ymin + cy * hy;
    for (int cx = 0; cx < rule.cells_x; ++cx) {
      const double x0 = rule.domain.xmin + cx * hx;
      for (std::size_t iy = 0; iy < g.points.size(); ++iy)
        for (std::size_t ix = 0; ix < g.points.size(); ++ix) {
          const Point2 p{x0 + 0.5 * hx * (1.0 + g.points[ix]),
                         y0 + 0.5 * hy * (1.0 + g.points[iy])};
          const double weight = jacobian * g.weights[ix] * g.weights[iy];
          visitor(p, weight);
        }
    }
  }
}

std::vector<Node> make_regular_cloud(Rectangle domain, int nx, int ny,
                                     double support_radius) {
  validate(domain);
  if (nx < 2 || ny < 2 || !std::isfinite(support_radius) ||
      support_radius <= 0.0)
    throw std::invalid_argument("Invalid regular cloud parameters");
  std::vector<Node> nodes;
  nodes.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      nodes.push_back({{domain.xmin +
                          (domain.xmax-domain.xmin) * (double(i)/(nx-1)),
                         domain.ymin +
                          (domain.ymax-domain.ymin) * (double(j)/(ny-1))},
                        support_radius});
  return nodes;
}

CartesianAssembly assemble_cartesian(const MLSEvaluator& mls,
                                      const BackgroundQuadrature& rule) {
  // M2/M4 reference prototype: accumulate directly into compact dense arrays.
  // The previous code emitted SIX sparse triplets for each active-node pair
  // at EVERY quadrature point. At a 17x15 cloud, this can allocate many
  // millions of temporary triplets even though the final K/M fit in MB.
  //
  // This is NOT yet a scalable sparse assembler; it is a memory-conscious
  // assembly path for small-to-medium clouds and a dense reference solver.
  const std::size_t n = mls.nodes().size();
  if (n > 600)
    throw std::invalid_argument(
        "Dense assembly workspace supports at most 600 nodes");
  if (n > static_cast<std::size_t>(std::numeric_limits<int>::max()/2))
    throw std::invalid_argument("Too many nodes for current DOF indexing");
  const int ndof = static_cast<int>(2*n);
  Eigen::MatrixXd denseK = Eigen::MatrixXd::Zero(ndof, ndof);
  Eigen::MatrixXd denseM = Eigen::MatrixXd::Zero(ndof, ndof);
  CartesianAssembly result;

  for_each_quadrature_point(rule, [&](Point2 p, double w) {
    const MLSResult v = mls.evaluate(p);
    const std::size_t m = v.active_nodes.size();
    if (v.phi.size()!=m || v.dx.size()!=m || v.dy.size()!=m)
      throw std::runtime_error("Inconsistent MLS support vectors");
    ++result.quadrature_points;
    result.max_active_nodes = std::max(result.max_active_nodes, m);
    result.max_moment_condition =
        std::max(result.max_moment_condition, v.condition_number);

    // K_ij = integral curl(N_i) curl(N_j)
    // curl([phi_i,0])=-phi_i,y; curl([0,phi_i])=phi_i,x.
    // M_ij = integral N_i dot N_j.
    for (std::size_t a = 0; a < m; ++a) {
      const int i = static_cast<int>(v.active_nodes[a]);
      const double phi_i = v.phi[a], dx_i = v.dx[a], dy_i = v.dy[a];
      for (std::size_t b = 0; b < m; ++b) {
        const int j = static_cast<int>(v.active_nodes[b]);
        denseK(i,   j)   += w * dy_i * v.dy[b];
        denseK(i,   j+static_cast<int>(n)) += -w * dy_i * v.dx[b];
        denseK(i+static_cast<int>(n), j)   += -w * dx_i * v.dy[b];
        denseK(i+static_cast<int>(n),
               j+static_cast<int>(n)) += w * dx_i * v.dx[b];
        const double mass = w * phi_i * v.phi[b];
        denseM(i,j) += mass;
        denseM(i+static_cast<int>(n),j+static_cast<int>(n)) += mass;
      }
    }
  });

  result.stiffness = denseK.sparseView(0.0, 0.0);
  result.mass = denseM.sparseView(0.0, 0.0);
  result.stiffness.makeCompressed();
  result.mass.makeCompressed();
  return result;
}

} // namespace al_maleh
