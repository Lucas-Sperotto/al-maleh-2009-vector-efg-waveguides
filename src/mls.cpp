// SPDX-License-Identifier: MIT
// Independently adapted from consistent first-derivative MLS in
// Lucas-Sperotto/efgmi-hcurl-diffuse-whitney (2026, MIT).
// Adds consistent second derivatives, Wendland C2 supports and the
// quartic spline weight of Al Maleh (2009), Eq. (4-2).
#include "al_maleh/mls.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace al_maleh {
namespace {

struct Basis {
  Eigen::VectorXd p, px, py, pxx, pxy, pyy;
};

Basis polynomial(Point2 x, int degree, double cx, double cy, double scale) {
  const int m = degree == 1 ? 3 : 6;
  Basis b{Eigen::VectorXd::Zero(m), Eigen::VectorXd::Zero(m),
          Eigen::VectorXd::Zero(m), Eigen::VectorXd::Zero(m),
          Eigen::VectorXd::Zero(m), Eigen::VectorXd::Zero(m)};
  const double u = (x.x - cx) / scale;
  const double v = (x.y - cy) / scale;
  const double s = 1.0 / scale;
  b.p(0) = 1.0;
  b.p(1) = u;
  b.p(2) = v;
  b.px(1) = s;
  b.py(2) = s;
  if (degree == 2) {
    b.p(3) = u * u;
    b.p(4) = u * v;
    b.p(5) = v * v;
    b.px(3) = 2.0 * u * s;
    b.px(4) = v * s;
    b.py(4) = u * s;
    b.py(5) = 2.0 * v * s;
    b.pxx(3) = 2.0 * s * s;
    b.pxy(4) = s * s;
    b.pyy(5) = 2.0 * s * s;
  }
  return b;
}

struct Weight {
  double w, wx, wy, wxx, wxy, wyy;
};

Weight wendland_c2(Point2 x, const Node& node) {
  const double dx = x.x - node.position.x;
  const double dy = x.y - node.position.y;
  const double r = std::hypot(dx, dy);
  const double h = node.support_radius;
  const double q = r / h;
  if (q >= 1.0) return {0, 0, 0, 0, 0, 0};
  const double t = 1.0 - q;
  const double t2 = t * t;
  const double t3 = t2 * t;
  const double factor = -20.0 * t3 / (h * h);
  // Hessian remains continuous at the center (r -> 0).
  const double radial = r == 0.0 ? 0.0 : 60.0 * t2 / (h * h * h * r);
  return {t2 * t2 * (4.0 * q + 1.0),
          factor * dx, factor * dy,
          factor + radial * dx * dx,
          radial * dx * dy,
          factor + radial * dy * dy};
}

// Al Maleh (2009), Eq. (4-2): w = 1-6q^2+8q^3-3q^4, dw/dq=-12q(1-q)^2.
// grad w = -12(1-q)^2/h^2 (x-x_i); the radial Hessian term vanishes at r=0.
Weight quartic_spline(Point2 x, const Node& node) {
  const double dx = x.x - node.position.x;
  const double dy = x.y - node.position.y;
  const double r = std::hypot(dx, dy);
  const double h = node.support_radius;
  const double q = r / h;
  if (q >= 1.0) return {0, 0, 0, 0, 0, 0};
  const double t = 1.0 - q;
  const double factor = -12.0 * t * t / (h * h);
  const double radial = r == 0.0 ? 0.0 : 24.0 * t / (h * h * h * r);
  // Factorization avoids cancellation (and spurious negative weights) near
  // q=1, where the expanded polynomial subtracts numbers of order one.
  return {t * t * t * (1.0 + 3.0 * q),
          factor * dx, factor * dy,
          factor + radial * dx * dx,
          radial * dx * dy,
          factor + radial * dy * dy};
}

struct Active {
  std::size_t index;
  Eigen::VectorXd p;
  Weight w;
};

bool finite_point(Point2 p) {
  return std::isfinite(p.x) && std::isfinite(p.y);
}

}  // namespace

const char* mls_weight_name(MLSWeight weight) noexcept {
  return weight == MLSWeight::QuarticSpline ? "quartic_spline" : "wendland_c2";
}

MLSWeight parse_mls_weight(const std::string& name) {
  if (name == "wendland_c2") return MLSWeight::WendlandC2;
  if (name == "quartic_spline") return MLSWeight::QuarticSpline;
  throw std::invalid_argument("MLS weight must be wendland_c2 or quartic_spline");
}

MLSEvaluator::MLSEvaluator(std::vector<Node> nodes, int degree,
                           double max_condition, MLSWeight weight)
    : nodes_(std::move(nodes)), degree_(degree),
      max_condition_(max_condition), weight_(weight) {
  if (weight_ != MLSWeight::WendlandC2 && weight_ != MLSWeight::QuarticSpline)
    throw std::invalid_argument("Unknown MLS weight function");
  if (degree_ != 1 && degree_ != 2)
    throw std::invalid_argument("MLS degree must be 1 or 2");
  if (!std::isfinite(max_condition_) || max_condition_ < 1.0)
    throw std::invalid_argument("Invalid moment condition limit");
  if (nodes_.size() < static_cast<std::size_t>(degree_ == 1 ? 3 : 6))
    throw std::invalid_argument("Insufficient MLS cloud size");
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = xmin;
  double xmax = -xmin;
  double ymax = -xmin;
  for (const Node& n : nodes_) {
    if (!finite_point(n.position) || !std::isfinite(n.support_radius) ||
        n.support_radius <= 0.0)
      throw std::invalid_argument("Invalid MLS node or support radius");
    xmin = std::min(xmin, n.position.x);
    xmax = std::max(xmax, n.position.x);
    ymin = std::min(ymin, n.position.y);
    ymax = std::max(ymax, n.position.y);
  }
  if (!(std::max(xmax - xmin, ymax - ymin) > 0.0))
    throw std::invalid_argument("Degenerate MLS cloud bounding box");
}

MLSResult MLSEvaluator::evaluate(Point2 x) const {
  if (!finite_point(x)) throw std::invalid_argument("Non-finite MLS query point");
  const int m = degree_ == 1 ? 3 : 6;
  using Matrix = Eigen::MatrixXd;
  Matrix A = Matrix::Zero(m, m), Ax = A, Ay = A;
  Matrix Axx = A, Axy = A, Ayy = A;
  std::vector<Active> active;
  active.reserve(nodes_.size());

  // Gather supports first: an inactive node must not affect normalization or
  // conditioning of the local fit, even if it has a much larger radius.
  double scale = 0.0;
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    const Weight w = weight_ == MLSWeight::QuarticSpline
                         ? quartic_spline(x, nodes_[i])
                         : wendland_c2(x, nodes_[i]);
    if (w.w <= 0.0) continue;
    scale = std::max(scale, nodes_[i].support_radius);
    active.push_back({i, {}, w});
  }
  if (active.size() < static_cast<std::size_t>(m))
    throw std::runtime_error("MLS support has fewer active nodes than basis dimension");

  // The basis is shifted to the query point and scaled by the largest active
  // support radius. MLS shapes are invariant under this affine change of a
  // fixed polynomial basis, so values and derivatives are unchanged in exact
  // arithmetic, while cond(A) no longer grows like (cloud size/h)^(2*degree).
  // Freeze this coordinate transformation at the current query while taking
  // derivatives: pi is then constant, and b below has the derivatives of the
  // evaluation polynomial in that fixed basis. Re-centering at the next query
  // is an equivalent change of basis, not a term to add to Ax or Axx.
  for (Active& a : active) {
    a.p = polynomial(nodes_[a.index].position, degree_, x.x, x.y, scale).p;
    const Eigen::VectorXd& pi = a.p;
    const Weight& w = a.w;
    const Matrix outer = pi * pi.transpose();
    A.noalias() += w.w * outer;
    Ax.noalias() += w.wx * outer;
    Ay.noalias() += w.wy * outer;
    Axx.noalias() += w.wxx * outer;
    Axy.noalias() += w.wxy * outer;
    Ayy.noalias() += w.wyy * outer;
  }

  Eigen::SelfAdjointEigenSolver<Matrix> spectrum(A, Eigen::EigenvaluesOnly);
  if (spectrum.info() != Eigen::Success)
    throw std::runtime_error("MLS moment eigensolver failed");
  const double amin = spectrum.eigenvalues().minCoeff();
  const double amax = spectrum.eigenvalues().maxCoeff();
  if (!(amin > 0.0) || !std::isfinite(amin) || !std::isfinite(amax))
    throw std::runtime_error("MLS moment matrix is not positive definite");
  const double condition = amax / amin;
  if (!std::isfinite(condition) || condition > max_condition_)
    throw std::runtime_error("MLS moment matrix is ill-conditioned");

  Eigen::LLT<Matrix> factor(A);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error("MLS Cholesky factorization failed");
  const Basis b = polynomial(x, degree_, x.x, x.y, scale);
  // Solve A gamma = p, then differentiate A gamma = p twice.
  const Eigen::VectorXd g = factor.solve(b.p);
  const Eigen::VectorXd gx = factor.solve(b.px - Ax * g);
  const Eigen::VectorXd gy = factor.solve(b.py - Ay * g);
  const Eigen::VectorXd gxx = factor.solve(b.pxx - Axx * g - 2.0 * Ax * gx);
  const Eigen::VectorXd gxy = factor.solve(b.pxy - Axy * g - Ax * gy - Ay * gx);
  const Eigen::VectorXd gyy = factor.solve(b.pyy - Ayy * g - 2.0 * Ay * gy);
  if (!g.allFinite() || !gx.allFinite() || !gy.allFinite() ||
      !gxx.allFinite() || !gxy.allFinite() || !gyy.allFinite())
    throw std::runtime_error("Non-finite MLS coefficient derivatives");

  MLSResult result;
  result.condition_number = condition;
  for (const Active& a : active) {
    const double z = a.p.dot(g), zx = a.p.dot(gx), zy = a.p.dot(gy);
    const Weight& w = a.w;
    result.active_nodes.push_back(a.index);
    result.phi.push_back(w.w * z);
    result.dx.push_back(w.wx * z + w.w * zx);
    result.dy.push_back(w.wy * z + w.w * zy);
    result.dxx.push_back(w.wxx * z + 2.0 * w.wx * zx + w.w * a.p.dot(gxx));
    result.dxy.push_back(w.wxy * z + w.wx * zy + w.wy * zx + w.w * a.p.dot(gxy));
    result.dyy.push_back(w.wyy * z + 2.0 * w.wy * zy + w.w * a.p.dot(gyy));
  }
  return result;
}

}  // namespace al_maleh
