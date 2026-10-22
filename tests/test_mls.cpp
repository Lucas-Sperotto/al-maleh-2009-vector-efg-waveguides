#include "al_maleh/mls.hpp"

#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using al_maleh::MLSEvaluator;
using al_maleh::MLSResult;
using al_maleh::Node;
using al_maleh::Point2;

namespace {

void check(double actual, double expected, double tol, const std::string& name) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tol) {
    std::ostringstream error;
    error << std::setprecision(17) << name << " differs: got=" << actual
          << " expected=" << expected << " tolerance=" << tol;
    throw std::runtime_error(error.str());
  }
}

void reject(const std::function<void()>& operation, const std::string& name) {
  try { operation(); } catch (const std::exception&) { return; }
  throw std::runtime_error(name + " should have raised an exception");
}

std::vector<Node> regular_cloud() {
  std::vector<Node> nodes;
  for (int j = -2; j <= 2; ++j)
    for (int i = -2; i <= 2; ++i)
      nodes.push_back({{0.5 * i, 0.5 * j}, 2.5});
  return nodes;
}

void check_reproduction(int degree, Point2 x,
                        al_maleh::MLSWeight weight = al_maleh::MLSWeight::WendlandC2) {
  const auto cloud = regular_cloud();
  const MLSEvaluator evaluator(cloud, degree, 1e12, weight);
  const auto s = evaluator.evaluate(x);
  if (s.active_nodes.size() != cloud.size())
    throw std::runtime_error("Unexpected support for polynomial patch");
  double pu = 0, sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
  double value = 0, dx = 0, dy = 0, dxx = 0, dxy = 0, dyy = 0;
  for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
    const Point2 p = cloud[s.active_nodes[a]].position;
    // General quadratic polynomial. Linear patch uses its linear part.
    const double u = degree == 1
      ? 2.4 - 1.2 * p.x + 0.8 * p.y
      : 2.4 - 1.2 * p.x + 0.8 * p.y +
          0.7 * p.x * p.x - 0.3 * p.x * p.y + 1.1 * p.y * p.y;
    pu += s.phi[a]; sx += s.dx[a]; sy += s.dy[a];
    sxx += s.dxx[a]; sxy += s.dxy[a]; syy += s.dyy[a];
    value += s.phi[a] * u;
    dx += s.dx[a] * u; dy += s.dy[a] * u;
    dxx += s.dxx[a] * u; dxy += s.dxy[a] * u; dyy += s.dyy[a] * u;
  }
  const double target = degree == 1
      ? 2.4 - 1.2 * x.x + 0.8 * x.y
      : 2.4 - 1.2 * x.x + 0.8 * x.y +
          0.7 * x.x * x.x - 0.3 * x.x * x.y + 1.1 * x.y * x.y;
  check(pu, 1.0, 2e-11, "partition of unity");
  check(sx, 0, 2e-11, "sum dx");
  check(sy, 0, 2e-11, "sum dy");
  check(sxx, 0, 3e-10, "sum dxx");
  check(sxy, 0, 3e-10, "sum dxy");
  check(syy, 0, 3e-10, "sum dyy");
  check(value, target, 3e-11, "polynomial value");
  check(dx, -1.2 + (degree == 1 ? 0 : 1.4 * x.x - 0.3 * x.y),
        3e-10, "polynomial dx");
  check(dy, 0.8 + (degree == 1 ? 0 : -0.3 * x.x + 2.2 * x.y),
        3e-10, "polynomial dy");
  check(dxx, degree == 1 ? 0 : 1.4, 5e-9, "polynomial dxx");
  check(dxy, degree == 1 ? 0 : -0.3, 5e-9, "polynomial dxy");
  check(dyy, degree == 1 ? 0 : 2.2, 5e-9, "polynomial dyy");
}

void check_finite_differences(int degree,
                              al_maleh::MLSWeight weight = al_maleh::MLSWeight::WendlandC2) {
  const MLSEvaluator e(regular_cloud(), degree, 1e12, weight);
  const Point2 p{0.17, -0.13};
  const double h = 1.0e-5;
  const auto s = e.evaluate(p), xp = e.evaluate({p.x+h,p.y}),
             xm = e.evaluate({p.x-h,p.y}),
             yp = e.evaluate({p.x,p.y+h}), ym = e.evaluate({p.x,p.y-h});
  if (s.active_nodes != xp.active_nodes || s.active_nodes != xm.active_nodes ||
      s.active_nodes != yp.active_nodes || s.active_nodes != ym.active_nodes)
    throw std::runtime_error("Finite-difference support is not stable");
  for (std::size_t i = 0; i < s.active_nodes.size(); ++i) {
    check(s.dx[i], (xp.phi[i] - xm.phi[i])/(2*h), 2e-7, "FD dx");
    check(s.dy[i], (yp.phi[i] - ym.phi[i])/(2*h), 2e-7, "FD dy");
    check(s.dxx[i], (xp.dx[i] - xm.dx[i])/(2*h), 5e-6, "FD dxx");
    check(s.dxy[i], (yp.dx[i] - ym.dx[i])/(2*h), 5e-6, "FD dxy");
    check(s.dyy[i], (yp.dy[i] - ym.dy[i])/(2*h), 5e-6, "FD dyy");
  }
}

void check_invalid_clouds() {
  reject([]{ MLSEvaluator e({{{0,0},1}, {{1,0},1}}); },
         "too few nodes");
  reject([]{ MLSEvaluator e({{{0,0},1}, {{1,0},-1}, {{0,1},1}}); },
         "negative support");
  const MLSEvaluator collinear({{{0,0},2},{{1,0},2},{{2,0},2}});
  reject([&]{ (void)collinear.evaluate({1,0}); }, "collinear moment");
  const MLSEvaluator regular(regular_cloud(), 2);
  reject([&]{ (void)regular.evaluate({100,100}); }, "empty support");
  reject([]{ MLSEvaluator e(regular_cloud(), 3); }, "invalid degree");
  // Exactly at a node: the radial second derivative has a well-defined limit.
  const auto center = regular.evaluate({0,0});
  for (double x : center.dxx)
    if (!std::isfinite(x)) throw std::runtime_error("Hessian undefined at node");
}

// Independent linear-MLS solution on a square: reproduction constraints give
// phi=b+t*z, b=(1-x-y,x,y,0), z=(1,-1,-1,1). Minimizing sum(phi_i^2/w_i)
// determines t without assembling or solving the production moment matrix.
// This exercises the implemented weight, including its small positive tail.
void check_quartic_weight_through_shapes() {
  const Point2 x{0.2, 0.3};
  for (double target_q : {0.25, 0.5, 0.9, 0.999, 1.0 - 1e-5,
                          1.0 - 1e-6, 1.0 - 1e-8, 1.0, 1.1}) {
    const double radius = std::hypot(1.0 - x.x, 1.0 - x.y) / target_q;
    const std::vector<Node> nodes{{{0,0},3}, {{1,0},3}, {{0,1},3}, {{1,1},radius}};
    double weights[4];
    for (int i = 0; i < 4; ++i) {
      const double q = std::hypot(x.x - nodes[i].position.x,
                                  x.y - nodes[i].position.y) / nodes[i].support_radius;
      const double t = 1 - q;
      weights[i] = q < 1 ? t * t * t * (1 + 3 * q) : 0;
    }
    const double t = weights[3] == 0 ? 0 :
      -(0.5 / weights[0] - x.x / weights[1] - x.y / weights[2]) * weights[3] /
       (1 + weights[3] / weights[0] + weights[3] / weights[1] + weights[3] / weights[2]);
    const auto s = MLSEvaluator(nodes, 1, 1e12, al_maleh::MLSWeight::QuarticSpline).evaluate(x);
    double phi[4] = {};
    for (std::size_t a = 0; a < s.active_nodes.size(); ++a) phi[s.active_nodes[a]] = s.phi[a];
    check(phi[0], 0.5 + t, 2e-13, "quartic square phi0");
    check(phi[1], x.x - t, 2e-13, "quartic square phi1");
    check(phi[2], x.y - t, 2e-13, "quartic square phi2");
    check(phi[3], t, std::abs(t) * 2e-11, "quartic square small-tail phi3");
    if ((weights[3] > 0) != (std::find(s.active_nodes.begin(), s.active_nodes.end(), 3) != s.active_nodes.end()))
      throw std::runtime_error("quartic support lost a positive weight");
  }
}

// Al Maleh (2009) Eq. (4-2), checked through the actual MLS evaluator.
void check_quartic_spline_profile() {
  const std::vector<Node> one{{{0, 0}, 2.0}, {{5, 0}, 1.0}, {{0, 5}, 1.0}};
  const MLSEvaluator e(one, 1, 1e12, al_maleh::MLSWeight::QuarticSpline);
  if (std::string(al_maleh::mls_weight_name(e.weight())) != "quartic_spline" ||
      al_maleh::parse_mls_weight("wendland_c2") != al_maleh::MLSWeight::WendlandC2)
    throw std::runtime_error("weight names do not round-trip");
  reject([] { (void)al_maleh::parse_mls_weight("gaussian"); }, "unknown weight name");
  check_quartic_weight_through_shapes();
  for (int degree : {1, 2}) {
    for (Point2 p : {Point2{0.20, -0.15}, Point2{-0.35, 0.30}, Point2{0, 0}})
      check_reproduction(degree, p, al_maleh::MLSWeight::QuarticSpline);
    check_finite_differences(degree, al_maleh::MLSWeight::QuarticSpline);
  }
}

// Reference uses a fixed physical-coordinate Vandermonde matrix and weighted
// least-squares QR in extended precision, without forming A or using any of
// the production derivative formulas.
std::vector<long double> qr_reference(const std::vector<Node>& nodes, int degree,
                                     al_maleh::MLSWeight weight, Point2 x) {
  using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
  const int m = degree == 1 ? 3 : 6;
  auto basis = [m](Point2 p) {
    Vector b(m);
    const long double u = p.x, v = p.y;
    b.head(3) << 1, u, v;
    if (m == 6) b.tail(3) << u * u, u * v, v * v;
    return b;
  };
  Matrix design(nodes.size(), m);
  Vector root_weight(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const long double q = std::hypot(static_cast<long double>(x.x) - nodes[i].position.x,
                                     static_cast<long double>(x.y) - nodes[i].position.y) /
                          nodes[i].support_radius;
    const long double t = 1 - q;
    const long double w = q >= 1 ? 0 : weight == al_maleh::MLSWeight::QuarticSpline
      ? t * t * t * (1 + 3 * q) : t * t * t * t * (1 + 4 * q);
    root_weight(i) = std::sqrt(w);
    design.row(i) = root_weight(i) * basis(nodes[i].position).transpose();
  }
  const Matrix rhs = root_weight.asDiagonal();
  const Matrix coefficients = design.colPivHouseholderQr().solve(rhs);
  const Vector shapes = coefficients.transpose() * basis(x);
  return {shapes.data(), shapes.data() + shapes.size()};
}

std::vector<Node> asymmetric_cloud() {
  auto nodes = regular_cloud();
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    nodes[i].position.x += 0.07 * std::sin(1.3 * i);
    nodes[i].position.y += 0.05 * std::cos(0.9 * i);
    nodes[i].support_radius = 2.1 + 0.01 * i;
  }
  return nodes;
}

void check_fixed_basis_qr_reference() {
  const auto nodes = asymmetric_cloud();
  const Point2 x{0.17, -0.13};
  const double h = 1e-3;
  for (auto weight : {al_maleh::MLSWeight::WendlandC2, al_maleh::MLSWeight::QuarticSpline})
    for (int degree : {1, 2}) {
      const auto s = MLSEvaluator(nodes, degree, 1e12, weight).evaluate(x);
      const auto v = qr_reference(nodes, degree, weight, x);
      const auto xp = qr_reference(nodes, degree, weight, {x.x + h, x.y});
      const auto xm = qr_reference(nodes, degree, weight, {x.x - h, x.y});
      const auto xpp = qr_reference(nodes, degree, weight, {x.x + 2*h, x.y});
      const auto xmm = qr_reference(nodes, degree, weight, {x.x - 2*h, x.y});
      const auto yp = qr_reference(nodes, degree, weight, {x.x, x.y + h});
      const auto ym = qr_reference(nodes, degree, weight, {x.x, x.y - h});
      const auto ypp = qr_reference(nodes, degree, weight, {x.x, x.y + 2*h});
      const auto ymm = qr_reference(nodes, degree, weight, {x.x, x.y - 2*h});
      const int offsets[] = {-2, -1, 1, 2};
      const int factors[] = {1, -8, 8, -1};
      std::vector<long double> mixed(nodes.size(), 0);
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          const auto ij = qr_reference(nodes, degree, weight,
                                       {x.x + offsets[i]*h, x.y + offsets[j]*h});
          for (std::size_t a = 0; a < nodes.size(); ++a)
            mixed[a] += factors[i] * factors[j] * ij[a] / (144*h*h);
        }
      for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
        const auto i = s.active_nodes[a];
        check(s.phi[a], v[i], 2e-12, "fixed-basis QR value");
        check(s.dx[a], (xmm[i] - 8*xm[i] + 8*xp[i] - xpp[i]) / (12*h), 2e-8,
              "fixed-basis QR dx");
        check(s.dy[a], (ymm[i] - 8*ym[i] + 8*yp[i] - ypp[i]) / (12*h), 2e-8,
              "fixed-basis QR dy");
        check(s.dxx[a], (-xpp[i] + 16*xp[i] - 30*v[i] + 16*xm[i] - xmm[i]) / (12*h*h),
              5e-7, "fixed-basis QR dxx");
        check(s.dxy[a], mixed[i], 5e-7, "fixed-basis QR dxy");
        check(s.dyy[a], (-ypp[i] + 16*yp[i] - 30*v[i] + 16*ym[i] - ymm[i]) / (12*h*h),
              5e-7, "fixed-basis QR dyy");
      }
    }
}

void check_rigid_motion_and_unit_invariance() {
  const auto nodes = asymmetric_cloud();
  const Point2 x{0.17, -0.13};
  const double c = 0.6, d = 0.8;
  for (auto weight : {al_maleh::MLSWeight::WendlandC2, al_maleh::MLSWeight::QuarticSpline})
    for (int degree : {1, 2}) {
      const auto s = MLSEvaluator(nodes, degree, 1e12, weight).evaluate(x);
      for (double scale : {0.125, 8.0}) {
        auto transform = [scale, c, d](Point2 p) {
          return Point2{13.0 + scale * (c*p.x - d*p.y),
                        -9.0 + scale * (d*p.x + c*p.y)};
        };
        auto transformed = nodes;
        for (auto& n : transformed) {
          n.position = transform(n.position);
          n.support_radius *= scale;
        }
        const auto t = MLSEvaluator(transformed, degree, 1e12, weight).evaluate(transform(x));
        if (s.active_nodes != t.active_nodes) throw std::runtime_error("geometric transform changed support");
        for (std::size_t i = 0; i < s.active_nodes.size(); ++i) {
          check(t.phi[i], s.phi[i], 2e-12, "translation/rotation/unit invariant value");
          check(scale * t.dx[i], c*s.dx[i] - d*s.dy[i], 2e-11, "transformed dx");
          check(scale * t.dy[i], d*s.dx[i] + c*s.dy[i], 2e-11, "transformed dy");
          check(scale*scale * t.dxx[i], c*c*s.dxx[i] - 2*c*d*s.dxy[i] + d*d*s.dyy[i],
                5e-10, "transformed dxx");
          check(scale*scale * t.dxy[i], c*d*(s.dxx[i] - s.dyy[i]) + (c*c - d*d)*s.dxy[i],
                5e-10, "transformed dxy");
          check(scale*scale * t.dyy[i], d*d*s.dxx[i] + 2*c*d*s.dxy[i] + c*c*s.dyy[i],
                5e-10, "transformed dyy");
        }
      }
    }
}

void check_inactive_nodes_do_not_change_conditioning() {
  auto nodes = regular_cloud();
  const Point2 x{0.17, -0.13};
  for (auto weight : {al_maleh::MLSWeight::WendlandC2, al_maleh::MLSWeight::QuarticSpline})
    for (int degree : {1, 2}) {
      const auto local = MLSEvaluator(nodes, degree, 1e12, weight).evaluate(x);
      auto extended = nodes;
      extended.push_back({{1e6, 1e6}, 1000});
      const auto remote = MLSEvaluator(extended, degree, 1e12, weight).evaluate(x);
      if (local.active_nodes != remote.active_nodes)
        throw std::runtime_error("inactive remote node changed support");
      check(remote.condition_number, local.condition_number, 0, "local moment conditioning");
      for (std::size_t i = 0; i < local.active_nodes.size(); ++i) {
        check(remote.phi[i], local.phi[i], 0, "local value");
        check(remote.dx[i], local.dx[i], 0, "local dx");
        check(remote.dy[i], local.dy[i], 0, "local dy");
        check(remote.dxx[i], local.dxx[i], 0, "local dxx");
        check(remote.dxy[i], local.dxy[i], 0, "local dxy");
        check(remote.dyy[i], local.dyy[i], 0, "local dyy");
      }
    }
}

// Central differences: the error must first fall like h^2 (truncation) and
// then rise again (round-off). Both regimes are asserted, not just one h.
void check_fd_step_study() {
  const MLSEvaluator e(regular_cloud(), 2);
  const Point2 p{0.17, -0.13};
  const auto s = e.evaluate(p);
  std::vector<double> err;
  for (double h : {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8}) {
    const auto xp = e.evaluate({p.x + h, p.y}), xm = e.evaluate({p.x - h, p.y});
    if (xp.active_nodes != s.active_nodes || xm.active_nodes != s.active_nodes)
      throw std::runtime_error("FD step changed the MLS support");
    double worst = 0;
    for (std::size_t i = 0; i < s.active_nodes.size(); ++i) {
      worst = std::max(worst, std::abs((xp.phi[i] - xm.phi[i]) / (2 * h) - s.dx[i]));
      worst = std::max(worst, std::abs((xp.dx[i] - xm.dx[i]) / (2 * h) - s.dxx[i]));
    }
    err.push_back(worst);
  }
  // Truncation regime: one decade in h gives ~two decades in error.
  const double rate = std::log10(err[0] / err[1]);
  if (!(rate > 1.7 && rate < 2.3))
    throw std::runtime_error("FD truncation rate is not second order: " + std::to_string(rate));
  const double best = *std::min_element(err.begin(), err.end());
  if (!(best < 1e-8)) throw std::runtime_error("FD never reaches 1e-8 agreement");
  if (!(err.back() > 10 * best))
    throw std::runtime_error("FD round-off regime not observed at h=1e-8");
}

// Corner, edge and interior points of the rectangular 5x3 cloud used by the
// thesis (support factor 1.5 and 3.3): reproduction must hold up to the wall.
void check_rectangular_boundary_reproduction() {
  for (auto weight : {al_maleh::MLSWeight::WendlandC2, al_maleh::MLSWeight::QuarticSpline})
    for (int degree : {1, 2})
      for (double factor : {1.5, 3.3}) {
        if (degree == 2 && factor < 2.0) continue;  // < 6 neighbours at corners
        std::vector<Node> nodes;
        for (int j = 0; j < 3; ++j)
          for (int i = 0; i < 5; ++i) nodes.push_back({{0.5 * i, 0.5 * j}, 0.5 * factor});
        const MLSEvaluator e(nodes, degree, 1e12, weight);
        for (Point2 x : {Point2{0, 0}, Point2{2, 1}, Point2{1, 0}, Point2{0, 0.5},
                         Point2{1.37, 0.81}, Point2{0.05, 0.97}}) {
          const auto s = e.evaluate(x);
          double v = 0, gx = 0, gy = 0, lap = 0, pu = 0;
          for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
            const Point2 n = nodes[s.active_nodes[a]].position;
            const double u = degree == 1 ? 0.3 - 1.1 * n.x + 0.7 * n.y
                                         : 0.3 - 1.1 * n.x + 0.7 * n.y + n.x * n.x + n.y * n.y;
            pu += s.phi[a]; v += s.phi[a] * u; gx += s.dx[a] * u; gy += s.dy[a] * u;
            lap += (s.dxx[a] + s.dyy[a]) * u;
          }
          const double exact = degree == 1 ? 0.3 - 1.1 * x.x + 0.7 * x.y
                                           : 0.3 - 1.1 * x.x + 0.7 * x.y + x.x * x.x + x.y * x.y;
          check(pu, 1, 1e-12, "boundary partition of unity");
          check(v, exact, 1e-11, "boundary reproduction");
          check(gx, -1.1 + (degree == 2 ? 2 * x.x : 0), 1e-9, "boundary dx");
          check(gy, 0.7 + (degree == 2 ? 2 * x.y : 0), 1e-9, "boundary dy");
          check(lap, degree == 2 ? 4.0 : 0.0, 1e-8, "boundary Laplacian of reproduced field");
        }
      }
}

// Regression for the basis-centering defect: with a basis centered on the
// cloud and scaled by its size, cond(A) grew like (L/h)^(2 degree) (3.9e3,
// 8.0e4, 1.5e6 for 5x3, 9x5, 17x9 degree-2 clouds). With the query-point
// shifted basis it must be refinement-independent.
void check_conditioning_independent_of_refinement() {
  double smallest = 1e300, largest = 0;
  for (int k : {1, 2, 4, 8}) {
    const int nx = 4 * k + 1, ny = 2 * k + 1;
    std::vector<Node> nodes;
    const double h = 2.0 / (nx - 1);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) nodes.push_back({{h * i, h * j}, 3.3 * h});
    const MLSEvaluator e(nodes, 2);
    double worst = 0;
    for (Point2 x : {Point2{0, 0}, Point2{1, 0.5}, Point2{0.3 * h, 1 - 0.2 * h},
                     Point2{2 - 0.5 * h, 0.5 * h}})
      worst = std::max(worst, e.evaluate(x).condition_number);
    smallest = std::min(smallest, worst);
    largest = std::max(largest, worst);
  }
  if (!(largest < 1e4) || !(largest / smallest < 1.5))
    throw std::runtime_error("MLS moment conditioning depends on refinement: " +
                             std::to_string(smallest) + " .. " + std::to_string(largest));
}

// Second-order consistency: degree 2 reproduces lap(x^2+y^2)=4; degree 1 does
// NOT (defect independent of h). The latter is a characterization test: it
// keeps visible the mechanism behind the non-convergent M5 (Ch.6) results.
void check_laplacian_consistency() {
  double defect_coarse = 0, defect_fine = 0;
  for (int k : {1, 4}) {
    const int nx = 4 * k + 1, ny = 2 * k + 1;
    const double h = 2.0 / (nx - 1);
    std::vector<Node> nodes;
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) nodes.push_back({{h * i, h * j}, 1.5 * h});
    const MLSEvaluator linear(nodes, 1, 1e12, al_maleh::MLSWeight::QuarticSpline);
    // Same relative position inside a central nodal cell at both resolutions.
    const Point2 x{1.0 + 0.3 * h, 0.5 + 0.2 * h};
    const auto s = linear.evaluate(x);
    double lap = 0;
    for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
      const Point2 n = nodes[s.active_nodes[a]].position;
      lap += (n.x * n.x + n.y * n.y) * (s.dxx[a] + s.dyy[a]);
    }
    (k == 1 ? defect_coarse : defect_fine) = std::abs(lap - 4.0);
  }
  if (!(defect_coarse > 0.1) || std::abs(defect_coarse - defect_fine) > 1e-6)
    throw std::runtime_error("degree-1 Laplacian defect should be O(1) and h-independent");
}

} // namespace

int main() {
  try {
    for (int degree : {1, 2}) {
      for (Point2 p : {Point2{0.20,-0.15}, Point2{-0.35,0.30}, Point2{0,0}})
        check_reproduction(degree, p);
      check_finite_differences(degree);
    }
    check_invalid_clouds();
    check_quartic_spline_profile();
    check_fixed_basis_qr_reference();
    check_rigid_motion_and_unit_invariance();
    check_inactive_nodes_do_not_change_conditioning();
    check_fd_step_study();
    check_rectangular_boundary_reproduction();
    check_conditioning_independent_of_refinement();
    check_laplacian_consistency();
    std::cout << "PASS: MLS partition, linear/quadratic patches, gradient, Hessian, FD and failure cases; "
                 "quartic spline Eq.(4-2), FD step study, rectangular corners, refinement-independent "
                 "conditioning, second-order consistency, independent fixed-basis QR, "
                 "rigid-motion/unit invariance, local support normalization\n";
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
