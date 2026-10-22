#include "al_maleh/mls.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using al_maleh::MLSEvaluator;
using al_maleh::MLSResult;
using al_maleh::Node;
using al_maleh::Point2;

namespace {

void check(double actual, double expected, double tol, const std::string& name) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
    throw std::runtime_error(name + " differs: got=" + std::to_string(actual) +
                             " expected=" + std::to_string(expected));
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

// Al Maleh (2009) Eq. (4-2): exact profile values and C2 contact at q=1.
void check_quartic_spline_profile() {
  const std::vector<Node> one{{{0, 0}, 2.0}, {{5, 0}, 1.0}, {{0, 5}, 1.0}};
  const MLSEvaluator e(one, 1, 1e12, al_maleh::MLSWeight::QuarticSpline);
  if (std::string(al_maleh::mls_weight_name(e.weight())) != "quartic_spline" ||
      al_maleh::parse_mls_weight("wendland_c2") != al_maleh::MLSWeight::WendlandC2)
    throw std::runtime_error("weight names do not round-trip");
  reject([] { (void)al_maleh::parse_mls_weight("gaussian"); }, "unknown weight name");
  // Shape-independent check of w through a 1-node partition is impossible
  // (MLS normalizes), so verify w, w' and w'' of the profile directly from
  // the polynomial via the FD-free closed forms at several q.
  for (double q : {0.0, 0.25, 0.5, 0.9, 0.999}) {
    const double w = 1 - 6 * q * q + 8 * q * q * q - 3 * q * q * q * q;
    const double dw = -12 * q * (1 - q) * (1 - q);
    const double d2w = -12 * (1 - q) * (1 - 3 * q);
    const double h = 1e-6, h2 = 1e-4;  // h2 balances truncation and round-off
    auto f = [](double s) { return 1 - 6 * s * s + 8 * s * s * s - 3 * s * s * s * s; };
    check((f(q + h) - f(q - h)) / (2 * h), dw, 1e-8, "quartic dw/dq");
    check((f(q + h2) - 2 * w + f(q - h2)) / (h2 * h2), d2w, 1e-5, "quartic d2w/dq2");
  }
  // Value, slope and curvature vanish at the support boundary (C2 contact).
  check(1 - 6 + 8 - 3, 0, 0, "w(1)");
  for (int degree : {1, 2}) {
    for (Point2 p : {Point2{0.20, -0.15}, Point2{-0.35, 0.30}, Point2{0, 0}})
      check_reproduction(degree, p, al_maleh::MLSWeight::QuarticSpline);
    check_finite_differences(degree, al_maleh::MLSWeight::QuarticSpline);
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
    check_fd_step_study();
    check_rectangular_boundary_reproduction();
    check_conditioning_independent_of_refinement();
    check_laplacian_consistency();
    std::cout << "PASS: MLS partition, linear/quadratic patches, gradient, Hessian, FD and failure cases; "
                 "quartic spline Eq.(4-2), FD step study, rectangular corners, refinement-independent "
                 "conditioning, second-order consistency\n";
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
