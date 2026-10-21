#include "al_maleh/mls.hpp"

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

void check_reproduction(int degree, Point2 x) {
  const auto cloud = regular_cloud();
  const MLSEvaluator evaluator(cloud, degree);
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

void check_finite_differences(int degree) {
  const MLSEvaluator e(regular_cloud(), degree);
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

} // namespace

int main() {
  try {
    for (int degree : {1, 2}) {
      for (Point2 p : {Point2{0.20,-0.15}, Point2{-0.35,0.30}, Point2{0,0}})
        check_reproduction(degree, p);
      check_finite_differences(degree);
    }
    check_invalid_clouds();
    std::cout << "PASS: MLS partition, linear/quadratic patches, gradient, Hessian, FD and failure cases\n";
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
