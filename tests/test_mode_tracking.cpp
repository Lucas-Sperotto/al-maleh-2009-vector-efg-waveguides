// SPDX-License-Identifier: MIT
// Synthetic contracts for M7 mode tracking: sign/scale/basis-rotation
// invariance, principal angles and optimal assignment.
#include "al_maleh/mode_tracking.hpp"

#include <Eigen/QR>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {

void near(double actual, double expected, double tol, const std::string& what) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
    throw std::runtime_error(what + " expected=" + std::to_string(expected) +
                             " got=" + std::to_string(actual));
}

void reject(const std::function<void()>& fn, const std::string& what) {
  try { fn(); } catch (const std::exception&) { return; }
  throw std::runtime_error(what + " should have been rejected");
}

// Deterministic smooth "fields" on 200 sample rows with uneven weights.
Eigen::MatrixXd fields(int columns) {
  const int rows = 200;
  Eigen::MatrixXd F(rows, columns);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < columns; ++j)
      F(i, j) = std::sin(0.031 * (i + 1) * (j + 1) + 0.2 * j) +
                0.1 * std::cos(0.007 * (i + 3) * (j + 2));
  return F;
}

Eigen::VectorXd weights(int rows) {
  Eigen::VectorXd w(rows);
  for (int i = 0; i < rows; ++i) w[i] = 0.5 + 0.4 * std::sin(0.05 * i) * std::sin(0.05 * i);
  return w;
}

void test_sign_and_scale_invariance() {
  const Eigen::MatrixXd A = fields(4), B = fields(5).rightCols(3);
  const Eigen::VectorXd w = weights(200);
  const Eigen::MatrixXd base = normalized_overlaps(A, B, w);
  Eigen::MatrixXd As = A, Bs = B;
  As.col(1) *= -3.7; As.col(3) *= 1e-6; Bs.col(0) *= -1.0; Bs.col(2) *= 42.0;
  near((normalized_overlaps(As, Bs, w) - base).cwiseAbs().maxCoeff(), 0, 1e-12,
       "overlap sign/scale invariance");
  near(normalized_overlaps(A, A, w).diagonal().minCoeff(), 1, 1e-12,
       "self overlap equals one");
  if (base.maxCoeff() > 1.0 || base.minCoeff() < 0.0)
    throw std::runtime_error("overlap outside [0,1]");
}

void test_weights_are_used() {
  // Two fields orthogonal in the weighted product, not in the plain sum.
  const int rows = 4;
  Eigen::VectorXd w(rows); w << 1, 3, 1, 3;
  Eigen::MatrixXd a(rows, 1), b(rows, 1);
  a << 3, 1, 0, 0;
  b << 1, -1, 0, 0;  // weighted: 3*1*1 + 1*(-1)*3 = 0
  near(normalized_overlaps(a, b, w)(0, 0), 0, 1e-14, "weighted orthogonality");
  near(normalized_overlaps(a, b, Eigen::VectorXd::Ones(rows))(0, 0),
       2.0 / std::sqrt(10.0 * 2.0), 1e-14, "unweighted overlap differs");
}

void test_principal_angles() {
  const Eigen::VectorXd w = weights(200);
  const Eigen::MatrixXd A = fields(2);
  // Rotation + scaling inside the subspace must not change the cosines.
  Eigen::Matrix2d R;
  const double t = 0.913;
  R << std::cos(t), -std::sin(t), std::sin(t), std::cos(t);
  Eigen::Matrix2d S; S << 2.0, 0.3, -0.1, -5.0;  // nonsingular, nonorthogonal
  const Eigen::VectorXd same = principal_cosines(A, (A * R * S).eval(), w);
  near(same.minCoeff(), 1, 1e-10, "identical subspace cosines");
  // Individual overlaps DO change under the same rotation.
  const Eigen::MatrixXd single = normalized_overlaps(A, (A * R).eval(), w);
  if (single(0, 0) > 0.999)
    throw std::runtime_error("synthetic rotation unexpectedly preserved columns");
  // Known angle: build W-orthonormal e1, e2 and b = cos(theta) e1 + sin(theta) e2.
  Eigen::MatrixXd E = fields(3);
  const Eigen::VectorXd sw = w.array().sqrt();
  Eigen::MatrixXd scaled = sw.asDiagonal() * E;
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(scaled);
  const Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(200, 3);
  const Eigen::MatrixXd basis = sw.cwiseInverse().asDiagonal() * Q;  // W-orthonormal
  const double theta = 0.37;
  Eigen::MatrixXd a = basis.col(0);
  Eigen::MatrixXd b = std::cos(theta) * basis.col(0) + std::sin(theta) * basis.col(1);
  near(principal_cosines(a, b, w)[0], std::cos(theta), 1e-12, "known principal angle");
  // Two-dimensional subspaces sharing one direction, second at 90 degrees.
  Eigen::MatrixXd P(200, 2), T(200, 2);
  P << basis.col(0), basis.col(1);
  T << basis.col(0), basis.col(2);
  const Eigen::VectorXd c = principal_cosines(P, T, w);
  near(c[0], 1, 1e-12, "shared direction");
  near(c[1], 0, 1e-12, "orthogonal direction");
  // Dependent columns are dropped instead of producing spurious cosines.
  Eigen::MatrixXd D(200, 3);
  D << basis.col(0), basis.col(1), basis.col(0) + 2 * basis.col(1);
  if (principal_cosines(D, P, w).size() != 2)
    throw std::runtime_error("dependent column not removed");
  reject([&] { (void)principal_cosines(Eigen::MatrixXd::Zero(200, 2), P, w); },
         "zero subspace");
  Eigen::VectorXd bad = w; bad[3] = 0.0;
  reject([&] { (void)principal_cosines(P, P, bad); }, "nonpositive weight");
}

void test_assignment() {
  // Permuted near-identity: optimal assignment recovers the permutation.
  const std::vector<int> perm{3, 0, 4, 1, 2};
  Eigen::MatrixXd O = Eigen::MatrixXd::Constant(5, 5, 0.05);
  for (int i = 0; i < 5; ++i) O(i, perm[static_cast<std::size_t>(i)]) = 0.9 - 0.01 * i;
  const auto a = max_overlap_assignment(O);
  for (int i = 0; i < 5; ++i)
    if (a[static_cast<std::size_t>(i)] != perm[static_cast<std::size_t>(i)])
      throw std::runtime_error("assignment did not recover permutation");
  // Greedy (row 0 takes its best column) is suboptimal here; Hungarian is not.
  Eigen::MatrixXd G(2, 2);
  G << 0.90, 0.85,
       0.80, 0.10;
  const auto g = max_overlap_assignment(G);
  if (g[0] != 1 || g[1] != 0)
    throw std::runtime_error("assignment is not globally optimal");
  // Rectangular: more previous modes than current ones.
  Eigen::MatrixXd Rr(3, 2);
  Rr << 0.1, 0.95,
        0.2, 0.1,
        0.9, 0.3;
  const auto r = max_overlap_assignment(Rr);
  if (r[0] != 1 || r[1] != -1 || r[2] != 0)
    throw std::runtime_error("rectangular assignment wrong");
}

void test_clusters() {
  Eigen::VectorXd v(6);
  v << 1.0, 1.0 + 1e-9, 2.0, 5.0, 5.0 + 4e-7, 5.0 + 9e-7;
  const auto c = eigenvalue_clusters(v, 1e-6);
  const std::vector<Eigen::Index> expected{0, 2, 3, 6};
  if (c != expected) throw std::runtime_error("eigenvalue clusters wrong");
  Eigen::VectorXd unsorted(2); unsorted << 2.0, 1.0;
  reject([&] { (void)eigenvalue_clusters(unsorted, 1e-6); }, "unsorted eigenvalues");
}

} // namespace

int main() {
  try {
    test_sign_and_scale_invariance();
    test_weights_are_used();
    test_principal_angles();
    test_assignment();
    test_clusters();
    std::cout << "PASS: M7 overlap sign/scale invariance, weighted product, "
                 "principal-angle basis invariance, Hungarian assignment, clusters\n";
  } catch (const std::exception& e) {
    std::cerr << "M7 FAIL: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
