// SPDX-License-Identifier: MIT
#include "al_maleh/mode_tracking.hpp"

#include <Eigen/QR>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace al_maleh {
namespace {

void check_samples(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
                   const Eigen::VectorXd& w) {
  if (A.rows() < 1 || A.rows() != B.rows() || A.rows() != w.size() ||
      A.cols() < 1 || B.cols() < 1)
    throw std::invalid_argument("Sampled fields must share the same rows");
  if (!A.allFinite() || !B.allFinite() || !w.allFinite() || w.minCoeff() <= 0.0)
    throw std::invalid_argument("Nonfinite samples or nonpositive weights");
}

// Orthonormal basis of span(sqrt(W) X), rank-revealing.
Eigen::MatrixXd weighted_orthonormal(const Eigen::MatrixXd& X,
                                     const Eigen::VectorXd& w,
                                     double tolerance) {
  const Eigen::MatrixXd scaled = w.array().sqrt().matrix().asDiagonal() * X;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(scaled);
  const Eigen::MatrixXd R = qr.matrixR().triangularView<Eigen::Upper>();
  const Eigen::Index diagonal = std::min(R.rows(), R.cols());
  const double largest = diagonal > 0 ? std::abs(R(0, 0)) : 0.0;
  if (!(largest > 0.0))
    throw std::invalid_argument("Sampled subspace has zero weighted norm");
  Eigen::Index rank = 0;
  while (rank < diagonal && std::abs(R(rank, rank)) > tolerance * largest)
    ++rank;
  const Eigen::MatrixXd Q = qr.householderQ() *
      Eigen::MatrixXd::Identity(scaled.rows(), rank);
  return Q;
}

} // namespace

Eigen::MatrixXd normalized_overlaps(const Eigen::MatrixXd& A,
                                    const Eigen::MatrixXd& B,
                                    const Eigen::VectorXd& w) {
  check_samples(A, B, w);
  const Eigen::MatrixXd WB = w.asDiagonal() * B;
  const Eigen::MatrixXd cross = A.transpose() * WB;
  const Eigen::VectorXd na =
      (A.array().square().colwise() * w.array()).colwise().sum().sqrt();
  const Eigen::VectorXd nb =
      (B.array().square().colwise() * w.array()).colwise().sum().sqrt();
  if (!(na.minCoeff() > 0.0) || !(nb.minCoeff() > 0.0))
    throw std::invalid_argument("A sampled field has zero weighted norm");
  Eigen::MatrixXd out(A.cols(), B.cols());
  for (Eigen::Index i = 0; i < A.cols(); ++i)
    for (Eigen::Index j = 0; j < B.cols(); ++j) {
      const double c = std::abs(cross(i, j)) / (na[i] * nb[j]);
      if (!std::isfinite(c) || c > 1.0 + 1e-8)
        throw std::runtime_error("Overlap violates Cauchy-Schwarz bound");
      out(i, j) = std::min(1.0, c);
    }
  return out;
}

Eigen::VectorXd principal_cosines(const Eigen::MatrixXd& A,
                                  const Eigen::MatrixXd& B,
                                  const Eigen::VectorXd& w,
                                  double rank_tolerance) {
  check_samples(A, B, w);
  if (!(rank_tolerance > 0.0) || !(rank_tolerance < 1.0))
    throw std::invalid_argument("Invalid principal-angle rank tolerance");
  const Eigen::MatrixXd QA = weighted_orthonormal(A, w, rank_tolerance);
  const Eigen::MatrixXd QB = weighted_orthonormal(B, w, rank_tolerance);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(QA.transpose() * QB);
  Eigen::VectorXd s = svd.singularValues();  // already descending
  if (!s.allFinite() || (s.size() > 0 && s.maxCoeff() > 1.0 + 1e-8))
    throw std::runtime_error("Invalid principal cosines");
  return s.cwiseMin(1.0).cwiseMax(0.0);
}

std::vector<int> max_overlap_assignment(const Eigen::MatrixXd& overlaps) {
  if (overlaps.rows() < 1 || overlaps.cols() < 1 || !overlaps.allFinite())
    throw std::invalid_argument("Invalid overlap matrix for assignment");
  const bool transpose = overlaps.rows() > overlaps.cols();
  const Eigen::MatrixXd S = transpose ? Eigen::MatrixXd(overlaps.transpose())
                                      : overlaps;
  const int n = static_cast<int>(S.rows()), m = static_cast<int>(S.cols());
  // Hungarian algorithm (minimization of 1-overlap), n <= m, 1-based.
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> u(n + 1, 0.0), v(m + 1, 0.0);
  std::vector<int> p(m + 1, 0), way(m + 1, 0);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> minv(m + 1, inf);
    std::vector<char> used(m + 1, 0);
    do {
      used[j0] = 1;
      const int i0 = p[j0];
      double delta = inf;
      int j1 = 0;
      for (int j = 1; j <= m; ++j) {
        if (used[j]) continue;
        const double cur = (1.0 - S(i0 - 1, j - 1)) - u[i0] - v[j];
        if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
        if (minv[j] < delta) { delta = minv[j]; j1 = j; }
      }
      for (int j = 0; j <= m; ++j) {
        if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
        else minv[j] -= delta;
      }
      j0 = j1;
    } while (p[j0] != 0);
    do {
      const int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0 != 0);
  }
  std::vector<int> row_to_col(static_cast<std::size_t>(n), -1);
  for (int j = 1; j <= m; ++j)
    if (p[j] != 0) row_to_col[static_cast<std::size_t>(p[j] - 1)] = j - 1;
  if (!transpose) return row_to_col;
  std::vector<int> result(static_cast<std::size_t>(overlaps.rows()), -1);
  for (int i = 0; i < n; ++i)
    if (row_to_col[static_cast<std::size_t>(i)] >= 0)
      result[static_cast<std::size_t>(row_to_col[static_cast<std::size_t>(i)])] = i;
  return result;
}

std::vector<Eigen::Index> eigenvalue_clusters(const Eigen::VectorXd& values,
                                              double relative_tolerance) {
  if (!values.allFinite() || !(relative_tolerance >= 0.0))
    throw std::invalid_argument("Invalid eigenvalue clustering input");
  std::vector<Eigen::Index> starts;
  for (Eigen::Index k = 0; k < values.size(); ++k) {
    if (k > 0 && values[k] < values[k - 1])
      throw std::invalid_argument("Eigenvalues must be sorted ascending");
    const double scale = k > 0 ? std::max({std::abs(values[k]),
                                           std::abs(values[k - 1]), 1e-300})
                               : 1.0;
    if (k == 0 || values[k] - values[k - 1] > relative_tolerance * scale)
      starts.push_back(k);
  }
  starts.push_back(values.size());
  return starts;
}

} // namespace al_maleh
