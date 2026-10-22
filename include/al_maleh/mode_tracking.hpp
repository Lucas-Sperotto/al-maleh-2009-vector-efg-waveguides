// SPDX-License-Identifier: MIT
// M7: mode identity across discretization parameters, measured on a common
// spatial sample rule (never by position in the sorted eigenvalue list).
#pragma once

#include <Eigen/Core>
#include <vector>

namespace al_maleh {

// Columns of A and B are sampled fields (any stacking of components) on the
// SAME points; w holds one positive weight per row (quadrature weights).
// Entry (i,j) = |<a_i,b_j>_w| / (||a_i||_w ||b_j||_w), in [0,1].
// Invariant to the sign and the scale of every column.
Eigen::MatrixXd normalized_overlaps(const Eigen::MatrixXd& A,
                                    const Eigen::MatrixXd& B,
                                    const Eigen::VectorXd& w);

// Cosines of the principal angles between span(A) and span(B) in the
// weighted L2 inner product, in descending order (min(rank A, rank B)
// values). Invariant under A->A R and B->B S for nonsingular R and S, so
// it is the correct comparison for (near-)degenerate eigenspaces.
// Columns are equilibrated before the numerical rank test; dependent columns
// (relative tolerance) are dropped. Ill-conditioned changes of basis can still
// change the numerical rank at the requested tolerance.
Eigen::VectorXd principal_cosines(const Eigen::MatrixXd& A,
                                  const Eigen::MatrixXd& B,
                                  const Eigen::VectorXd& w,
                                  double rank_tolerance = 1e-10);

// Maximum-total-overlap one-to-one assignment (Hungarian algorithm).
// result[i] = matched column of row i, or -1 when rows exceed columns.
std::vector<int> max_overlap_assignment(const Eigen::MatrixXd& overlaps);

// Assigned overlap minus the largest competing entry in its row OR column.
// A large row-only margin cannot establish identity when another source mode
// competes for the same destination. This margin is a confidence diagnostic,
// not the difference between the two best global assignment objectives.
double overlap_assignment_margin(const Eigen::MatrixXd& overlaps,
                                Eigen::Index row, Eigen::Index col);

// Groups ascending eigenvalues into clusters whose consecutive relative
// gaps are <= relative_tolerance; returns the first index of each cluster
// plus a final sentinel equal to values.size().
std::vector<Eigen::Index> eigenvalue_clusters(const Eigen::VectorXd& values,
                                              double relative_tolerance);

} // namespace al_maleh
