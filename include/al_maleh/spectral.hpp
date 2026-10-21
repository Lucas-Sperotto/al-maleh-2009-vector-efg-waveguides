// SPDX-License-Identifier: MIT
#pragma once
#include "al_maleh/cartesian.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
namespace al_maleh {
struct SpectralOptions {
  double zero_absolute=1e-9, zero_relative=1e-9;
  double max_mass_condition=1e13, symmetry_tolerance=1e-10;
  int max_dense_dofs=600;
};
struct Spectrum {
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors; // M-orthonormal columns
  std::vector<double> residuals;
  int numerical_kernel_dimension=0;
  double zero_cutoff=0;
  double mass_condition=0, min_mass_eigenvalue=0;
};
Spectrum solve_generalized_symmetric(
  const Eigen::SparseMatrix<double>& K,
  const Eigen::SparseMatrix<double>& M,
  const SpectralOptions& options={});
struct TMMode { int m=0,n=0; double lambda=0,kc=0; };
std::vector<TMMode> rectangular_tm_reference(
  Rectangle domain,int mmax,int nmax);
struct FieldMetrics {
  double field_l2=0,curl_l2=0,div_l2=0;
  double divergence_ratio=0,normal_ratio=0;
};
FieldMetrics modal_diagnostics(const MLSEvaluator& mls,
  const BackgroundQuadrature& rule,
  const Eigen::VectorXd& coefficients,double eigenvalue);
} // namespace al_maleh
