// SPDX-License-Identifier: MIT
#pragma once
#include "al_maleh/cartesian.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <filesystem>
#include <string>
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

// Persist every column without filtering or changing the solver's original
// M normalization. CSV spectrum indices are zero based; Matrix Market rows
// and columns follow the format's one-based convention. Hx DOFs precede Hy.
struct SpectralExportOptions {
  int nx=0,ny=0;
  double support_factor=0;
  std::string origin="cartesian_EFG";
  int grid_nx=80,grid_ny=40;
};
void export_spectrum(const std::filesystem::path& directory,
  const MLSEvaluator& mls,const BackgroundQuadrature& rule,
  const CartesianAssembly& matrices,const Spectrum& spectrum,
  const SpectralExportOptions& options);
} // namespace al_maleh
