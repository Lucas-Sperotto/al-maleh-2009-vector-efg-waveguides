// SPDX-License-Identifier: MIT
// Independent theory/code audit (2026-10-22): reproducible numerical studies.
// Every subcommand writes ONE CSV; nothing modifies K/M or historical data.
//
//   independent_audit ch4-scan    weight degree fmin fmax fstep out.csv
//   independent_audit laplacian   out.csv
//   independent_audit convergence out.csv
//   independent_audit rank-sweep  out.csv
//   independent_audit doublet     nx ny out.csv
#include "al_maleh/constraints.hpp"
#include "al_maleh/gradient_audit.hpp"
#include "al_maleh/rot_irrot.hpp"
#include "al_maleh/spectral.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {

const Rectangle kRect{0, 2, 0, 1};
const double kPi = std::acos(-1.0);

struct Profile { MLSWeight weight; double factor; int degree; };

MLSEvaluator cloud(int nx, int ny, const Profile& p) {
  const double h = std::max(2.0 / (nx - 1), 1.0 / (ny - 1));
  return MLSEvaluator(make_regular_cloud(kRect, nx, ny, p.factor * h), p.degree,
                      1e12, p.weight);
}

std::ofstream open_csv(const std::string& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot create " + path);
  out << std::setprecision(17);
  return out;
}

// curl(H) = sum_i c_i lap(N_i) for M5. On a fixed reference rule we store
// L = int lap(N_i) lap(N_j), l_t = int lap(N_i) psi_t and r_t = int psi_t^2,
// so the curl correlation of any coefficient vector is |c.l_t|/sqrt(c.L.c r_t)
// (same definition as rot_irrot_rect, on a finer, converged rule).
struct CurlForms {
  Eigen::MatrixXd L;
  Eigen::MatrixXd l;   // n x 2 (TM11, TM21)
  Eigen::Vector2d r = Eigen::Vector2d::Zero();
};

CurlForms curl_forms(const MLSEvaluator& mls, const BackgroundQuadrature& q) {
  const auto n = static_cast<Eigen::Index>(mls.nodes().size());
  CurlForms f;
  f.L = Eigen::MatrixXd::Zero(n, n);
  f.l = Eigen::MatrixXd::Zero(n, 2);
  for_each_quadrature_point(q, [&](Point2 p, double w) {
    const auto s = mls.evaluate(p);
    const double psi[2] = {std::sin(kPi * p.x / 2) * std::sin(kPi * p.y),
                           std::sin(kPi * p.x) * std::sin(kPi * p.y)};
    for (int t = 0; t < 2; ++t) f.r[t] += w * psi[t] * psi[t];
    for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
      const auto i = static_cast<Eigen::Index>(s.active_nodes[a]);
      const double li = s.dxx[a] + s.dyy[a];
      for (int t = 0; t < 2; ++t) f.l(i, t) += w * li * psi[t];
      for (std::size_t b = 0; b < s.active_nodes.size(); ++b)
        f.L(i, static_cast<Eigen::Index>(s.active_nodes[b])) +=
            w * li * (s.dxx[b] + s.dyy[b]);
    }
  });
  return f;
}

Eigen::VectorXd m5_curl_rho(const CurlForms& f, const RankReducedSpectrum& s, int t) {
  const Eigen::Index n = f.L.rows();
  const Eigen::MatrixXd C = s.eigenvectors.bottomRows(n);  // rotational coefficients
  const Eigen::VectorXd cross = C.transpose() * f.l.col(t);
  const Eigen::VectorXd energy = (C.array() * (f.L * C).array()).colwise().sum();
  Eigen::VectorXd rho(C.cols());
  for (Eigen::Index j = 0; j < rho.size(); ++j)
    rho[j] = energy[j] > 1e-30 ? std::min(1.0, std::abs(cross[j]) / std::sqrt(energy[j] * f.r[t]))
                               : 0.0;
  return rho;
}

// ---------------------------------------------------------------- ch4-scan
// First Cartesian and rotational/irrotational eigenvalues of the 3x5 grid
// (author naming: 3 rows x 5 columns) for a range of radius factors.
void ch4_scan(const std::vector<std::string>& a) {
  if (a.size() != 6) throw std::invalid_argument("ch4-scan weight degree fmin fmax fstep out.csv");
  const MLSWeight weight = parse_mls_weight(a[0]);
  const int degree = std::stoi(a[1]);
  const double f0 = std::stod(a[2]), f1 = std::stod(a[3]), df = std::stod(a[4]);
  if (!(df > 0) || !(f1 >= f0) || !(f0 > 0)) throw std::invalid_argument("Invalid factor range");
  auto out = open_csv(a[5]);
  out << "formulation,weight,degree,factor,gauss,position,lambda,status\n";
  const int steps = static_cast<int>(std::floor((f1 - f0) / df + 1e-9));
  for (int k = 0; k <= steps; ++k) {
    const double f = f0 + k * df;
    for (int g = 2; g <= 4; ++g) {
      const BackgroundQuadrature q{kRect, 8, 4, g};
      try {
        const auto mls = cloud(5, 3, {weight, f, degree});
        const auto c = assemble_cartesian(mls, q);
        const auto sc = solve_generalized_symmetric(c.stiffness, c.mass);
        for (Eigen::Index j = 0; j < std::min<Eigen::Index>(20, sc.eigenvalues.size()); ++j)
          out << "cartesian," << mls_weight_name(weight) << ',' << degree << ',' << f << ','
              << g << ',' << j + 1 << ',' << sc.eigenvalues[j] << ",ok\n";
        const auto r = assemble_rot_irrot(mls, q);
        const auto sr = solve_rot_irrot_rank_reduced(r);
        for (Eigen::Index j = 0; j < std::min<Eigen::Index>(20, sr.eigenvalues.size()); ++j)
          out << "rot_irrot," << mls_weight_name(weight) << ',' << degree << ',' << f << ','
              << g << ',' << j + 1 << ',' << sr.eigenvalues[j] << ",ok\n";
      } catch (const std::exception& e) {
        out << "failed," << mls_weight_name(weight) << ',' << degree << ',' << f << ',' << g
            << ",0,,\"" << e.what() << "\"\n";
      }
    }
  }
}

// --------------------------------------------------------------- laplacian
// Reproduction defect of the MLS Laplacian for q=x^2+y^2 (exact value 4)
// and of the gradient of q, at all points of a fine rule. Degree 1 has no
// second-order consistency; degree 2 must reproduce it to round-off.
void laplacian(const std::vector<std::string>& a) {
  if (a.size() != 1) throw std::invalid_argument("laplacian out.csv");
  auto out = open_csv(a[0]);
  out << "weight,degree,factor,nodes_y,nodes_x,max_moment_condition,laplacian_defect_max,"
         "laplacian_defect_rms,interior_laplacian_defect_rms,gradient_defect_rms,status\n";
  const std::vector<Profile> profiles{
      {MLSWeight::QuarticSpline, 1.5, 1}, {MLSWeight::QuarticSpline, 2.5, 1},
      {MLSWeight::QuarticSpline, 3.3, 1}, {MLSWeight::WendlandC2, 3.3, 1},
      {MLSWeight::QuarticSpline, 2.5, 2}, {MLSWeight::QuarticSpline, 3.3, 2},
      {MLSWeight::WendlandC2, 3.3, 2}};
  const std::vector<std::pair<int, int>> grids{{5, 3}, {9, 5}, {17, 9}, {33, 17}};
  for (const auto& p : profiles)
    for (const auto& g : grids) {
      out << mls_weight_name(p.weight) << ',' << p.degree << ',' << p.factor << ','
          << g.second << ',' << g.first << ',';
      try {
        const auto mls = cloud(g.first, g.second, p);
        const BackgroundQuadrature q{kRect, (g.first - 1) * 2, (g.second - 1) * 2, 4};
        const double margin = p.factor * 2.0 / (g.first - 1);
        double cmax = 0, dmax = 0, d2 = 0, area = 0, i2 = 0, iarea = 0, g2 = 0;
        for_each_quadrature_point(q, [&](Point2 x, double w) {
          const auto s = mls.evaluate(x);
          cmax = std::max(cmax, s.condition_number);
          double lap = 0, gx = 0, gy = 0;
          for (std::size_t k = 0; k < s.active_nodes.size(); ++k) {
            const Point2 n = mls.nodes()[s.active_nodes[k]].position;
            const double v = n.x * n.x + n.y * n.y;
            lap += v * (s.dxx[k] + s.dyy[k]);
            gx += v * s.dx[k];
            gy += v * s.dy[k];
          }
          const double d = lap - 4.0;
          dmax = std::max(dmax, std::abs(d));
          d2 += w * d * d; area += w;
          g2 += w * ((gx - 2 * x.x) * (gx - 2 * x.x) + (gy - 2 * x.y) * (gy - 2 * x.y));
          if (x.x > margin && x.x < 2 - margin && x.y > margin && x.y < 1 - margin) {
            i2 += w * d * d; iarea += w;
          }
        });
        out << cmax << ',' << dmax << ',' << std::sqrt(d2 / area) << ',';
        if (iarea > 0) out << std::sqrt(i2 / iarea);
        out << ',' << std::sqrt(g2 / area) << ",ok\n";
      } catch (const std::exception& e) {
        out << ",,,,,\"" << e.what() << "\"\n";
      }
    }
}

// ------------------------------------------------------------- convergence
// M5 TM11/TM21 with the identified author profile and alternatives, at a
// production rule (4 cells/rectangle, Gauss g) and a converged rule
// (16 cells/rectangle, Gauss 5), separating quadrature from approximation.
void convergence(const std::vector<std::string>& a) {
  if (a.size() != 1) throw std::invalid_argument("convergence out.csv");
  auto out = open_csv(a[0]);
  out << "formulation,weight,degree,factor,nodes_y,nodes_x,subdivision,gauss,mode,"
         "best_curl_lambda,best_curl_signed_error_pct,best_curl_rho,"
         "nearest_lambda,nearest_signed_error_pct,positive_below_0p98_TM11,"
         "retained_rank,mass_nullity,curl_zero_count,mass_kernel_stiffness_coupling,"
         "max_full_residual\n";
  const std::vector<Profile> profiles{
      {MLSWeight::QuarticSpline, 1.5, 1}, {MLSWeight::QuarticSpline, 2.5, 1},
      {MLSWeight::QuarticSpline, 3.3, 1}, {MLSWeight::WendlandC2, 3.3, 1},
      {MLSWeight::QuarticSpline, 2.5, 2}, {MLSWeight::WendlandC2, 3.3, 2}};
  const std::vector<std::pair<int, int>> grids{{5, 3}, {9, 5}, {13, 7}, {17, 9}, {25, 13}};
  const std::vector<std::pair<int, int>> rules{{2, 1}, {2, 2}, {2, 4}, {4, 5}};
  const double exact[2] = {1.25 * kPi * kPi, 2 * kPi * kPi};
  for (const auto& p : profiles)
    for (const auto& g : grids) {
      const auto mls = cloud(g.first, g.second, p);
      const BackgroundQuadrature ref{kRect, (g.first - 1) * 4, (g.second - 1) * 4, 5};
      const CurlForms forms_ref = curl_forms(mls, ref);
      for (const auto& r : rules) {
        const BackgroundQuadrature q{kRect, (g.first - 1) * r.first,
                                     (g.second - 1) * r.first, r.second};
        const auto forms = assemble_rot_irrot(mls, q);
        const auto s = solve_rot_irrot_rank_reduced(forms);
        int below = 0;
        for (Eigen::Index j = 0; j < s.eigenvalues.size(); ++j)
          if (s.eigenvalues[j] > s.eigenvalue_zero_cutoff && s.eigenvalues[j] < 0.98 * exact[0])
            ++below;
        for (int t = 0; t < 2; ++t) {
          const Eigen::VectorXd rho = m5_curl_rho(forms_ref, s, t);
          Eigen::Index best = -1, nearest = -1;
          for (Eigen::Index j = 0; j < s.eigenvalues.size(); ++j) {
            if (s.eigenvalues[j] <= s.eigenvalue_zero_cutoff) continue;
            if (best < 0 || rho[j] > rho[best]) best = j;
            if (nearest < 0 || std::abs(s.eigenvalues[j] - exact[t]) <
                                   std::abs(s.eigenvalues[nearest] - exact[t]))
              nearest = j;
          }
          out << "m5," << mls_weight_name(p.weight) << ',' << p.degree << ',' << p.factor << ','
              << g.second << ',' << g.first << ',' << r.first << ',' << r.second
              << ",TM" << t + 1 << "1," << s.eigenvalues[best] << ','
              << 100 * (s.eigenvalues[best] / exact[t] - 1) << ',' << rho[best] << ','
              << s.eigenvalues[nearest] << ',' << 100 * (s.eigenvalues[nearest] / exact[t] - 1)
              << ',' << below << ',' << s.retained_rank << ',' << s.removed_mass_nullity << ','
              << s.numerical_curl_kernel << ',' << s.mass_kernel_stiffness_coupling << ','
              << s.max_full_residual << "\n";
        }
      }
      std::cerr << "convergence " << mls_weight_name(p.weight) << " f=" << p.factor
                << " deg=" << p.degree << ' ' << g.second << 'x' << g.first << " done\n";
    }
}

// -------------------------------------------------------------- rank-sweep
void rank_sweep(const std::vector<std::string>& a) {
  if (a.size() != 1) throw std::invalid_argument("rank-sweep out.csv");
  auto out = open_csv(a[0]);
  out << "weight,degree,factor,nodes_y,nodes_x,gauss,relative_cutoff,status,retained_rank,"
         "mass_nullity,mass_eigen_max_removed_rel,mass_eigen_min_retained_rel,curl_zero_count,"
         "first_positive_lambda,max_lambda,mass_kernel_stiffness_coupling,"
         "max_reduced_residual,max_full_residual\n";
  struct Case { int nx, ny, g; Profile p; };
  const std::vector<Case> cases{
      {5, 3, 2, {MLSWeight::WendlandC2, 3.3, 1}}, {5, 3, 2, {MLSWeight::QuarticSpline, 1.5, 1}},
      {9, 5, 2, {MLSWeight::WendlandC2, 3.3, 1}}, {9, 5, 2, {MLSWeight::WendlandC2, 3.3, 2}},
      {13, 7, 1, {MLSWeight::WendlandC2, 3.3, 1}}, {13, 7, 2, {MLSWeight::QuarticSpline, 1.5, 1}}};
  for (const auto& c : cases) {
    const auto mls = cloud(c.nx, c.ny, c.p);
    const BackgroundQuadrature q{kRect, (c.nx - 1) * 2, (c.ny - 1) * 2, c.g};
    const auto forms = assemble_rot_irrot(mls, q);
    for (double tol : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-15}) {
      out << mls_weight_name(c.p.weight) << ',' << c.p.degree << ',' << c.p.factor << ','
          << c.ny << ',' << c.nx << ',' << c.g << ',' << tol << ',';
      try {
        const auto s = solve_rot_irrot_rank_reduced(forms, tol, 0.0);
        double first = 0;
        for (Eigen::Index j = 0; j < s.eigenvalues.size(); ++j)
          if (s.eigenvalues[j] > s.eigenvalue_zero_cutoff) { first = s.eigenvalues[j]; break; }
        out << "ok," << s.retained_rank << ',' << s.removed_mass_nullity << ','
            << s.mass_eigen_max_removed / s.mass_eigen_max << ','
            << s.mass_eigen_min_retained / s.mass_eigen_max << ',' << s.numerical_curl_kernel
            << ',' << first << ',' << s.eigenvalues.maxCoeff() << ','
            << s.mass_kernel_stiffness_coupling << ',' << s.max_reduced_residual << ','
            << s.max_full_residual << "\n";
      } catch (const std::exception& e) {
        out << "\"failed: " << e.what() << "\",,,,,,,,,,\n";
      }
    }
  }
}

// ----------------------------------------------------------------- doublet
// Cartesian TM22/TM41 doublet (M6a/M6b setup). Compares the computed
// eigen-subspace, the L2-best approximation of the analytic doublet in the
// SAME MLS space, and the best 2D subspace of a spectral window.
struct AnalyticField { int m, n; };

double rayleigh_max(const Eigen::MatrixXd& V, const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B, double* min_out) {
  Eigen::MatrixXd a = V.transpose() * A * V, b = V.transpose() * B * V;
  a = 0.5 * (a + a.transpose()); b = 0.5 * (b + b.transpose());
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(a, b);
  if (es.info() != Eigen::Success) throw std::runtime_error("Rayleigh range failed");
  if (min_out) *min_out = std::max(0.0, es.eigenvalues().minCoeff());
  return std::max(0.0, es.eigenvalues().maxCoeff());
}

void doublet(const std::vector<std::string>& args) {
  if (args.size() != 3) throw std::invalid_argument("doublet nx ny out.csv");
  const int nx = std::stoi(args[0]), ny = std::stoi(args[1]);
  const auto mls = cloud(nx, ny, {MLSWeight::WendlandC2, 3.3, 1});
  const BackgroundQuadrature q{kRect, nx - 1, ny - 1, 3};
  const auto cart = assemble_cartesian(mls, q);
  const auto cons = assemble_constraints(mls, q, 3);
  const auto spec = solve_generalized_symmetric(cart.stiffness, cart.mass);
  const auto pairing = assemble_gradient_pairing(mls, q);
  const auto gp = build_gradient_projection(cart.mass, pairing);
  const Eigen::MatrixXd M(cart.mass), K(cart.stiffness), D(cons.divergence),
      Bn(cons.normal_boundary);
  const Eigen::MatrixXd G = M * gp.basis;          // P_g form: M Q Q^T M
  const Eigen::MatrixXd Pg = G * G.transpose();
  const auto n = static_cast<Eigen::Index>(mls.nodes().size());
  const double lam_ref = 5 * kPi * kPi;
  const std::vector<AnalyticField> doublet_fields{{2, 2}, {4, 1}};
  // Load vectors b_k = int N . H_k and exact Gram of the analytic fields.
  Eigen::MatrixXd load = Eigen::MatrixXd::Zero(2 * n, 2);
  Eigen::Matrix2d exact_gram = Eigen::Matrix2d::Zero();
  for_each_quadrature_point(q, [&](Point2 p, double w) {
    const auto s = mls.evaluate(p);
    double hx[2], hy[2];
    for (int k = 0; k < 2; ++k) {
      const double kx = doublet_fields[static_cast<std::size_t>(k)].m * kPi / 2;
      const double ky = doublet_fields[static_cast<std::size_t>(k)].n * kPi;
      hx[k] = -ky * std::sin(kx * p.x) * std::cos(ky * p.y);
      hy[k] = kx * std::cos(kx * p.x) * std::sin(ky * p.y);
    }
    for (int k = 0; k < 2; ++k)
      for (int l = 0; l < 2; ++l) exact_gram(k, l) += w * (hx[k] * hx[l] + hy[k] * hy[l]);
    for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
      const auto i = static_cast<Eigen::Index>(s.active_nodes[a]);
      for (int k = 0; k < 2; ++k) {
        load(i, k) += w * s.phi[a] * hx[k];
        load(n + i, k) += w * s.phi[a] * hy[k];
      }
    }
  });
  const Eigen::LLT<Eigen::MatrixXd> mchol(M);
  const Eigen::MatrixXd projected = mchol.solve(load);  // L2-best approximations
  // Principal cosines between a discrete subspace V and the analytic doublet.
  auto cosines_to_analytic = [&](const Eigen::MatrixXd& V) {
    Eigen::MatrixXd gv = V.transpose() * M * V;
    const Eigen::MatrixXd lv = Eigen::LLT<Eigen::MatrixXd>(0.5 * (gv + gv.transpose())).matrixL();
    const Eigen::Matrix2d la = Eigen::LLT<Eigen::Matrix2d>(exact_gram).matrixL();
    const Eigen::MatrixXd cross = lv.triangularView<Eigen::Lower>().solve(
        Eigen::MatrixXd(V.transpose() * load));
    const Eigen::MatrixXd c = la.triangularView<Eigen::Lower>()
                                  .solve(cross.transpose()).transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(c);
    return Eigen::VectorXd(svd.singularValues().cwiseMin(1.0));
  };
  auto out = open_csv(args[2]);
  out << "nx,ny,kind,indices,lambda_min,lambda_max,dimension,analytic_cos_1,analytic_cos_2,"
         "div_ratio_min,div_ratio_max,normal_ratio_min,normal_ratio_max,"
         "gradient_fraction_min,gradient_fraction_max,curl_rayleigh_min,curl_rayleigh_max\n";
  const double perimeter_scale = 2.0 / 6.0;  // area/perimeter of the 2x1 guide
  auto report = [&](const std::string& kind, const std::string& ids, const Eigen::MatrixXd& V,
                    double lmin, double lmax) {
    double dmin, nmin, gmin, kmin;
    const double dmax = rayleigh_max(V, D, M, &dmin);
    const double nmax = rayleigh_max(V, Bn, M, &nmin);
    const double gmax = rayleigh_max(V, Pg, M, &gmin);
    const double kmax = rayleigh_max(V, K, M, &kmin);
    const Eigen::VectorXd c = cosines_to_analytic(V);
    out << nx << ',' << ny << ',' << kind << ',' << ids << ',' << lmin << ',' << lmax << ','
        << V.cols() << ',' << c[0] << ',' << (c.size() > 1 ? c[1] : 0.0) << ','
        << std::sqrt(dmin / lam_ref) << ',' << std::sqrt(dmax / lam_ref) << ','
        << std::sqrt(nmin * perimeter_scale) << ',' << std::sqrt(nmax * perimeter_scale) << ','
        << gmin << ',' << gmax << ',' << kmin << ',' << kmax << "\n";
  };
  // Window of computed eigenvectors around 5 pi^2 (within +-3 percent).
  std::vector<Eigen::Index> window;
  for (Eigen::Index j = 0; j < spec.eigenvalues.size(); ++j)
    if (std::abs(spec.eigenvalues[j] / lam_ref - 1) <= 0.03) window.push_back(j);
  if (window.size() < 2) throw std::runtime_error("Fewer than two modes near 5 pi^2");
  // Individual eigenvectors (basis-dependent quantities, one per column).
  for (Eigen::Index j : window) {
    report("eigenvector", std::to_string(j), spec.eigenvectors.col(j),
           spec.eigenvalues[j], spec.eigenvalues[j]);
  }
  // The computed doublet = the two window modes closest to the analytic doublet.
  std::vector<std::pair<double, Eigen::Index>> score;
  for (Eigen::Index j : window)
    score.push_back({cosines_to_analytic(spec.eigenvectors.col(j))[0], j});
  std::sort(score.rbegin(), score.rend());
  Eigen::Index i0 = std::min(score[0].second, score[1].second),
               i1 = std::max(score[0].second, score[1].second);
  Eigen::MatrixXd pair(2 * n, 2);
  pair << spec.eigenvectors.col(i0), spec.eigenvectors.col(i1);
  report("computed_doublet_span", std::to_string(i0) + ";" + std::to_string(i1), pair,
         spec.eigenvalues[i0], spec.eigenvalues[i1]);
  report("l2_projection_of_analytic_doublet", "P(TM22);P(TM41)", projected, lam_ref, lam_ref);
  for (int k = 0; k < 2; ++k)
    report(k == 0 ? "l2_projection_TM22" : "l2_projection_TM41", "single",
           projected.col(k), lam_ref, lam_ref);
  // Best 2D subspace of span(window) w.r.t. the analytic doublet:
  // principal vectors obtained from the SVD of the M/L2 cross matrix.
  Eigen::MatrixXd W(2 * n, static_cast<Eigen::Index>(window.size()));
  for (std::size_t c = 0; c < window.size(); ++c) W.col(static_cast<Eigen::Index>(c)) =
      spec.eigenvectors.col(window[c]);
  const Eigen::Matrix2d la = Eigen::LLT<Eigen::Matrix2d>(exact_gram).matrixL();
  const Eigen::MatrixXd cross = W.transpose() * load;  // W is M-orthonormal
  const Eigen::MatrixXd c = la.triangularView<Eigen::Lower>().solve(cross.transpose()).transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(c, Eigen::ComputeThinU);
  const Eigen::MatrixXd best = W * svd.matrixU().leftCols(2);
  std::string ids;
  for (std::size_t k = 0; k < window.size(); ++k) ids += (k ? ";" : "") + std::to_string(window[k]);
  report("best_2d_subspace_of_window", ids, best, spec.eigenvalues[window.front()],
         spec.eigenvalues[window.back()]);
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::invalid_argument(
        "Usage: independent_audit ch4-scan|laplacian|convergence|rank-sweep|doublet ...");
    const std::string cmd = argv[1];
    const std::vector<std::string> args(argv + 2, argv + argc);
    if (cmd == "ch4-scan") ch4_scan(args);
    else if (cmd == "laplacian") laplacian(args);
    else if (cmd == "convergence") convergence(args);
    else if (cmd == "rank-sweep") rank_sweep(args);
    else if (cmd == "doublet") doublet(args);
    else throw std::invalid_argument("Unknown subcommand " + cmd);
  } catch (const std::exception& e) {
    std::cerr << "AUDIT ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
