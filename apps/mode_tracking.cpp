// SPDX-License-Identifier: MIT
// M7: follow eigenvectors across MLS support factors by their physical fields
// on ONE common sample rule (independent of the cloud), not by sorted index.
#include "al_maleh/mode_tracking.hpp"
#include "al_maleh/rot_irrot.hpp"
#include "al_maleh/spectral.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {

struct Target { int m, n; };
const Target kTargets[] = {{1,1},{2,1},{3,1},{1,2},{2,2},{4,1}};
constexpr int kTargetCount = 6;

struct Sampled {
  std::vector<Eigen::Index> indices;  // positions in the solver's spectrum
  Eigen::VectorXd lambda;
  Eigen::MatrixXd field;   // rows: hx at P points, then hy at P points
  Eigen::MatrixXd curl;    // rows: P points
  double zero_cutoff = 0;
  double max_full_residual = 0;
};

Sampled solve_and_sample(const std::string& formulation, const MLSEvaluator& mls,
                         const BackgroundQuadrature& assembly_rule,
                         const std::vector<Point2>& points, double lambda_max) {
  const auto n = static_cast<Eigen::Index>(mls.nodes().size());
  Eigen::VectorXd values;
  Eigen::MatrixXd vectors;
  Sampled out;
  if (formulation == "m5") {
    const auto forms = assemble_rot_irrot(mls, assembly_rule);
    const auto s = solve_rot_irrot_rank_reduced(forms);
    values = s.eigenvalues; vectors = s.eigenvectors;
    out.zero_cutoff = s.eigenvalue_zero_cutoff;
    out.max_full_residual = s.max_full_residual;
  } else {
    const auto a = assemble_cartesian(mls, assembly_rule);
    const auto s = solve_generalized_symmetric(a.stiffness, a.mass);
    values = s.eigenvalues; vectors = s.eigenvectors;
    out.zero_cutoff = s.zero_cutoff;
    out.max_full_residual = *std::max_element(s.residuals.begin(), s.residuals.end());
  }
  for (Eigen::Index j = 0; j < values.size(); ++j)
    if (values[j] > out.zero_cutoff && values[j] <= lambda_max) out.indices.push_back(j);
  const auto k = static_cast<Eigen::Index>(out.indices.size());
  if (k < 1) throw std::runtime_error("No positive mode below lambda_max");
  Eigen::MatrixXd V(vectors.rows(), k);
  out.lambda.resize(k);
  for (Eigen::Index c = 0; c < k; ++c) {
    V.col(c) = vectors.col(out.indices[static_cast<std::size_t>(c)]);
    out.lambda[c] = values[out.indices[static_cast<std::size_t>(c)]];
  }
  const auto P = static_cast<Eigen::Index>(points.size());
  out.field = Eigen::MatrixXd::Zero(2 * P, k);
  out.curl = Eigen::MatrixXd::Zero(P, k);
  for (Eigen::Index q = 0; q < P; ++q) {
    const auto s = mls.evaluate(points[static_cast<std::size_t>(q)]);
    for (std::size_t a = 0; a < s.active_nodes.size(); ++a) {
      const auto i = static_cast<Eigen::Index>(s.active_nodes[a]);
      if (formulation == "m5") {
        const double lap = s.dxx[a] + s.dyy[a];
        // H = sum phi_i grad N_i + c_i z x grad N_i ; curl H = sum c_i lap N_i
        out.field.row(q).noalias() += s.dx[a] * V.row(i) - s.dy[a] * V.row(n + i);
        out.field.row(P + q).noalias() += s.dy[a] * V.row(i) + s.dx[a] * V.row(n + i);
        out.curl.row(q).noalias() += lap * V.row(n + i);
      } else {
        out.field.row(q).noalias() += s.phi[a] * V.row(i);
        out.field.row(P + q).noalias() += s.phi[a] * V.row(n + i);
        out.curl.row(q).noalias() += s.dx[a] * V.row(n + i) - s.dy[a] * V.row(i);
      }
    }
  }
  if (!out.field.allFinite() || !out.curl.allFinite())
    throw std::runtime_error("Nonfinite sampled modal field");
  return out;
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 11)
      throw std::invalid_argument(
          "Usage: mode_tracking m2|m5 nx ny gauss subdivision degree "
          "wendland_c2|quartic_spline lambda_max output_dir factor_1 [factor_2 ...]");
    const std::string formulation = argv[1];
    const int nx = std::stoi(argv[2]), ny = std::stoi(argv[3]);
    const int gauss = std::stoi(argv[4]), sub = std::stoi(argv[5]);
    const int degree = std::stoi(argv[6]);
    const MLSWeight weight = parse_mls_weight(argv[7]);
    const double lambda_max = std::stod(argv[8]);
    const std::filesystem::path out = argv[9];
    std::vector<double> factors;
    for (int k = 10; k < argc; ++k) factors.push_back(std::stod(argv[k]));
    if ((formulation != "m2" && formulation != "m5") || nx < 3 || ny < 3 ||
        nx * ny > 600 || gauss < 1 || gauss > 5 || sub < 1 || sub > 4 ||
        (degree != 1 && degree != 2) || !std::isfinite(lambda_max) || lambda_max <= 0)
      throw std::invalid_argument("Invalid mode-tracking parameters");
    for (double f : factors)
      if (!std::isfinite(f) || f <= 0) throw std::invalid_argument("Invalid support factor");

    const Rectangle rect{0, 2, 0, 1};
    const double h = std::max(2.0 / (nx - 1), 1.0 / (ny - 1));
    const BackgroundQuadrature assembly{rect, (nx - 1) * sub, (ny - 1) * sub, gauss};
    // Common spatial rule, identical for every factor (and every cloud).
    const BackgroundQuadrature common{rect, 32, 16, 3};
    std::vector<Point2> points;
    std::vector<double> pw;
    for_each_quadrature_point(common, [&](Point2 p, double w) {
      points.push_back(p); pw.push_back(w);
    });
    const auto P = static_cast<Eigen::Index>(points.size());
    Eigen::VectorXd w_curl(P), w_field(2 * P);
    for (Eigen::Index q = 0; q < P; ++q) {
      w_curl[q] = pw[static_cast<std::size_t>(q)];
      w_field[q] = w_field[P + q] = pw[static_cast<std::size_t>(q)];
    }
    const double pi = std::acos(-1.0);
    Eigen::MatrixXd ref_field(2 * P, kTargetCount), ref_curl(P, kTargetCount);
    Eigen::VectorXd ref_lambda(kTargetCount);
    for (int t = 0; t < kTargetCount; ++t) {
      const double kx = kTargets[t].m * pi / 2.0, ky = kTargets[t].n * pi;
      ref_lambda[t] = kx * kx + ky * ky;
      for (Eigen::Index q = 0; q < P; ++q) {
        const Point2 p = points[static_cast<std::size_t>(q)];
        const double psi = std::sin(kx * p.x) * std::sin(ky * p.y);
        ref_curl(q, t) = psi;  // curl(z x grad psi) = lap psi = -k^2 psi
        ref_field(q, t) = -ky * std::sin(kx * p.x) * std::cos(ky * p.y);
        ref_field(P + q, t) = kx * std::cos(kx * p.x) * std::sin(ky * p.y);
      }
    }

    std::filesystem::create_directories(out);
    std::ofstream spectrum(out / "spectrum.csv"), tracking(out / "tracking.csv"),
        selection(out / "selection.csv");
    if (!spectrum || !tracking || !selection)
      throw std::runtime_error("Cannot create mode-tracking outputs");
    for (auto* s : {&spectrum, &tracking, &selection}) *s << std::setprecision(17);
    spectrum << "factor,index,lambda";
    for (const auto& t : kTargets)
      spectrum << ",h_overlap_TM" << t.m << t.n << ",curl_rho_TM" << t.m << t.n;
    spectrum << "\n";
    tracking << "factor_from,factor_to,index_from,lambda_from,index_to,lambda_to,"
                "h_overlap,runner_up_overlap,status\n";
    selection << "factor,target,lambda_exact,best_curl_index,best_curl_lambda,"
                 "best_curl_signed_error_pct,best_curl_rho,nearest_index,nearest_lambda,"
                 "nearest_signed_error_pct,best_curl_from_tracked_previous_best,"
                 "two_best_curl_span_cosine_with_reference_H,max_full_residual\n";

    Sampled previous;
    std::vector<int> previous_best(kTargetCount, -1);
    for (std::size_t step = 0; step < factors.size(); ++step) {
      const double f = factors[step];
      const MLSEvaluator mls(make_regular_cloud(rect, nx, ny, f * h), degree, 1e12, weight);
      const Sampled cur = solve_and_sample(formulation, mls, assembly, points, lambda_max);
      const Eigen::MatrixXd hover = normalized_overlaps(cur.field, ref_field, w_field);
      const Eigen::MatrixXd crho = normalized_overlaps(cur.curl, ref_curl, w_curl);
      for (Eigen::Index c = 0; c < cur.lambda.size(); ++c) {
        spectrum << f << ',' << cur.indices[static_cast<std::size_t>(c)] << ',' << cur.lambda[c];
        for (int t = 0; t < kTargetCount; ++t)
          spectrum << ',' << hover(c, t) << ',' << crho(c, t);
        spectrum << "\n";
      }
      std::vector<int> map_from_previous;
      if (step > 0) {
        const Eigen::MatrixXd O = normalized_overlaps(previous.field, cur.field, w_field);
        map_from_previous = max_overlap_assignment(O);
        for (Eigen::Index r = 0; r < O.rows(); ++r) {
          const int to = map_from_previous[static_cast<std::size_t>(r)];
          double runner = 0;
          for (Eigen::Index c = 0; c < O.cols(); ++c)
            if (c != to) runner = std::max(runner, O(r, c));
          const double ov = to >= 0 ? O(r, to) : 0.0;
          const char* status = to < 0 ? "unmatched"
                              : (ov >= 0.9 && ov - runner >= 0.2 ? "tracked" : "ambiguous");
          tracking << factors[step - 1] << ',' << f << ','
                   << previous.indices[static_cast<std::size_t>(r)] << ',' << previous.lambda[r] << ',';
          if (to >= 0)
            tracking << cur.indices[static_cast<std::size_t>(to)] << ',' << cur.lambda[to];
          else
            tracking << ',';
          tracking << ',' << ov << ',' << runner << ',' << status << "\n";
        }
      }
      for (int t = 0; t < kTargetCount; ++t) {
        Eigen::Index best = 0, nearest = 0, second = -1;
        for (Eigen::Index c = 0; c < cur.lambda.size(); ++c) {
          if (crho(c, t) > crho(best, t)) best = c;
          if (std::abs(cur.lambda[c] - ref_lambda[t]) <
              std::abs(cur.lambda[nearest] - ref_lambda[t])) nearest = c;
        }
        for (Eigen::Index c = 0; c < cur.lambda.size(); ++c)
          if (c != best && (second < 0 || crho(c, t) > crho(second, t))) second = c;
        double span_cos = hover(best, t);
        if (second >= 0) {
          Eigen::MatrixXd pair(2 * P, 2);
          pair << cur.field.col(best), cur.field.col(second);
          span_cos = principal_cosines(ref_field.col(t), pair, w_field)[0];
        }
        std::string tracked = "";
        if (step > 0 && previous_best[static_cast<std::size_t>(t)] >= 0)
          tracked = map_from_previous[static_cast<std::size_t>(
                        previous_best[static_cast<std::size_t>(t)])] == best ? "1" : "0";
        selection << f << ",TM" << kTargets[t].m << kTargets[t].n << ',' << ref_lambda[t] << ','
                  << cur.indices[static_cast<std::size_t>(best)] << ',' << cur.lambda[best] << ','
                  << 100 * (cur.lambda[best] / ref_lambda[t] - 1) << ',' << crho(best, t) << ','
                  << cur.indices[static_cast<std::size_t>(nearest)] << ',' << cur.lambda[nearest] << ','
                  << 100 * (cur.lambda[nearest] / ref_lambda[t] - 1) << ',' << tracked << ','
                  << span_cos << ',' << cur.max_full_residual << "\n";
        previous_best[static_cast<std::size_t>(t)] = static_cast<int>(best);
      }
      previous = cur;
      std::cout << "M7 " << formulation << ' ' << ny << 'x' << nx << " factor=" << f
                << " modes<=" << lambda_max << ": " << cur.lambda.size() << "\n";
    }
    std::ofstream meta(out / "run_metadata.txt");
    meta << std::setprecision(17) << "formulation=" << formulation << "\nnx=" << nx
         << "\nny=" << ny << "\ngauss=" << gauss << "\nsubdivision=" << sub
         << "\ndegree=" << degree << "\nweight=" << mls_weight_name(weight)
         << "\nlambda_max=" << lambda_max
         << "\ncommon_rule=rectangle_2x1_cells_32x16_gauss3_points=" << P
         << "\ntracking_norm=L2_of_H_on_common_rule"
         << "\nassignment=Hungarian_max_total_overlap"
         << "\nstatus_rule=tracked_if_overlap>=0.9_and_margin>=0.2_else_ambiguous\n";
  } catch (const std::exception& e) {
    std::cerr << "M7 ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
