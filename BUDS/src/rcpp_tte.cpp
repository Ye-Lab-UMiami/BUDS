// rcpp_tte.cpp — TTE design engine
// Both single-point and BUDS use the paper's iterative zoom-in grid search.
// BUDS: generate local feasible candidates -> evaluate analytically over the
// interval -> robust filter -> pool -> select by interval BUDS objective.

#include <Rcpp.h>
#include "tte_math.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <cctype>

using namespace Rcpp;

// ============================================================
// Core: iterative zoom-in grid search for a single (hz0, hz1, hr)
// This is a direct C++ translation of the paper's R algorithm.
// Returns the optimal (n1, c1, n2, c2, DA1, DA2, EN_H0, PET).
// ============================================================
struct TTEResult {
    int n1, n2;
    double c1, c2, DA1, DA2, EN_H0, PET;
    int iters;
    bool valid;
};

static TTEResult tte_zoom_search(
    double hz0, double hz1, double hr, double x0,
    double alpha, double beta, double rate,
    int nsingle, double tasingle,
    double ceps, double alphaeps, int nbmaxiter)
{
    TTEResult res = {0,0,0,0,0,0,0,0,0,false};

    // Start from a reasonable center: the corresponding single-stage design
    // gives a practical scale for the two-stage optimum.
    double best_n = (double)nsingle;
    double best_DA1 = tasingle;
    double best_c1 = 0.25;

    // nbpt: number of grid points per dimension in the local search window.
    int nbpt = 11;
    // pascote: window shrink factor carried over from the existing
    // twostage_tte.R implementation in this project.
    // Each iteration reduces the search radius by about 1 / 1.26 = 0.79,
    // which is a mild contraction: fast enough to converge, but not so
    // aggressive that we easily jump over a good neighborhood.
    double pascote = 1.26;
    // cote: current half-width scale of the local search window.
    double cote = pascote;
    // c1lo/c1hi track the current stage-1 cutoff search interval.
    double c1lo = -nbpt, c1hi = nbpt;
    // EnH0: best expected sample size under H0 found so far.
    double EnH0 = 1e6;
    int iter = 0;

    while (iter < nbmaxiter && (c1hi - c1lo) / nbpt > 0.001) {
        iter++;
        if (iter % 2 == 0) nbpt++;
        cote /= pascote;

        // Build a smaller box centered at the current incumbent.
        // nlo/nhi: total sample-size window
        // dlo/dhi: stage-1 analysis time window
        // c1lo/c1hi: stage-1 critical value window
        double nlo = best_n - nsingle * cote;
        double nhi = best_n + nsingle * cote;
        double dlo = std::max(best_DA1 - tasingle * cote, 1e-4);
        double dhi = best_DA1 + tasingle * cote;
        c1lo = best_c1 - cote;
        c1hi = best_c1 + cote;

        // Build local grid vectors (like expand.grid in R).
        // tvec is expressed in calendar time ta = n / rate.
        std::vector<double> tvec, dvec, cvec;
        for (int i = 0; i < nbpt; i++) {
            int nc = (int)std::ceil(nlo + (nhi - nlo) * i / (nbpt - 1.0));
            if (nc > 0 && (tvec.empty() || tvec.back() != nc / rate))
                tvec.push_back(nc / rate);
        }
        for (int i = 0; i < nbpt; i++) {
            double v = dlo + (dhi - dlo) * i / (nbpt - 1.0);
            if (v > 0) dvec.push_back(v);
        }
        for (int i = 0; i < nbpt; i++)
            cvec.push_back(c1lo + (c1hi - c1lo) * i / (nbpt - 1.0));

        // Best candidate inside the current zoom box.
        double le = EnH0, ln = best_n, ld = best_DA1, lc = best_c1;

        // Evaluate all grid combinations (like apply in R)
        for (auto ta : tvec)
          for (auto DA1 : dvec) {
            if (ta <= DA1) continue;
            for (auto c1v : cvec) {
                // pap: probability of continuing past stage 1 under H0.
                double pap = 1.0 - R::pnorm(c1v, 0, 1, 1, 0);
                // enh0: expected sample size under H0 for this candidate.
                double enh0 = (ta - std::max(0.0, ta - DA1) * pap) * rate;
                if (enh0 > EnH0) continue;

                // Variance / covariance components from Appendix A.
                // hzb: midpoint hazard used in local power approximation.
                double hzb = (hz0 + hz1) / 2.0;
                // B: stage-1 follow-up is restricted by the design window x0.
                double B = (DA1 < x0) ? DA1 : x0;
                // v1, v: null-information terms for stage 1 and final analysis.
                double v1 = hz0 * tte::fct(DA1, ta, hz0, B);
                double v  = 1.0 - std::exp(-hz0 * x0);
                if (v <= 0 || v1 <= 0) continue;

                // s11, s1, s01, s0: alternative-information terms used to
                // evaluate power and the stage-1 / stage-2 correlation.
                double s11 = hzb * tte::fct(DA1, ta, hzb, B);
                double s1  = 1.0 - std::exp(-hzb * x0);
                double s01 = hz0 * tte::fct(DA1, ta, hz1, B);
                double s0  = hr * (1.0 - std::exp(-hz1 * x0));
                // om1, om: mean shifts under H1 at stage 1 and final analysis.
                double om1 = (hz1 - hz0) * tte::fct(DA1, ta, hz1, B);
                double om  = (1.0 - hr) * (1.0 - std::exp(-hz1 * x0));
                if (s11 <= 0 || s1 <= 0 || s01 <= 0 || s0 <= 0) continue;

                // rho0: correlation between stage-1 and final Z statistics
                // under H0.
                double rho0 = std::sqrt(v1 / v);

                // Solve for c2 by bisection so that type-I error is alpha.
                // cd/cf are the lower/upper brackets; ac is the achieved alpha.
                double cd = -10, cf = 10, ac = 100;
                int it2 = 0;
                while ((std::abs(ac - alpha) > alphaeps || cf - cd > ceps) && it2 < nbmaxiter) {
                    it2++;
                    double cc = (cd + cf) / 2.0;
                    ac = tte::bvn_integral(cc, c1v, rho0);
                    if (ac > alpha) cf = cc; else cd = cc;
                }
                if (cf - cd > ceps * 10) continue;
                double c2v = (cd + cf) / 2.0;

                // rho1: stage-1/final correlation under the working alternative.
                double rho1 = std::sqrt(s11 / s1);
                // cb1/cb: transformed critical values used in the alternative
                // bivariate normal power calculation.
                double cb1 = std::sqrt(s01 / s11) * (c1v - om1 * std::sqrt(rate * DA1) / std::sqrt(s01));
                double cb  = std::sqrt(s0 / s1) * (c2v - om * std::sqrt(rate * ta) / std::sqrt(s0));
                if (tte::bvn_integral(cb, cb1, rho1) < 1.0 - beta) continue;

                if (enh0 < le) {
                    le = enh0; ln = ta * rate; ld = DA1; lc = c1v;
                }
            }
          }

        if (le < EnH0) {
            EnH0 = le; best_n = ln; best_DA1 = ld; best_c1 = lc;
        }
    }

    // Final pass: recompute c2 for the best point found in the last zoom box.
    double bta = best_n / rate;
    double B2 = (best_DA1 < x0) ? best_DA1 : x0;
    double v1 = hz0 * tte::fct(best_DA1, bta, hz0, B2);
    double v  = 1.0 - std::exp(-hz0 * x0);
    if (v <= 0 || v1 <= 0) return res;
    double rho0 = std::sqrt(v1 / v);
    double cd = -10, cf = 10, ac = 100;
    int it2 = 0;
    while ((std::abs(ac - alpha) > alphaeps || cf - cd > ceps) && it2 < nbmaxiter) {
        it2++;
        double cc = (cd + cf) / 2.0;
        ac = tte::bvn_integral(cc, best_c1, rho0);
        if (ac > alpha) cf = cc; else cd = cc;
    }

    res.n1 = (int)std::ceil(best_DA1 * rate);
    res.n2 = (int)std::round(best_n);
    res.c1 = best_c1;
    res.c2 = (cd + cf) / 2.0;
    res.DA1 = best_DA1;
    res.DA2 = bta + x0;
    res.EN_H0 = EnH0;
    res.PET = 1.0 - R::pnorm(best_c1, 0, 1, 1, 0);
    res.iters = iter;
    res.valid = true;
    return res;
}

// Single-stage sample size (closed-form)
static int tte_single_n(double hz0, double hz1, double hr, double x0,
                         double alpha, double beta) {
    double za = R::qnorm(1 - alpha, 0, 1, 1, 0);
    double zb = R::qnorm(1 - beta, 0, 1, 1, 0);
    double s0 = std::sqrt(hr * (1 - std::exp(-hz1 * x0)));
    double hzb = (hz0 + hz1) / 2.0;
    double s1 = std::sqrt(1 - std::exp(-hzb * x0));
    double om = (1 - std::exp(-hz1 * x0)) * (1 - hr);
    if (std::abs(om) < 1e-15) return 9999;
    return (int)std::ceil(std::pow(s0 * za + s1 * zb, 2) / (om * om));
}


// ============================================================
// Exported: single-point grid search
// ============================================================
// [[Rcpp::export]]
List tte_grid_search_cpp(
    double hz0, double hz1, double hr, double x0,
    double alpha, double beta, double rate,
    int nsingle, double tasingle,
    double ceps, double alphaeps, int nbmaxiter,
    bool verbose)
{
    // Just call the zoom search with verbose printing
    if (verbose) Rprintf("  Starting zoom-in grid search (nsingle=%d)...\n", nsingle);

    TTEResult r = tte_zoom_search(hz0, hz1, hr, x0, alpha, beta, rate,
                                   nsingle, tasingle, ceps, alphaeps, nbmaxiter);

    if (verbose) Rprintf("  Converged in %d iterations.\n", r.iters);

    return List::create(
        Named("n1") = r.n1, Named("c1") = r.c1,
        Named("n2") = r.n2, Named("c2") = r.c2,
        Named("DA1") = r.DA1, Named("DA2") = r.DA2,
        Named("EN_H0") = r.EN_H0, Named("PET") = r.PET,
        Named("iters") = r.iters);
}


// ============================================================
// BUDS SEARCH
// Step 1: Generate locally feasible candidates at each design point
// Step 2: Evaluate each candidate analytically across the null interval
// Step 3: Filter candidates by robust type-I error control
// Step 4: Build the local-KJ benchmark at each evaluation point
// Step 5: Pool robust candidates across design points
// Step 6: Select final winners by interval BUDS objectives
// ============================================================

struct RKJCand {
    int source_idx;
    double p0_d;
    int n1, n2;
    double c1, c2, DA1, DA2, en, pet, lambda0_d;
};

struct LocalBenchmark {
    double p0_eval;
    double lambda0_eval;
    int n1;
    int n2;
    double c1;
    double c2;
    double DA1;
    double DA2;
    double EN_local;
    bool valid;
};

static LocalBenchmark invalid_local_benchmark(double p0_eval, double lambda0_eval) {
    LocalBenchmark lb = {
        p0_eval, lambda0_eval, NA_INTEGER, NA_INTEGER,
        NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, false
    };
    return lb;
}

static LocalBenchmark benchmark_from_zoom_result(
    double p0_eval,
    double lambda0_eval,
    const TTEResult& r)
{
    LocalBenchmark lb = invalid_local_benchmark(p0_eval, lambda0_eval);
    if (!r.valid) return lb;
    lb.n1 = r.n1;
    lb.n2 = r.n2;
    lb.c1 = r.c1;
    lb.c2 = r.c2;
    lb.DA1 = r.DA1;
    lb.DA2 = r.DA2;
    lb.EN_local = r.n1 + R::pnorm(r.c1, 0, 1, 1, 0) * (r.n2 - r.n1);
    lb.valid = true;
    return lb;
}

static bool is_same_candidate(const RKJCand& a, const RKJCand& b) {
    return a.n1 == b.n1 && a.n2 == b.n2 &&
           std::abs(a.c1 - b.c1) < 0.005 &&
           std::abs(a.c2 - b.c2) < 0.005 &&
           std::abs(a.DA1 - b.DA1) < 0.005 &&
           std::abs(a.DA2 - b.DA2) < 0.005 &&
           std::abs(a.lambda0_d - b.lambda0_d) < 1e-10;
}

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return (char)std::tolower(c);
    });
    return s;
}

static double stage_integral(double lambda_eval, double DAk, double ta, double x0) {
    if (lambda_eval <= 0.0) return NA_REAL;
    double B = std::min(DAk, x0);
    if (std::abs(DAk - (ta + x0)) < 1e-10) {
        return (1.0 - std::exp(-lambda_eval * x0)) / lambda_eval;
    }
    return tte::fct(DAk, ta, lambda_eval, B);
}

struct IntervalEvalResult {
    double alpha;
    double pet;
    double en;
    bool valid;
};

static IntervalEvalResult evaluate_candidate_interval(
    const RKJCand& cand,
    double lambda_eval,
    double x0,
    double rate)
{
    IntervalEvalResult out = {NA_REAL, NA_REAL, NA_REAL, false};

    double ta = (double)cand.n2 / rate;
    double lambda_d = cand.lambda0_d;

    double int1 = stage_integral(lambda_eval, cand.DA1, ta, x0);
    double int2 = stage_integral(lambda_eval, cand.DA2, ta, x0);
    if (!R_finite(int1) || !R_finite(int2) || int1 <= 0.0 || int2 <= 0.0) {
        return out;
    }

    double sigma01_sq = lambda_d * int1;
    double sigma02_sq = lambda_d * int2;
    double sigmae1_sq = lambda_eval * int1;
    double sigmae2_sq = lambda_eval * int2;
    if (sigma01_sq <= 0.0 || sigma02_sq <= 0.0 || sigmae1_sq <= 0.0 || sigmae2_sq <= 0.0) {
        return out;
    }

    double sigma01 = std::sqrt(sigma01_sq);
    double sigma02 = std::sqrt(sigma02_sq);
    double sigmae1 = std::sqrt(sigmae1_sq);
    double sigmae2 = std::sqrt(sigmae2_sq);
    double rho = sigmae1 / sigmae2;
    if (!R_finite(rho) || rho <= 0.0 || rho >= 1.0) {
        return out;
    }

    double omega1 = (lambda_eval - lambda_d) * int1;
    double omega2 = (lambda_eval - lambda_d) * int2;
    double cbar1 =
        (sigma01 * cand.c1 - omega1 * std::sqrt((double)cand.n1)) / sigmae1;
    double cbar2 =
        (sigma02 * cand.c2 - omega2 * std::sqrt((double)cand.n2)) / sigmae2;
    if (!R_finite(cbar1) || !R_finite(cbar2)) {
        return out;
    }

    double alpha_pt = tte::bvn_integral(cbar2, cbar1, rho);
    double pet = 1.0 - R::pnorm(cbar1, 0, 1, 1, 0);
    double en = cand.n1 + (1.0 - pet) * (cand.n2 - cand.n1);

    out.alpha = alpha_pt;
    out.pet = pet;
    out.en = en;
    out.valid = R_finite(alpha_pt) && R_finite(pet) && R_finite(en);
    return out;
}

// [[Rcpp::export]]
List tte_BUDS_search_cpp(
    double p0L, double p0U, double p1, double x0, double rate,
    double alpha, double beta,
    int n_design_points, int n_eval_grid, int nsim_eval,
    int n_ub,
    double ceps, double alphaeps, int nbmaxiter,
    Rcpp::StringVector BUDS_objective,
    Rcpp::String regret_basis,
    bool verbose)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    (void)n_eval_grid;  // interval evaluation uses the same grid as candidate generation
    (void)nsim_eval;  // kept for API compatibility; interval evaluation is closed-form

    std::vector<std::string> objectives;
    for (int i = 0; i < BUDS_objective.size(); ++i) {
        std::string obj = to_lower(Rcpp::as<std::string>(BUDS_objective[i]));
        if (obj == "avg_en" || obj == "least_regret" ||
            obj == "avg_regret" || obj == "least_en" || obj == "min_n") {
            bool seen = false;
            for (const auto& s : objectives) if (s == obj) { seen = true; break; }
            if (!seen) objectives.push_back(obj);
        }
    }
    if (objectives.empty()) objectives = {"avg_en", "least_regret"};

    std::string regret_basis_s = to_lower(std::string(regret_basis));
    if (regret_basis_s != "en" && regret_basis_s != "n")
        regret_basis_s = "en";

    // Use one grid across the null interval for all roles:
    // - candidate generation (design points)
    // - interval evaluation
    // - local KJ benchmark construction
    std::vector<double> grid(std::max(2, n_design_points));
    for (int i = 0; i < (int)grid.size(); i++)
        grid[i] = p0L + (p0U - p0L) * i / std::max(1, (int)grid.size() - 1);

    // ===== Step 1: Generate locally feasible candidate designs at each design-point p0 =====
    std::vector<RKJCand> candidates;
    candidates.reserve((int)grid.size() * 128);
    std::vector<LocalBenchmark> cached_local_bench(grid.size());

    double hz1 = -std::log(p1) / x0;

    for (int di = 0; di < (int)grid.size(); di++) {
        double p0d = grid[di];
        double hz0 = -std::log(p0d) / x0;
        double hr  = hz0 / hz1;
        cached_local_bench[di] = invalid_local_benchmark(p0d, hz0);
        if (hr <= 1.0) continue;

        int ns = tte_single_n(hz0, hz1, hr, x0, alpha, beta);
        if (ns > n_ub) {
            if (verbose) Rprintf("p0=%.3f (HR=%.3f): ns=%d > n.ub=%d, skipped\n", p0d, hr, ns, n_ub);
            continue;
        }
        double tasingle = (double)ns / rate;
        std::vector<RKJCand> source_candidates;
        source_candidates.reserve(256);

        auto add_source_candidate = [&](const RKJCand& cand) {
            for (const auto& kept : source_candidates) {
                if (is_same_candidate(kept, cand)) return;
            }
            source_candidates.push_back(cand);
        };

        // --- (A) Add the zoom-search local optimum as one feasible seed ---
        TTEResult r = tte_zoom_search(hz0, hz1, hr, x0, alpha, beta, rate,
                                       ns, tasingle, ceps, alphaeps, nbmaxiter);
        if (r.valid) {
            // Cache the local optimum on this shared-grid point so the
            // benchmark construction step can reuse the same zoom result
            // instead of recomputing it.
            cached_local_bench[di] = benchmark_from_zoom_result(p0d, hz0, r);
        }
        if (r.valid && r.n2 <= n_ub) {
            add_source_candidate({di, p0d, r.n1, r.n2, r.c1, r.c2, r.DA1, r.DA2,
                                  r.EN_H0, r.PET, hz0});
        }

        // --- (B) Systematically scan a broader feasible set for this source p0 ---
        // Keep the locally feasible source-specific set and let the
        // interval-level BUDS objectives compare designs after interval
        // evaluation and robust filtering.
        double hzb = (hz0 + hz1) / 2.0;
        int nlo = std::max(4, (int)(ns * 0.5));
        int nhi = std::min(n_ub, (int)(ns * 1.5) + 2);
        const int NC1_SCAN = 15;
        const double C1_LO = -0.5, C1_HI = 1.5;

        // Track best-min-N and best-max-PET for this p0
        RKJCand best_minN = {di, p0d, 0,99999,0,0,0,0,1e18,0,hz0};
        RKJCand best_maxPET = {di, p0d, 0,0,0,0,0,0,1e18,0,hz0};
        double max_pet_seen = -1;

        for (int n = nlo; n <= nhi; n++) {
            double ta = (double)n / rate;
            // Try a few DA1 fractions
            for (int frac : {2, 3, 4, 5, 7, 10}) {
                double DA1 = ta * frac / 15.0;
                if (DA1 <= 0 || DA1 >= ta) continue;

                for (int ic = 0; ic < NC1_SCAN; ic++) {
                    double c1v = C1_LO + (C1_HI - C1_LO) * ic / (NC1_SCAN - 1.0);

                    double B = (DA1 < x0) ? DA1 : x0;
                    double v1 = hz0 * tte::fct(DA1, ta, hz0, B);
                    double v  = 1.0 - std::exp(-hz0 * x0);
                    if (v <= 0 || v1 <= 0) continue;
                    double s11 = hzb * tte::fct(DA1, ta, hzb, B);
                    double s1  = 1.0 - std::exp(-hzb * x0);
                    double s01 = hz0 * tte::fct(DA1, ta, hz1, B);
                    double s0  = hr * (1.0 - std::exp(-hz1 * x0));
                    double om1 = (hz1 - hz0) * tte::fct(DA1, ta, hz1, B);
                    double om  = (1.0 - hr) * (1.0 - std::exp(-hz1 * x0));
                    if (s11<=0||s1<=0||s01<=0||s0<=0) continue;

                    double rho0 = std::sqrt(v1/v);
                    double cd=-10,cf=10,ac=100; int it2=0;
                    while ((std::abs(ac-alpha)>alphaeps||cf-cd>ceps)&&it2<nbmaxiter) {
                        it2++; double cc=(cd+cf)/2;
                        ac=tte::bvn_integral(cc,c1v,rho0);
                        if (ac>alpha) cf=cc; else cd=cc;
                    }
                    if (cf-cd>ceps*10) continue;
                    double c2v=(cd+cf)/2;

                    double rho1=std::sqrt(s11/s1);
                    double cb1=std::sqrt(s01/s11)*(c1v-om1*std::sqrt(rate*DA1)/std::sqrt(s01));
                    double cb=std::sqrt(s0/s1)*(c2v-om*std::sqrt(rate*ta)/std::sqrt(s0));
                    if (tte::bvn_integral(cb,cb1,rho1)<1.0-beta) continue;

                    double pet=1.0-R::pnorm(c1v,0,1,1,0);
                    double en=ta*rate-std::max(0.0,ta-DA1)*rate*pet;
                    int n1_int=(int)std::ceil(DA1*rate);
                    RKJCand feasible = {di, p0d, n1_int, n, c1v, c2v, DA1, ta+x0, en, pet, hz0};
                    add_source_candidate(feasible);

                    // Track min-N
                    if (n < best_minN.n2 || (n == best_minN.n2 && en < best_minN.en)) {
                        best_minN = feasible;
                    }
                    // Track max-PET (highest prob of early stop)
                    if (pet > max_pet_seen || (pet == max_pet_seen && en < best_maxPET.en)) {
                        max_pet_seen = pet;
                        best_maxPET = feasible;
                    }
                }
            }
        }

        candidates.insert(candidates.end(), source_candidates.begin(), source_candidates.end());

        if (verbose) {
            int nc_before = (int)candidates.size();
            Rprintf("p0=%.3f (HR=%.3f): %d candidates so far",
                    p0d, hr, nc_before);
            if (r.valid) Rprintf(" [zoom: n2=%d, EN=%.1f]", r.n2, r.EN_H0);
            if (best_minN.n2 < 99999) Rprintf(" [minN: n2=%d]", best_minN.n2);
            if (max_pet_seen > 0) Rprintf(" [maxPET: %.2f]", max_pet_seen);
            Rprintf(" [feasible: %d]", (int)source_candidates.size());
            Rprintf("\n");
        }
        Rcpp::checkUserInterrupt();
    }

    // Count distinct full designs for reporting only. We do not deduplicate
    // before interval evaluation because the operating characteristics depend
    // on the design-time lambda0, so early collapsing can be unsafe.
    std::vector<RKJCand> distinct_report;
    for (auto& c : candidates) {
        bool dup = false;
        for (auto& u : distinct_report) {
            if (is_same_candidate(u, c)) { dup = true; break; }
        }
        if (!dup) distinct_report.push_back(c);
    }
    int nc_unique = (int)distinct_report.size();

    int nc = (int)candidates.size();
    if (nc == 0)
        return List::create(Named("error") = "No feasible designs found within n.ub.");

    if (verbose)
        Rprintf("\nGenerated %d candidates (%d distinct full designs).\n"
                "Evaluating on [%.3f, %.3f] (%d shared grid points, closed-form interval evaluation)...\n",
                nc, nc_unique, p0L, p0U, (int)grid.size());

    // ===== Step 2: Closed-form interval evaluation across the null interval =====
    int n_grid = (int)grid.size();

    std::vector<double> sup_a(nc), w_EN(nc), a_EN(nc);
    std::vector<std::vector<double>> EN_mat(nc, std::vector<double>(n_grid));

    for (int j = 0; j < nc; j++) {
        auto& c = candidates[j];
        double sa = 0, se = 0, me = -1e18;

        for (int i = 0; i < n_grid; i++) {
            double lt = -std::log(grid[i]) / x0;
            IntervalEvalResult ev = evaluate_candidate_interval(c, lt, x0, rate);
            if (!ev.valid) {
                return List::create(Named("error") =
                    "Failed to evaluate interval operating characteristics analytically.");
            }
            EN_mat[j][i] = ev.en;
            sa = std::max(sa, ev.alpha);
            se += ev.en;
            me = std::max(me, ev.en);
        }
        sup_a[j] = sa;
        w_EN[j]  = me;
        a_EN[j]  = se / n_grid;

        Rcpp::checkUserInterrupt();
    }

    // ===== Step 3: Keep only interval-robust candidates =====
    std::vector<int> robust_idx;
    for (int j = 0; j < nc; j++)
        if (sup_a[j] <= alpha) robust_idx.push_back(j);

    if (verbose)
        Rprintf("\nRetained robust candidates (sup_alpha <= %.3f): %d / %d\n",
                alpha, (int)robust_idx.size(), nc);

    if (robust_idx.empty())
        return List::create(Named("error") =
            "No BUDS-feasible designs. Try narrower [S0L,S0U] or larger n.ub.");

    // ===== Step 4: Build local-KJ benchmark at each evaluation grid point =====
    std::vector<LocalBenchmark> local_bench(n_grid);
    for (int i = 0; i < n_grid; i++) {
        double p0e = grid[i];
        double hz0e = -std::log(p0e) / x0;
        double hre = hz0e / hz1;
        LocalBenchmark lb = invalid_local_benchmark(p0e, hz0e);
        if (hre <= 1.0) {
            local_bench[i] = lb;
            continue;
        }
        if (cached_local_bench[i].valid) {
            local_bench[i] = cached_local_bench[i];
            continue;
        }
        int ns = tte_single_n(hz0e, hz1, hre, x0, alpha, beta);
        double tasingle = (double)ns / rate;
        TTEResult rloc = tte_zoom_search(
            hz0e, hz1, hre, x0, alpha, beta, rate,
            ns, tasingle, ceps, alphaeps, nbmaxiter
        );
        local_bench[i] = benchmark_from_zoom_result(p0e, hz0e, rloc);
    }

    for (int i = 0; i < n_grid; i++) {
        if (!local_bench[i].valid) {
            return List::create(Named("error") =
                "Failed to build local KJ benchmark on evaluation grid.");
        }
    }

    // Compute interval-level BUDS summaries for every generated candidate.
    std::vector<double> vwEN(nc), vaEN(nc), vwR_EN(nc), vaR_EN(nc), vwR_N(nc), vaR_N(nc);
    std::vector<int> vN(nc);

    for (int j = 0; j < nc; j++) {
        vN[j] = candidates[j].n2;
        double se = 0, me = -1e18;
        double sr_en = 0, mr_en = -1e18;
        double sr_n = 0, mr_n = -1e18;
        for (int i = 0; i < n_grid; i++) {
            double e  = EN_mat[j][i];
            double rg_en = e - local_bench[i].EN_local;
            double rg_n = (double)candidates[j].n2 - (double)local_bench[i].n2;

            se += e;
            sr_en += rg_en;
            sr_n += rg_n;
            me = std::max(me, e);
            mr_en = std::max(mr_en, rg_en);
            mr_n = std::max(mr_n, rg_n);
        }
        vwEN[j] = me; vaEN[j] = se / n_grid;
        vwR_EN[j] = mr_en;  vaR_EN[j] = sr_en / n_grid;
        vwR_N[j] = mr_n;    vaR_N[j] = sr_n / n_grid;
    }

    auto objective_score = [&](const std::string& obj, int idx) {
        if      (obj == "least_regret") return vwR_EN[idx];
        else if (obj == "least_en")     return vwEN[idx];
        else if (obj == "avg_en")       return vaEN[idx];
        else if (obj == "avg_regret")   return vaR_EN[idx];
        return (double)vN[idx];
    };

    auto choose_best = [&](const std::vector<int>& idx_set, const std::string& obj) {
        std::vector<int> idx = idx_set;
        double best = 1e18;
        for (int j : idx) best = std::min(best, objective_score(obj, j));
        std::vector<int> keep;
        for (int j : idx) if (std::abs(objective_score(obj, j) - best) < 1e-12) keep.push_back(j);
        idx.swap(keep);

        auto tiebreak = [&](const std::vector<double>& metric) {
            if (idx.size() <= 1) return;
            double b = 1e18;
            for (int j : idx) b = std::min(b, metric[j]);
            std::vector<int> t;
            for (int j : idx) if (std::abs(metric[j] - b) < 1e-12) t.push_back(j);
            idx.swap(t);
        };

        tiebreak(vwEN);
        if (idx.size() > 1) {
            std::vector<double> n_metric(nc);
            for (int j = 0; j < nc; j++) n_metric[j] = (double)vN[j];
            tiebreak(n_metric);
        }
        tiebreak(vaEN);
        return idx.front();
    };

    // ===== Step 5: Pool robust candidates across design points =====
    std::vector<int> pooled_idx;
    for (int idx : robust_idx) {
        bool dup = false;
        for (size_t pi = 0; pi < pooled_idx.size(); ++pi) {
            int kept = pooled_idx[pi];
            if (is_same_candidate(candidates[kept], candidates[idx])) {
                dup = true;
                if (vaEN[idx] < vaEN[kept]) pooled_idx[pi] = idx;
                break;
            }
        }
        if (!dup) pooled_idx.push_back(idx);
    }

    int nr = (int)pooled_idx.size();
    if (nr == 0)
        return List::create(Named("error") = "No pooled robust candidates remain.");

    if (verbose)
        Rprintf("Pooled robust candidate set: %d designs.\n", nr);

    // ===== Step 6: Select final winners by requested BUDS objectives =====
    std::vector<std::string> ol, oobj, obasis;
    std::vector<int> on1, on2;
    std::vector<double> oc1, oc2, oDA1, oDA2, opet, oen, osup, owEN, oaEN, op0d, owR, oaR;

    for (const auto& obj : objectives) {
        int ch = choose_best(pooled_idx, obj);
        auto& c = candidates[ch];
        std::string label = "BUDS";
        if (obj == "avg_en") label = "BUDS (Average EN)";
        else if (obj == "least_en") label = "BUDS (Least EN)";
        else if (obj == "least_regret") label = "BUDS (Least Regret)";
        else if (obj == "avg_regret")   label = "BUDS (Average Regret)";
        else if (obj == "min_n") label = "BUDS (Min N)";

        ol.push_back(label);
        oobj.push_back(obj);
        obasis.push_back(obj.find("regret") != std::string::npos ? std::string("en") : "none");
        on1.push_back(c.n1); on2.push_back(c.n2);
        oc1.push_back(c.c1); oc2.push_back(c.c2);
        oDA1.push_back(c.DA1); oDA2.push_back(c.DA2);
        opet.push_back(c.pet); oen.push_back(c.en);
        osup.push_back(sup_a[ch]);
        owEN.push_back(vwEN[ch]); oaEN.push_back(vaEN[ch]);
        owR.push_back(vwR_EN[ch]);
        oaR.push_back(vaR_EN[ch]);
        op0d.push_back(c.p0_d);
    }

    // Candidate pool table: retained robust candidates after pooling
    std::vector<int> cp_source, cp_n1, cp_n2;
    std::vector<double> cp_source_p0, cp_c1, cp_c2, cp_DA1, cp_DA2, cp_EN_design, cp_PET_design;
    std::vector<double> cp_sup, cp_wEN, cp_aEN, cp_wR_EN, cp_aR_EN, cp_wR_N, cp_aR_N;
    cp_source.reserve(nr);
    for (int k = 0; k < nr; k++) {
        int j = pooled_idx[k];
        auto& c = candidates[j];
        cp_source.push_back(c.source_idx + 1);
        cp_source_p0.push_back(c.p0_d);
        cp_n1.push_back(c.n1);
        cp_n2.push_back(c.n2);
        cp_c1.push_back(c.c1);
        cp_c2.push_back(c.c2);
        cp_DA1.push_back(c.DA1);
        cp_DA2.push_back(c.DA2);
        cp_EN_design.push_back(c.en);
        cp_PET_design.push_back(c.pet);
        cp_sup.push_back(sup_a[j]);
        cp_wEN.push_back(vwEN[j]);
        cp_aEN.push_back(vaEN[j]);
        cp_wR_EN.push_back(vwR_EN[j]);
        cp_aR_EN.push_back(vaR_EN[j]);
        cp_wR_N.push_back(vwR_N[j]);
        cp_aR_N.push_back(vaR_N[j]);
    }

    // Local benchmark table
    std::vector<double> lb_p0, lb_lambda, lb_c1, lb_c2, lb_DA1, lb_DA2, lb_EN;
    std::vector<int> lb_n1, lb_n2;
    std::vector<int> lb_valid;
    lb_p0.reserve(n_grid);
    for (int i = 0; i < n_grid; i++) {
        lb_p0.push_back(local_bench[i].p0_eval);
        lb_lambda.push_back(local_bench[i].lambda0_eval);
        lb_n1.push_back(local_bench[i].valid ? local_bench[i].n1 : NA_INTEGER);
        lb_n2.push_back(local_bench[i].valid ? local_bench[i].n2 : NA_INTEGER);
        lb_c1.push_back(local_bench[i].valid ? local_bench[i].c1 : NA_REAL);
        lb_c2.push_back(local_bench[i].valid ? local_bench[i].c2 : NA_REAL);
        lb_DA1.push_back(local_bench[i].valid ? local_bench[i].DA1 : NA_REAL);
        lb_DA2.push_back(local_bench[i].valid ? local_bench[i].DA2 : NA_REAL);
        lb_EN.push_back(local_bench[i].valid ? local_bench[i].EN_local : NA_REAL);
        lb_valid.push_back(local_bench[i].valid ? 1 : 0);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    NumericVector eg(n_grid);
    for (int i = 0; i < n_grid; i++) eg[i] = grid[i];

    return List::create(
        Named("designs") = DataFrame::create(
            Named("label") = ol,
            Named("objective") = oobj,
            Named("regret_basis") = obasis,
            Named("p0_designed") = op0d,
            Named("n1") = on1, Named("c1") = oc1,
            Named("n2") = on2, Named("c2") = oc2,
            Named("DA1") = oDA1, Named("DA2") = oDA2,
            Named("EN") = oen, Named("PET") = opet,
            Named("sup_alpha") = osup,
            Named("least_EN") = owEN,
            Named("avg_EN") = oaEN,
            Named("least_regret") = owR,
            Named("avg_regret") = oaR
        ),
        Named("candidate_pool") = DataFrame::create(
            Named("source_design_idx") = cp_source,
            Named("source_p0") = cp_source_p0,
            Named("n1") = cp_n1, Named("c1") = cp_c1,
            Named("n2") = cp_n2, Named("c2") = cp_c2,
            Named("DA1") = cp_DA1, Named("DA2") = cp_DA2,
            Named("EN_design") = cp_EN_design,
            Named("PET_design") = cp_PET_design,
            Named("sup_alpha") = cp_sup,
            Named("least_EN") = cp_wEN,
            Named("avg_EN") = cp_aEN,
            Named("least_regret_EN") = cp_wR_EN,
            Named("avg_regret_EN") = cp_aR_EN,
            Named("least_regret_N") = cp_wR_N,
            Named("avg_regret_N") = cp_aR_N
        ),
        Named("local_benchmarks") = DataFrame::create(
            Named("p0_eval") = lb_p0,
            Named("lambda0_eval") = lb_lambda,
            Named("n1_local") = lb_n1,
            Named("c1_local") = lb_c1,
            Named("n2_local") = lb_n2,
            Named("c2_local") = lb_c2,
            Named("DA1_local") = lb_DA1,
            Named("DA2_local") = lb_DA2,
            Named("EN_local") = lb_EN,
            Named("valid") = lb_valid
        ),
        Named("eval_grid") = eg,
        Named("n_total") = (int)candidates.size(),
        Named("n_unique") = nc_unique,
        Named("n_robust") = (int)robust_idx.size(),
        Named("n_pooled") = nr,
        Named("p0L") = p0L, Named("p0U") = p0U,
        Named("regret_basis") = regret_basis_s,
        Named("elapsed") = elapsed);
}
