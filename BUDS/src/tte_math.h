#ifndef TTE_MATH_H
#define TTE_MATH_H

#include <cmath>
#include <Rcpp.h>

namespace tte {

// 20-point Gauss-Legendre quadrature
static const int NQUAD = 20;
static const double GL_NODES[NQUAD] = {
    -0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
    -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
    -0.5108670019508271, -0.3737060887154195, -0.2277858511416451,
    -0.0765265211334973,  0.0765265211334973,  0.2277858511416451,
     0.3737060887154195,  0.5108670019508271,  0.6360536807265150,
     0.7463319064601508,  0.8391169718222188,  0.9122344282513259,
     0.9639719272779138,  0.9931285991850949
};
static const double GL_WEIGHTS[NQUAD] = {
     0.0176140071391521,  0.0406014298003869,  0.0626720483341091,
     0.0832767415767048,  0.1019301198172404,  0.1181945319615184,
     0.1316886384491766,  0.1420961093183820,  0.1491729864726037,
     0.1527533871307258,  0.1527533871307258,  0.1491729864726037,
     0.1420961093183820,  0.1316886384491766,  0.1181945319615184,
     0.1019301198172404,  0.0832767415767048,  0.0626720483341091,
     0.0406014298003869,  0.0176140071391521
};

// Bivariate normal tail probability used by the two-stage TTE calculations.
// This evaluates
//   integral_{-inf}^{upper} phi(z) * Phi((c1 - rho * z) / sqrt(1-rho^2)) dz
// which is the joint probability needed for:
// - type-I error calibration when solving for c2
// - power evaluation under the working alternative
//
// Implementation details:
// - the left tail is truncated to a finite interval because the standard
//   normal density is negligible far into the tail
// - the truncated interval is mapped to [-1, 1]
// - Gauss-Legendre quadrature then approximates the integral efficiently
inline double bvn_integral(double upper, double c1, double rho) {
    // Truncate the far left tail for numerical integration.
    double lower = std::max(upper - 8.0, -8.0);
    // Affine map from [-1, 1] to [lower, upper].
    double half  = (upper - lower) / 2.0;
    double mid   = (upper + lower) / 2.0;
    // sqrt(1-rho^2) appears in the conditional normal representation.
    double sr    = std::sqrt(1.0 - rho * rho);
    double result = 0.0;
    for (int i = 0; i < NQUAD; i++) {
        double z = mid + half * GL_NODES[i];
        // Evaluate the integrand at the mapped quadrature node z.
        result += GL_WEIGHTS[i] * R::dnorm(z,0,1,0) * R::pnorm((c1-rho*z)/sr,0,1,1,0);
    }
    // Multiply by the Jacobian term from the interval transformation.
    return result * half;
}

// Variance / information component from Appendix A of Belin et al. (2017).
// This helper is reused in the H0/H1 information terms inside the TTE design
// calculations. The input BORNE is the relevant follow-up bound:
// typically min(DA1, x0) at stage 1.
inline double fct(double DA1, double ta, double hz, double BORNE) {
    double denom = ta * hz;
    double e     = std::exp(-hz * BORNE);
    return (DA1*(1.0-e) + BORNE*e + e/hz - 1.0/hz) / denom;
}

// -------------------------------------------------------------------
// Simulate rejection probability for a FIXED design under a given
// true hazard rate (lambda_true).
// The design was built with lambda0_design as the null reference.
// Returns: (reject_rate, PET, EN)
// -------------------------------------------------------------------
struct SimResult {
    double reject_rate;
    double pet;
    double en;
};

inline SimResult simulate_design(
    int n1, double c1, int n2, double c2,
    double DA1, double DA2, double x0, double rate,
    double lambda0_design,  // lambda used in the test statistic (null reference)
    double lambda_true,     // true hazard rate for data generation
    int nsim, unsigned int seed)
{
    // ta: total accrual duration implied by the maximum sample size n2.
    double ta = (double)n2 / rate;
    int n_reject = 0, n_early = 0;
    double total_n = 0.0;

    // Simple fast RNG for the inner simulation loop.
    // We avoid the R RNG here because this function is called very often.
    auto next_uniform = [&seed]() -> double {
        seed = seed * 1664525u + 1013904223u;
        return (double)(seed >> 1) / 2147483648.0;
    };
    // Generate an exponential event time by inverse transform.
    auto next_exp = [&](double rate_param) -> double {
        double u = next_uniform();
        if (u < 1e-15) u = 1e-15;
        return -std::log(u) / rate_param;
    };

    for (int sim = 0; sim < nsim; sim++) {
        // Generate accrual times uniformly over the accrual window [0, ta].
        // Event times follow an exponential model with hazard lambda_true.
        // The code generates accruals unsorted, then sorts them.
        std::vector<double> A(n2), Ti(n2);
        for (int j = 0; j < n2; j++) {
            A[j]  = next_uniform() * ta;
            Ti[j] = next_exp(lambda_true);
        }
        std::sort(A.begin(), A.end());

        // Stage 1 log-rank style statistic under the design null lambda0_design.
        // O1: observed events by stage 1
        // E1: expected events under lambda0_design
        double O1 = 0, E1 = 0;
        for (int j = 0; j < n1; j++) {
            double C_j = std::min(DA1 - A[j], x0);
            double X_j = std::min(Ti[j], C_j);
            if (Ti[j] <= C_j) O1 += 1.0;
            E1 += lambda0_design * X_j;
        }
        double Z1 = (E1 > 0) ? (O1 - E1) / std::sqrt(E1) : 0.0;

        // If stage 1 crosses c1, the design stops early for futility.
        if (Z1 > c1) {
            n_early++;
            total_n += n1;
            continue;
        }

        // Final analysis using all n2 patients and the stage-2 analysis time DA2.
        double O2 = 0, E2 = 0;
        for (int j = 0; j < n2; j++) {
            double C_j = std::min(DA2 - A[j], x0);
            double X_j = std::min(Ti[j], C_j);
            if (Ti[j] <= C_j) O2 += 1.0;
            E2 += lambda0_design * X_j;
        }
        double Z2 = (E2 > 0) ? (O2 - E2) / std::sqrt(E2) : 0.0;
        total_n += n2;
        // In this codebase, rejection occurs for sufficiently small Z2.
        if (Z2 <= c2) n_reject++;
    }

    // Return Monte Carlo estimates of rejection rate, PET, and EN.
    SimResult res;
    res.reject_rate = (double)n_reject / nsim;
    res.pet = (double)n_early / nsim;
    res.en  = total_n / nsim;
    return res;
}

} // namespace tte
#endif
