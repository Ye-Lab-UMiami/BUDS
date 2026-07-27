// rcpp_binary.cpp -- Rcpp bridge to the C++ binary-endpoint design engine
//
// This file exposes binary_BUDS_cpp() to R via Rcpp. The R wrapper
// binary_BUDS() calls this function and wraps the result in the
// same S3 "binary_BUDS" class that the R interface produces.

#include <Rcpp.h>
#include "binary_core.h"

// [[Rcpp::export]]
Rcpp::List binary_BUDS_cpp(
    double alpha, double power, double p0, double p1,
    int n_lb, int n_ub,
    double p0L, double p0U,
    int pgrid_points,
    Rcpp::StringVector robust_objective,
    Rcpp::String robust_ref,
    double robust_max_inflation,
    int robust_N_cap,
    bool robust_auto_expand,
    double robust_expand_step,
    double robust_expand_max)
{
    // Convert R string vector to std::vector<std::string>
    std::vector<std::string> obj_vec;
    for (int i = 0; i < robust_objective.size(); ++i) {
        obj_vec.push_back(Rcpp::as<std::string>(robust_objective[i]));
    }

    // Call the C++ engine
    TwoStageDesigns result = binary_BUDS_engine(
        alpha, power, p0, p1,
        n_lb, n_ub, p0L, p0U, pgrid_points,
        obj_vec,
        std::string(robust_ref),
        robust_max_inflation,
        robust_N_cap,
        robust_auto_expand,
        robust_expand_step,
        robust_expand_max);

    // Build the design matrix (same layout as the pure-R version)
    // Row names: r1, n1, PET, r, EN, n, ratio, alpha1, alpha2, alpha, beta1, beta2, beta
    const int NROWS = 13;

    // Collect all columns: Optimal, Minimax, Balanced, robust...
    struct ColInfo { std::string name; DesignParams dp; };
    std::vector<ColInfo> cols;
    cols.push_back({"Optimal",  result.optimal});
    cols.push_back({"Minimax",  result.minimax});
    cols.push_back({"Balanced", result.balanced});
    for (auto& r : result.robust) {
        cols.push_back({r.label, r.params});
    }

    int ncols = (int)cols.size();
    Rcpp::NumericMatrix design_mat(NROWS, ncols);

    for (int j = 0; j < ncols; ++j) {
        auto& dp = cols[j].dp;
        if (!dp.valid()) {
            for (int i = 0; i < NROWS; ++i) design_mat(i, j) = NA_REAL;
        } else {
            design_mat(0,  j) = dp.r1;
            design_mat(1,  j) = dp.n1;
            design_mat(2,  j) = dp.pet;
            design_mat(3,  j) = dp.r;
            design_mat(4,  j) = dp.en;
            design_mat(5,  j) = dp.n;
            design_mat(6,  j) = dp.ratio;
            design_mat(7,  j) = dp.alpha1;
            design_mat(8,  j) = dp.alpha2;
            design_mat(9,  j) = dp.alpha;
            design_mat(10, j) = dp.beta1;
            design_mat(11, j) = dp.beta2;
            design_mat(12, j) = dp.beta;
        }
    }

    // Set row/col names
    Rcpp::CharacterVector rnames = Rcpp::CharacterVector::create(
        "r1", "n1", "PET", "r", "EN", "n", "ratio (n1:n2)",
        "alpha1", "alpha2", "alpha", "beta1", "beta2", "beta");
    Rcpp::CharacterVector cnames(ncols);
    for (int j = 0; j < ncols; ++j) cnames[j] = cols[j].name;

    design_mat.attr("dimnames") = Rcpp::List::create(rnames, cnames);

    // Return the same structure as the R version
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("designs") = design_mat,
        Rcpp::Named("inputs") = Rcpp::List::create(
            Rcpp::Named("alpha") = alpha,
            Rcpp::Named("power") = power,
            Rcpp::Named("beta")  = 1.0 - power,
            Rcpp::Named("p0")    = p0,
            Rcpp::Named("p1")    = p1,
            Rcpp::Named("p0L")   = p0L,
            Rcpp::Named("p0U")   = p0U,
            Rcpp::Named("n.lb")  = n_lb,
            Rcpp::Named("n.ub")  = n_ub,
            Rcpp::Named("pgrid_points")        = pgrid_points,
            Rcpp::Named("robust_objective")    = robust_objective,
            Rcpp::Named("robust_ref")          = robust_ref,
            Rcpp::Named("robust_max_inflation") = robust_max_inflation,
            Rcpp::Named("robust_N_cap")        = robust_N_cap
        ),
        Rcpp::Named("meta") = Rcpp::List::create(
            Rcpp::Named("execution_time") = Rcpp::String(
                std::to_string(result.elapsed_sec).substr(0,4) + " seconds"),
            Rcpp::Named("engine") = "cpp"
        )
    );
    out.attr("class") = "binary_BUDS";
    return out;
}
