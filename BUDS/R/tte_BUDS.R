#' Convert Hazard Rate to Survival Probability
#'
#' Internal helper for resolving TTE inputs supplied on the hazard-rate scale.
#'
#' @param lambda Hazard rate, or \code{NULL}.
#' @param x0 Restricted follow-up time.
#' @param label Input label used in error messages.
#'
#' @return Survival probability \eqn{\exp(-\lambda x0)}, or \code{NULL}.
#' @keywords internal
.tte_survival_from_lambda <- function(lambda, x0, label) {
  if (is.null(lambda)) return(NULL)
  if (!is.finite(x0) || x0 <= 0) stop("x0 must be positive when hazard rates are supplied.")
  if (!is.finite(lambda) || lambda <= 0) stop(label, " must be positive.")
  exp(-lambda * x0)
}

#' Resolve Time-to-Event Inputs to Survival Scale
#'
#' Internal helper that reconciles survival-probability, hazard-rate, and hazard
#' ratio inputs before running the r-KJ search.
#'
#' @param S0,S0L,S0U Null survival inputs at \code{x0}.
#' @param S1 Alternative survival probability at \code{x0}.
#' @param HR Hazard ratio \eqn{\lambda_0 / \lambda_1}.
#' @param lambda1 Alternative hazard rate.
#' @param lambda0,lambda0L,lambda0U Null hazard inputs.
#' @param require_interval Logical; require a proper interval.
#' @param x0 Restricted follow-up time.
#'
#' @return A list of resolved survival, hazard-ratio, and hazard-rate inputs.
#' @keywords internal
.tte_resolve_survival_inputs <- function(S0 = NULL, S0L = NULL, S0U = NULL, S1 = NULL,
                                         HR = NULL, lambda1 = NULL,
                                         lambda0 = NULL, lambda0L = NULL, lambda0U = NULL,
                                         require_interval = FALSE, x0 = NULL) {
  S0_from_lambda <- .tte_survival_from_lambda(lambda0, x0, "lambda0")
  S0L_from_lambda <- .tte_survival_from_lambda(lambda0U, x0, "lambda0U")
  S0U_from_lambda <- .tte_survival_from_lambda(lambda0L, x0, "lambda0L")

  if (!is.null(S0) && !is.null(S0_from_lambda) &&
      !isTRUE(all.equal(S0, S0_from_lambda, tolerance = 1e-8))) {
    stop("S0 and lambda0 are not numerically consistent.")
  }
  if (!is.null(S0L) && !is.null(S0L_from_lambda) &&
      !isTRUE(all.equal(S0L, S0L_from_lambda, tolerance = 1e-8))) {
    stop("S0L and lambda0L are not numerically consistent.")
  }
  if (!is.null(S0U) && !is.null(S0U_from_lambda) &&
      !isTRUE(all.equal(S0U, S0U_from_lambda, tolerance = 1e-8))) {
    stop("S0U and lambda0U are not numerically consistent.")
  }

  if (is.null(S0)) S0 <- S0_from_lambda
  if (is.null(S0L)) S0L <- S0L_from_lambda
  if (is.null(S0U)) S0U <- S0U_from_lambda

  if (is.null(S0) && !is.null(lambda0L) && !is.null(lambda0U)) {
    S0 <- exp(-((lambda0L + lambda0U) / 2) * x0)
  }
  if (is.null(S0) && !is.null(S0L) && !is.null(S0U)) {
    S0 <- (S0L + S0U) / 2
  }
  if (is.null(S0)) {
    stop("Supply S0 or lambda0, or supply both S0L/S0U (or lambda0L/lambda0U) so S0 can be derived as the midpoint.")
  }
  if (!is.finite(S0) || S0 <= 0 || S0 >= 1) {
    stop("S0 must be in (0,1).")
  }

  has_interval <- !(is.null(S0L) && is.null(S0U))
  if (!has_interval) {
    S0L <- S0
    S0U <- S0
  } else {
    if (is.null(S0L) || is.null(S0U)) {
      stop("Supply both S0L and S0U, both lambda0L and lambda0U, or leave both NULL.")
    }
    if (!is.finite(S0L) || S0L <= 0 || S0L >= 1) stop("S0L must be in (0,1).")
    if (!is.finite(S0U) || S0U <= 0 || S0U >= 1) stop("S0U must be in (0,1).")
    if (S0L > S0U) stop("S0L must be <= S0U.")
  }

  S1_from_lambda1 <- .tte_survival_from_lambda(lambda1, x0, "lambda1")

  if (!is.null(S1) && !is.null(S1_from_lambda1) &&
      !isTRUE(all.equal(S1, S1_from_lambda1, tolerance = 1e-8))) {
    stop("S1 and lambda1 are not numerically consistent.")
  }
  if (is.null(S1)) S1 <- S1_from_lambda1

  if (!is.null(HR) && !is.null(S1)) {
    S1_from_hr <- exp(log(S0) / HR)
    if (!isTRUE(all.equal(S1, S1_from_hr, tolerance = 1e-8))) {
      stop("Provide either S1 or HR, or provide numerically consistent values.")
    }
  }
  if (is.null(HR) && is.null(S1)) {
    stop("Supply either S1 or HR.")
  }
  if (!is.null(HR)) {
    if (!is.finite(HR) || HR <= 1) stop("HR must be > 1.")
    S1 <- exp(log(S0) / HR)
  } else {
    if (!is.finite(S1) || S1 <= 0 || S1 >= 1) stop("S1 must be in (0,1).")
    if (S1 <= S0) stop("S1 must be > S0.")
    HR <- log(S0) / log(S1)
  }

  if (require_interval && identical(S0L, S0U)) {
    stop("A proper interval requires S0L < S0U.")
  }
  if (S1 <= S0U) {
    stop("S1 must be > S0U.")
  }

  lambda1_resolved <- .tte_lambda_from_survival(S1, x0, "S1")

  list(S0 = S0, S0L = S0L, S0U = S0U, S1 = S1, HR = HR,
       lambda1 = lambda1_resolved)
}

#' Convert Survival Probability to Hazard Rate
#'
#' Internal helper for resolving TTE inputs supplied on the survival scale.
#'
#' @param survival Survival probability at \code{x0}, or \code{NULL}.
#' @param x0 Restricted follow-up time.
#' @param label Input label used in error messages.
#'
#' @return Hazard rate \eqn{-\log(S) / x0}, or \code{NULL}.
#' @keywords internal
.tte_lambda_from_survival <- function(survival, x0, label) {
  if (is.null(survival)) return(NULL)
  if (!is.finite(x0) || x0 <= 0) stop("x0 must be positive.")
  if (!is.finite(survival) || survival <= 0 || survival >= 1) {
    stop(label, " must be in (0,1).")
  }
  -log(survival) / x0
}

#' Compute a Classic Restricted Kwak-Jung Design
#'
#' Internal wrapper around the C++ r-KJ search used to construct the classic
#' comparator design at a single planning null value.
#'
#' @param S0,S1 Null and alternative survival probabilities at \code{x0}.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#' @param alpha,beta Type I and type II error constraints.
#' @param ceps,alphaeps,nbmaxiter Numerical search controls.
#' @param verbose Logical; print search diagnostics.
#'
#' @return A list containing planning parameters, single-stage quantities,
#'   two-stage design parameters, and iteration metadata.
#' @keywords internal
.tte_classic_design <- function(S0, S1, x0, rate,
                                alpha = 0.05, beta = 0.10,
                                ceps = 0.001, alphaeps = 0.001,
                                nbmaxiter = 100L, verbose = FALSE) {
  if (!is.finite(S0) || S0 <= 0 || S0 >= 1) stop("S0 must be in (0,1).")
  if (!is.finite(S1) || S1 <= 0 || S1 >= 1) stop("S1 must be in (0,1).")
  if (S1 <= S0) stop("S1 must be > S0.")
  if (x0 <= 0) stop("x0 must be positive.")
  if (rate <= 0) stop("rate must be positive.")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")
  if (beta <= 0 || beta >= 1) stop("beta must be between 0 and 1.")

  hr <- log(S0) / log(S1)
  hz0 <- .tte_lambda_from_survival(S0, x0, "S0")
  hz1 <- hz0 / hr

  if (verbose) {
    message(sprintf("lambda0 = %.4f, HR = %.3f", hz0, hr))
    message(sprintf("Median H0 = %.2f, Median H1 = %.2f",
                    log(2) / hz0, log(2) / hz1))
  }

  za <- qnorm(1 - alpha)
  zb <- qnorm(1 - beta)
  s0 <- sqrt(hr * (1 - exp(-hz1 * x0)))
  hzb <- (hz0 + hz1) / 2
  s1 <- sqrt(1 - exp(-hzb * x0))
  om <- (1 - exp(-hz1 * x0)) * (1 - hr)
  nsingle <- ceiling((s0 * za + s1 * zb)^2 / om^2)
  tasingle <- nsingle / rate

  search <- tte_grid_search_cpp(
    hz0 = hz0, hz1 = hz1, hr = hr, x0 = x0,
    alpha = alpha, beta = beta, rate = rate,
    nsingle = nsingle, tasingle = tasingle,
    ceps = ceps, alphaeps = alphaeps,
    nbmaxiter = as.integer(nbmaxiter), verbose = verbose
  )

  list(
    param = data.frame(
      S0 = S0,
      S1 = S1,
      HR = hr,
      alpha = alpha,
      beta = beta,
      rate = rate,
      x0 = x0,
      hz0 = hz0,
      hz1 = hz1
    ),
    single_stage = data.frame(
      n = nsingle,
      ta = tasingle,
      c = qnorm(1 - alpha)
    ),
    two_stage = data.frame(
      n1 = search$n1,
      c1 = search$c1,
      n2 = search$n2,
      c2 = search$c2,
      DA1 = search$DA1,
      DA2 = search$DA2,
      EN_H0 = search$EN_H0,
      PET = search$PET
    ),
    iters = search$iters
  )
}

#' Convert r-KJ Design Data Frame to Matrix
#'
#' Internal helper that stores selected TTE designs in the matrix format used by
#' \code{"tte_BUDS"} objects.
#'
#' @param design_df Data frame of selected design rows.
#'
#' @return Numeric matrix with one column per selected design.
#' @keywords internal
.tte_design_df_to_matrix <- function(design_df) {
  row_order <- c(
    "S0_designed", "c1", "n1", "DA1", "PET",
    "c2", "n2", "DA2", "EN",
    "sup_alpha"
  )
  if (!nrow(design_df)) {
      out <- matrix(numeric(0), nrow = length(row_order), ncol = 0)
    rownames(out) <- row_order
    return(out)
  }

  mat <- t(as.matrix(design_df[, row_order, drop = FALSE]))
  storage.mode(mat) <- "double"
  colnames(mat) <- design_df$design
  rownames(mat) <- row_order
  mat
}

#' Convert r-KJ Design Matrix to Data Frame
#'
#' Internal helper that converts the matrix stored in \code{"tte_BUDS"} objects
#' back to a row-wise data frame.
#'
#' @param design_matrix Numeric design matrix.
#'
#' @return Data frame with one row per selected design.
#' @keywords internal
.tte_design_matrix_to_df <- function(design_matrix) {
  if (is.null(dim(design_matrix)) || !length(design_matrix)) {
      return(data.frame(
        design = character(0),
        S0_designed = numeric(0),
        c1 = numeric(0),
      n1 = integer(0),
      DA1 = numeric(0),
      PET = numeric(0),
      c2 = numeric(0),
      n2 = integer(0),
      DA2 = numeric(0),
      EN = numeric(0),
      sup_alpha = numeric(0)
    ))
  }

  df <- as.data.frame(t(design_matrix), check.names = FALSE)
  df$design <- rownames(df)
  rownames(df) <- NULL
  df$n1 <- as.integer(round(df$n1))
  df$n2 <- as.integer(round(df$n2))
  df <- df[, c(
    "design", "S0_designed", "n1", "c1", "DA1", "PET",
    "n2", "c2", "DA2", "EN", "sup_alpha"
  )]
  names(df)[names(df) == "DA1"] <- "DA1"
  names(df)[names(df) == "DA2"] <- "DA2"
  names(df)[names(df) == "c1"] <- "c1"
  names(df)[names(df) == "c2"] <- "c2"
  df
}

#' Build the Classic r-KJ Summary Row
#'
#' Internal helper that evaluates the classic r-KJ comparator over the requested
#' benchmark range so it can be reported alongside BUDS designs.
#'
#' @param classic_design List returned by \code{.tte_classic_design()}.
#' @param S0,S0L,S0U Survival values under null hypothesis.
#' @param S1 Alternative survival probability.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#' @param grid_points Number of grid points for null hypothesis evaluation under BUDS.
#'
#' @return One-row data frame containing the classic comparator summary.
#' @keywords internal
.tte_build_classic_summary_row <- function(classic_design, S0, S0L, S0U, S1, x0, rate,
                                           grid_points = 101L) {
  ts <- classic_design$two_stage[1, , drop = FALSE]
  grid <- seq(S0L, S0U, length.out = max(1L, as.integer(grid_points)))

  alpha_grid <- vapply(grid, function(pp) {
    pr_reject_tte_BUDS(
      n1 = as.integer(ts$n1), c1 = ts$c1,
      n2 = as.integer(ts$n2), c2 = ts$c2,
      DA1 = ts$DA1, DA2 = ts$DA2,
      S0_design = S0, s = pp, x0 = x0, rate = rate
    )
  }, numeric(1))
  pet_grid <- vapply(grid, function(pp) {
    pet_tte_BUDS(
      n1 = as.integer(ts$n1), c1 = ts$c1,
      n2 = as.integer(ts$n2), c2 = ts$c2,
      DA1 = ts$DA1, DA2 = ts$DA2,
      S0_design = S0, s = pp, x0 = x0, rate = rate
    )
  }, numeric(1))
  en_grid <- vapply(grid, function(pp) {
    en_tte_BUDS(
      n1 = as.integer(ts$n1), c1 = ts$c1,
      n2 = as.integer(ts$n2), c2 = ts$c2,
      DA1 = ts$DA1, DA2 = ts$DA2,
      S0_design = S0, s = pp, x0 = x0, rate = rate
    )
  }, numeric(1))

    data.frame(
      design = "r-KJ",
      S0_designed = S0,
      c1 = ts$c1,
    n1 = as.integer(ts$n1),
    DA1 = ts$DA1,
    PET = pet_tte_BUDS(
      n1 = as.integer(ts$n1), c1 = ts$c1,
      n2 = as.integer(ts$n2), c2 = ts$c2,
      DA1 = ts$DA1, DA2 = ts$DA2,
      S0_design = S0, s = S0, x0 = x0, rate = rate
    ),
    c2 = ts$c2,
    n2 = as.integer(ts$n2),
    DA2 = ts$DA2,
    EN = en_tte_BUDS(
      n1 = as.integer(ts$n1), c1 = ts$c1,
      n2 = as.integer(ts$n2), c2 = ts$c2,
      DA1 = ts$DA1, DA2 = ts$DA2,
      S0_design = S0, s = S0, x0 = x0, rate = rate
    ),
    sup_alpha = max(alpha_grid, na.rm = TRUE),
    least_EN = max(en_grid, na.rm = TRUE),
    avg_EN = mean(en_grid, na.rm = TRUE),
    least_regret = 0,
    avg_regret = 0
  )
}

#'  Robust Restricted Kwak-Jung Design under BUDS
#'
#' Searches for r-KJ designs robust across a benchmark range for survival at
#' \code{x0}. The returned object mirrors the top-level structure of
#' \code{\link{binary_BUDS}}, with \code{designs}, \code{inputs}, and
#' \code{meta} components.
#'
#' The public interface uses survival notation \code{S0}, \code{S0L},
#' \code{S0U}, and either \code{S1} or \code{HR}.
#'
#' @param S0 Planning survival probability at \code{x0} used to construct the
#'   classic comparator \code{"r-KJ"} design.
#' @param S0L,S0U Null survival interval at \code{x0}. If both are omitted, the
#'   function returns only the classic r-KJ design evaluated at \code{S0}.
#' @param lambda0 Optional planning hazard rate at \code{x0}; converted
#'   internally to \code{S0 = exp(-lambda0 * x0)}.
#' @param lambda0L,lambda0U Optional null hazard-rate benchmark range at \code{x0}. If
#'   supplied, these are converted internally to \code{S0L/S0U}; when
#'   \code{S0} is omitted, the planning value is derived from the midpoint
#'   hazard rate.
#' @param S1 Target survival at \code{x0}. Supply either \code{S1} or
#'   \code{HR}.
#' @param x0 Follow-up window.
#' @param rate Accrual rate.
#' @param HR Hazard ratio \eqn{\lambda_0 / \lambda_1}. If provided, \code{S1}
#'   is derived internally as \eqn{\exp\{\log(S0) / HR\}}.
#' @param alpha One-sided type I error. Default 0.05.
#' @param beta Type II error. Default 0.10.
#' @param n.ub Maximum total sample size considered.
#' @param n_design_points Number of shared grid points used for candidate
#'   generation and interval evaluation.
#' @param n_eval_grid Legacy compatibility argument; ignored in favor of
#'   \code{n_design_points}.
#' @param nsim_eval Legacy compatibility argument retained for API stability.
#' @param BUDS_objective BUDS objective(s) used to select winner(s). Publicly
#'   supported defaults are \code{"avg_en"} and \code{"least_regret"}.
#' @param verbose Print progress. Default FALSE.
#' @param ... Additional arguments passed to the internal local r-KJ search and the
#'   C++ search engine (\code{ceps}, \code{alphaeps}, \code{nbmaxiter}).
#'
#' @return An object of class \code{"tte_BUDS"} with:
#' \describe{
#'   \item{designs}{Numeric matrix of selected designs, one column per design.}
#'   \item{inputs}{List of user inputs and resolved defaults.}
#'   \item{meta}{Execution metadata.}
#' }
#'
#' @examples
#' \dontrun{
#' d <- tte_BUDS(
#'   S0 = 0.49, S0L = 0.48, S0U = 0.52, S1 = 0.707,
#'   x0 = 1, rate = 15, n.ub = 100
#' )
#' d
#' }
#'
#' @export
tte_BUDS <- function(S0L = NULL, S0U = NULL, S1 = NULL, x0, rate, S0 = NULL,
                    lambda0 = NULL, lambda0L = NULL, lambda0U = NULL,
                    HR = NULL, lambda1 = NULL,
                    alpha = 0.05, beta = 0.10,
                    n.ub = 150L,
                    n_design_points = 21L, n_eval_grid = 21L,
                    nsim_eval = 5000L,
                    BUDS_objective = c("avg_en", "least_regret"),
                    verbose = FALSE, ...) {
  input_scale_null <- if (any(!vapply(list(S0, S0L, S0U), is.null, logical(1)))) {
    "survival"
  } else {
    "hazard"
  }
  input_scale_alt <- if (!is.null(lambda1)) {
    "lambda1"
  } else if (!is.null(S1)) {
    "S1"
  } else {
    "HR"
  }

  if (x0 <= 0) stop("x0 must be positive.")
  if (rate <= 0) stop("rate must be positive.")

  surv <- .tte_resolve_survival_inputs(
    S0 = S0, S0L = S0L, S0U = S0U, S1 = S1, HR = HR, lambda1 = lambda1,
    lambda0 = lambda0, lambda0L = lambda0L, lambda0U = lambda0U,
    x0 = x0
  )
  S0 <- surv$S0
  S0L <- surv$S0L
  S0U <- surv$S0U
  S1 <- surv$S1
  HR <- surv$HR
  lambda1_resolved <- surv$lambda1
  if (n.ub < 2) stop("n.ub must be >= 2.")
  if (length(BUDS_objective) == 0) stop("BUDS_objective must contain at least one objective.")

  valid_obj <- c("avg_en", "least_regret", "avg_regret", "least_en", "min_n")
  BUDS_objective <- unique(tolower(BUDS_objective))
  BUDS_objective <- BUDS_objective[BUDS_objective %in% valid_obj]
  if (length(BUDS_objective) == 0) {
    stop("BUDS_objective must contain at least one of: ", paste(valid_obj, collapse = ", "))
  }

  dots <- list(...)
  ceps <- if (!is.null(dots$ceps)) dots$ceps else 0.001
  alphaeps <- if (!is.null(dots$alphaeps)) dots$alphaeps else 0.001
  nbmaxiter <- if (!is.null(dots$nbmaxiter)) dots$nbmaxiter else 100L
  lambda0_resolved <- .tte_lambda_from_survival(S0, x0, "S0")
  lambda_interval_resolved <- sort(c(
    .tte_lambda_from_survival(S0L, x0, "S0L"),
    .tte_lambda_from_survival(S0U, x0, "S0U")
  ))
  lambda0L_resolved <- lambda_interval_resolved[1]
  lambda0U_resolved <- lambda_interval_resolved[2]

  if (!identical(as.integer(n_eval_grid), as.integer(n_design_points)) && isTRUE(verbose)) {
    message("Using a single shared grid: n_eval_grid is ignored and n_design_points is used for both generation and evaluation.")
  }
  n_eval_grid <- as.integer(n_design_points)

  classic_design <- .tte_classic_design(
    S0 = S0, S1 = S1, x0 = x0, rate = rate,
    alpha = alpha, beta = beta,
    ceps = ceps, alphaeps = alphaeps,
    nbmaxiter = as.integer(nbmaxiter),
    verbose = FALSE
  )

  is_single_point <- identical(S0L, S0U)
  if (is_single_point) {
    classic_row <- .tte_build_classic_summary_row(
      classic_design = classic_design,
      S0 = S0, S0L = S0, S0U = S0, S1 = S1,
      x0 = x0, rate = rate, grid_points = 1L
    )
    objective_map <- c(
      "avg_en" = "BUDS (Average EN)",
      "least_regret" = "BUDS (Least Regret)"
    )
    objective_keep <- BUDS_objective[BUDS_objective %in% names(objective_map)]
    if (length(objective_keep) > 0) {
      BUDS_rows <- do.call(rbind, lapply(objective_keep, function(obj) {
        row <- classic_row
        row$design <- unname(objective_map[obj])
        row
      }))
      classic_row <- rbind(classic_row, BUDS_rows)
    }
    out <- list(
      designs = .tte_design_df_to_matrix(classic_row),
      inputs = list(
        alpha = alpha,
        power = 1 - beta,
        beta = beta,
        S0 = S0,
        S1 = S1,
        S0L = S0,
        S0U = S0,
        lambda0 = lambda0,
        lambda0L = lambda0L,
        lambda0U = lambda0U,
        lambda1 = lambda1_resolved,
        lambda0_resolved = lambda0_resolved,
        lambda0L_resolved = lambda0L_resolved,
        lambda0U_resolved = lambda0U_resolved,
        HR = HR,
        x0 = x0,
        rate = rate,
        n.ub = as.integer(n.ub),
        grid_points = 1L,
        robust_objective = objective_keep,
        input_scale_null = input_scale_null,
        input_scale_alt = input_scale_alt
      ),
      meta = list(
        execution_time = "0.0 seconds",
        engine = "cpp",
        n_total = 1L,
        n_unique = 1L,
        n_robust = 1L,
        n_pooled = 1L
      )
    )
    class(out) <- c("tte_BUDS", "BUDS_design")
    return(out)
  }

  result <- tte_BUDS_search_cpp(
    S0L, S0U, S1, x0, rate,
    alpha, beta,
    as.integer(n_design_points),
    as.integer(n_eval_grid),
    as.integer(nsim_eval),
    as.integer(n.ub),
    ceps, alphaeps,
    as.integer(nbmaxiter),
    BUDS_objective,
    "EN",
    verbose
  )
  if (!is.null(result$error)) stop(result$error)

  if ("label" %in% names(result$designs) && !"design" %in% names(result$designs)) {
    names(result$designs)[names(result$designs) == "label"] <- "design"
  }
  if ("p0_designed" %in% names(result$designs) && !"S0_designed" %in% names(result$designs)) {
    names(result$designs)[names(result$designs) == "p0_designed"] <- "S0_designed"
  }
  result$designs$objective <- NULL
  result$designs$regret_basis <- NULL

  if (nrow(result$designs) > 0) {
    result$designs$PET <- vapply(seq_len(nrow(result$designs)), function(i) {
      row <- result$designs[i, , drop = FALSE]
      pet_tte_BUDS(
        n1 = as.integer(row$n1), c1 = row$c1,
        n2 = as.integer(row$n2), c2 = row$c2,
        DA1 = row$DA1, DA2 = row$DA2,
        S0_design = row$S0_designed, s = S0U, x0 = x0, rate = rate
      )
    }, numeric(1))
    result$designs$EN <- vapply(seq_len(nrow(result$designs)), function(i) {
      row <- result$designs[i, , drop = FALSE]
      en_tte_BUDS(
        n1 = as.integer(row$n1), c1 = row$c1,
        n2 = as.integer(row$n2), c2 = row$c2,
        DA1 = row$DA1, DA2 = row$DA2,
        S0_design = row$S0_designed, s = S0U, x0 = x0, rate = rate
      )
    }, numeric(1))
  }

  classic_row <- .tte_build_classic_summary_row(
    classic_design = classic_design,
    S0 = S0, S0L = S0L, S0U = S0U, S1 = S1,
    x0 = x0, rate = rate, grid_points = n_eval_grid
  )
  bench <- result$local_benchmarks
  bench_valid <- as.logical(bench$valid)
  grid <- seq(S0L, S0U, length.out = n_eval_grid)
  classic_en_grid <- vapply(grid, function(pp) {
    en_tte_BUDS(
      n1 = as.integer(classic_design$two_stage$n1),
      c1 = classic_design$two_stage$c1,
      n2 = as.integer(classic_design$two_stage$n2),
      c2 = classic_design$two_stage$c2,
      DA1 = classic_design$two_stage$DA1,
      DA2 = classic_design$two_stage$DA2,
      S0_design = S0,
      s = pp,
      x0 = x0,
      rate = rate
    )
  }, numeric(1))
  classic_rg_grid <- classic_en_grid - bench$EN_local
  classic_rg_grid[!bench_valid] <- NA_real_
  classic_row$least_regret <- max(classic_rg_grid, na.rm = TRUE)
  classic_row$avg_regret <- mean(classic_rg_grid, na.rm = TRUE)

  selected_df <- rbind(
    classic_row,
    result$designs[, names(classic_row), drop = FALSE]
  )
  design_matrix <- .tte_design_df_to_matrix(selected_df)

  out <- list(
    designs = design_matrix,
    inputs = list(
      alpha = alpha,
      power = 1 - beta,
      beta = beta,
      S0 = S0,
      S1 = S1,
      S0L = S0L,
      S0U = S0U,
      lambda0 = lambda0,
      lambda0L = lambda0L,
      lambda0U = lambda0U,
      lambda1 = lambda1_resolved,
      lambda0_resolved = lambda0_resolved,
      lambda0L_resolved = lambda0L_resolved,
      lambda0U_resolved = lambda0U_resolved,
      HR = HR,
      x0 = x0,
      rate = rate,
      n.ub = as.integer(n.ub),
      grid_points = as.integer(n_design_points),
      robust_objective = BUDS_objective,
      input_scale_null = input_scale_null,
      input_scale_alt = input_scale_alt
    ),
    meta = list(
      execution_time = sprintf("%.1f seconds", result$elapsed),
      engine = "cpp",
      n_total = result$n_total,
      n_unique = result$n_unique,
      n_robust = result$n_robust,
      n_pooled = result$n_pooled
    )
  )
  class(out) <- c("tte_BUDS", "BUDS_design")
  out
}

#' Print an r-KJ Design Object under BUDS
#'
#' Displays the input settings and selected TTE designs stored in an
#' \code{"tte_BUDS"} object.
#'
#' @param x Object of class \code{"tte_BUDS"}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.tte_BUDS <- function(x, ...) {
  inputs <- x$inputs
  meta <- x$meta

  has_interval <- !identical(inputs$S0L, inputs$S0U)
  cat("\n",
      "Two-Stage Phase II Designs for Time-to-Event Endpoints\n",
      sep = "")
  if (identical(inputs$input_scale_null, "hazard")) {
    cat(sprintf("\n Null hazard at x0 (lambda0): %g\n", inputs$lambda0_resolved))
    if (has_interval) {
      cat(sprintf(" Benchmark range [lambda0L, lambda0U]: [%g, %g]\n",
                  inputs$lambda0L_resolved, inputs$lambda0U_resolved))
    }
  } else {
    cat(sprintf("\n Null survival at x0 (S0): %g\n", inputs$S0))
    if (has_interval) {
      cat(sprintf(" Benchmark range [S0L, S0U]: [%g, %g]\n", inputs$S0L, inputs$S0U))
    }
  }
  if (identical(inputs$input_scale_alt, "lambda1") && !is.null(inputs$lambda1)) {
    cat(sprintf(" Target hazard at x0 (lambda1): %g\n", inputs$lambda1))
  } else {
    cat(sprintf(" Target survival at x0 (S1): %g\n", inputs$S1))
    cat(sprintf(" Hazard ratio (H0/H1): %.3f\n", inputs$HR))
  }
  cat(sprintf(" Error rates: alpha = %.2f; beta = %.2f\n", inputs$alpha, inputs$beta))
  cat(sprintf(" Follow-up / accrual: x0 = %.2f; rate = %.2f\n", inputs$x0, inputs$rate))
  cat(sprintf(" Maximum allowable sample size: %d\n\n", inputs$n.ub))
  if (has_interval) {
    cat(" Reported EN/PET use S0 for the classic r-KJ comparator and S0U for BUDS rows.\n\n")
  }
  if (length(inputs$robust_objective) > 0) {
    cat(sprintf(" BUDS objectives: %s\n", paste(inputs$robust_objective, collapse = ", ")))
  }
  if (!is.null(meta$n_total)) {
    cat(sprintf(" Candidate counts: %d admissible, %d unique, %d robust, %d pooled\n",
                meta$n_total, meta$n_unique, meta$n_robust, meta$n_pooled))
  }
  if (!is.null(meta$execution_time)) {
    cat(sprintf(" Elapsed: %s\n\n", meta$execution_time))
  } else {
    cat("\n")
  }
  df_t <- as.data.frame(t(x$designs), check.names = FALSE)
  int_cols <- c("n1", "n2")
  int_cols <- intersect(int_cols, names(df_t))
  prob_cols <- c("PET", "sup_alpha")
  prob_cols <- intersect(prob_cols, names(df_t))
  en_cols <- c("EN")
  en_cols <- intersect(en_cols, names(df_t))
  other_num_cols <- c("S0_designed", "c1", "DA1", "c2", "DA2")
  other_num_cols <- intersect(other_num_cols, names(df_t))

  df_t[int_cols] <- lapply(df_t[int_cols], function(col) as.integer(round(col)))
  df_t[prob_cols] <- lapply(df_t[prob_cols], function(col) formatC(col, format = "f", digits = 4))
  df_t[en_cols] <- lapply(df_t[en_cols], function(col) formatC(col, format = "f", digits = 2))
  df_t[other_num_cols] <- lapply(df_t[other_num_cols], function(col) formatC(col, format = "f", digits = 3))
  print(df_t, row.names = TRUE)
  invisible(x)
}
