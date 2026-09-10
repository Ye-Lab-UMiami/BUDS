#' BUDS for Simon Two-Stage Single-Arm Phase II Trials 
#' 
#' Enumerate admissible Simon two-stage designs under type-I error and power
#' constraints, including classical Optimal, Minimax, and Balanced designs,
#' with optional benchmark uncertainty design selection (BUDS).
#'
#' The core search engine is implemented in C++ for high performance.
#'
#' @param alpha Type-I error rate.
#' @param power Desired power (1 - beta).
#' @param p0 Null response probability.
#' @param p1 Alternative response probability.
#' @param n.lb Lower bound for total sample size search. If NULL, auto-derived.
#' @param n.ub Upper bound for total sample size search.
#' @param p0L Lower bound of interval null.
#' @param p0U Upper bound of interval null.
#' @param pgrid_points Number of grid points for interval null evaluation.
#' @param robust_objective Robust design criterion. Publicly supported defaults
#'   are \code{"least_regret"} and \code{"avg_en"}. Additional internal
#'   criteria may still be available.
#' @param robust_ref Reference design for inflation constraints: "Minimax" or "Optimal".
#' @param robust_max_inflation Maximum allowed inflation over reference N.
#' @param robust_N_cap Absolute cap on allowable sample size.
#' @param robust_auto_expand Logical; expand N cap if infeasible.
#' @param robust_expand_step Inflation step size for auto-expansion.
#' @param robust_expand_max Maximum allowable inflation.
#'
#' @return An object of class \code{"binary_BUDS"} containing:
#' \describe{
#'   \item{designs}{Matrix of design parameters.}
#'   \item{inputs}{List of input parameters.}
#'   \item{meta}{Execution metadata.}
#' }
#'
#' @examples
#' design <- binary_BUDS(alpha = 0.05, power = 0.8,
#'              p0 = 0.10, p1 = 0.25,
#'              p0L = 0.05, p0U = 0.15,
#'              n.ub = 150,
#'              robust_objective = c("least_regret", "avg_en"))
#' design
#' summary(design, "Optimal")
#'
#' @importFrom stats dbinom pbinom qnorm
#' @export
binary_BUDS <- function(alpha, power, p0 = NULL, p1,
                      n.lb = NULL, n.ub = 150,
                      p0L = NULL, p0U = NULL,
                      pgrid_points = NULL,
                      robust_objective = c("least_regret", "avg_en"),
                      robust_ref = c("Minimax", "Optimal"),
                      robust_max_inflation = NULL,
                      robust_N_cap = NULL,
                      robust_auto_expand = FALSE,
                      robust_expand_step = 0.05,
                      robust_expand_max = 2.0) {

  valid_objectives <- c("least_regret", "least_en", "avg_en", "avg_regret", "min_N")
  default_display  <- c("least_regret", "avg_en")

  if (is.null(robust_objective) || identical(robust_objective, valid_objectives)) {
    robust_objective <- default_display
  } else {
    robust_objective <- intersect(robust_objective, valid_objectives)
    if (length(robust_objective) == 0) {
      stop("No valid robust_objective specified")
    }
  }
  robust_objective <- unique(valid_objectives[valid_objectives %in% robust_objective])
  robust_ref <- match.arg(robust_ref, c("Minimax", "Optimal"), several.ok = FALSE)

  if (is.null(p0) && !is.null(p0L) && !is.null(p0U)) {
    p0 <- (p0L + p0U) / 2
  }
  if (is.null(p0)) {
    stop("Supply p0, or supply both p0L and p0U so p0 can be derived as the midpoint.")
  }
  if (is.null(p0L)) p0L <- p0
  if (is.null(p0U)) p0U <- p0
  if (p0L > p0U) stop("Require p0L <= p0U")

  if (is.null(pgrid_points)) {
    if (p0L == p0U) {
      pgrid_points <- 1L
    } else {
      width <- p0U - p0L
      target_step <- if (width <= 0.05) 0.0025 else if (width <= 0.10) 0.005 else 0.01
      pgrid_points <- as.integer(max(11, min(201, floor(width / target_step) + 1)))
    }
  }

  if (is.null(n.lb)) n.lb <- -1L
  if (is.null(robust_max_inflation)) robust_max_inflation <- -1
  if (is.null(robust_N_cap)) robust_N_cap <- -1L

  result <- binary_BUDS_cpp(
    alpha = alpha,
    power = power,
    p0 = p0,
    p1 = p1,
    n_lb = as.integer(n.lb),
    n_ub = as.integer(n.ub),
    p0L = p0L,
    p0U = p0U,
    pgrid_points = as.integer(pgrid_points),
    robust_objective = robust_objective,
    robust_ref = robust_ref,
    robust_max_inflation = robust_max_inflation,
    robust_N_cap = as.integer(robust_N_cap),
    robust_auto_expand = robust_auto_expand,
    robust_expand_step = robust_expand_step,
    robust_expand_max = robust_expand_max
  )
  class(result) <- c("binary_BUDS", "BUDS")

  is_single_point <- identical(p0L, p0U)
  if (is_single_point && !is.null(result$designs)) {
    keep <- intersect(c("Optimal", "Minimax", "Balanced"), colnames(result$designs))
    result$designs <- result$designs[, keep, drop = FALSE]
    objective_map <- c(
      "least_regret" = "BUDS (Least Regret)",
      "avg_en" = "BUDS (Average EN)"
    )
    objective_keep <- robust_objective[robust_objective %in% names(objective_map)]
    if ("Optimal" %in% colnames(result$designs) && length(objective_keep) > 0) {
      BUDS_block <- sapply(objective_keep, function(obj) result$designs[, "Optimal"], simplify = "matrix")
      colnames(BUDS_block) <- unname(objective_map[objective_keep])
      result$designs <- cbind(result$designs, BUDS_block)
    }
    if (!is.null(result$inputs)) {
      result$inputs$robust_objective <- objective_keep
      result$inputs$p0L <- p0
      result$inputs$p0U <- p0
    }
  }

  result
}

#' Print Binary Simon BUDS Designs
#'
#' Displays the input settings and selected binary two-stage Simon designs
#' stored in a \code{"binary_BUDS"} object.
#'
#' @param x Object of class \code{"binary_BUDS"}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.binary_BUDS <- function(x, ...) {
  inputs <- x$inputs
  has_interval <- !identical(inputs$p0L, inputs$p0U)
  cat("\n",
      "Two-Stage Phase II Designs for Tests of One Proportion\n",
      sep = "")
  cat(sprintf("\n Null response rate (p0): %g\n", inputs$p0))
  if (has_interval) {
    cat(sprintf(" Benchmark range [p0L, p0U]: [%g, %g]\n", inputs$p0L, inputs$p0U))
  }
  cat(sprintf(" Target response rate (p1): %g\n", inputs$p1))
  cat(sprintf(" Error rates: alpha = %.2f; beta = %.2f\n", inputs$alpha, inputs$beta))
  cat(sprintf(" Maximum allowable sample size: %d\n\n", inputs$n.ub))
  
  df_t <- as.data.frame(t(x$designs), check.names = FALSE)
  int_cols <- c("r1", "n1", "r", "n")
  int_cols <- intersect(int_cols, names(df_t))
  prob_cols <- c("PET", "alpha", "beta", "alpha1", "alpha2", "beta1", "beta2")
  prob_cols <- intersect(prob_cols, names(df_t))
  en_cols <- c("EN")
  en_cols <- intersect(en_cols, names(df_t))
  
  df_t[int_cols] <- lapply(df_t[int_cols], function(col) as.integer(col))
  df_t[prob_cols] <- lapply(df_t[prob_cols], function(col) {
    formatC(col, format = "f", digits = 4)
  })
  df_t[en_cols] <- lapply(df_t[en_cols], function(col) {
    formatC(col, format = "f", digits = 2)
  })
  print(df_t, row.names = TRUE)
  invisible(x)
}
