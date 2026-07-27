#' Operating Characteristics for r-KJ and TTE BUDS
#'
#' Calculates operating characteristics either from an `tte_BUDS` object or
#' from manually supplied restricted-Kwak and Jung decision boundaries.
#'
#' @param design Optional object of class `tte_BUDS`.
#' @param n1,c1,n2,c2,DA1,DA2 Manual time-to-event decision boundaries used
#'   when `design` is `NULL`.
#' @param S0_designed Planning null survival used to build the design.
#' @param S0,S1 Survival probabilities at `x0` under the null and alternative.
#' @param x0 Clinically meaningful follow-up time.
#' @param rate Accrual rate.
#' @param alpha_target Target one-sided type I error rate.
#' @param S0L,S0U Optional benchmark bounds. When omitted, a fixed benchmark is
#'   used.
#' @param grid_points Number of grid points used to evaluate operating characteristics.
#'  under BUDS
#' @param design_name Optional character vector of design names to keep when a
#'   `design` object is supplied.
#'
#' @return An object of class `opchar`.
#' @export
opchar_tte_BUDS <- function(design = NULL,
                           n1 = NULL, c1 = NULL, n2 = NULL, c2 = NULL, DA1 = NULL, DA2 = NULL,
                           S0_designed = NULL, S0 = NULL, S1 = NULL,
                           x0 = NULL, rate = NULL, alpha_target = NULL,
                           S0L = NULL, S0U = NULL,
                           grid_points = 101, design_name = NULL) {
  
  # =============================
  # CASE 1: design object
  # =============================
  if (!is.null(design)) {
    if (!inherits(design, "tte_BUDS")) {
      stop("design must be an tte_BUDS object")
    }
    design_df <- .tte_design_matrix_to_df(design$designs)
    inputs <- design$inputs
    if (is.null(S0)) S0 <- inputs$S0
    if (is.null(S1)) S1 <- inputs$S1
    if (is.null(S0L)) S0L <- inputs$S0L
    if (is.null(S0U)) S0U <- inputs$S0U
    if (is.null(alpha_target)) alpha_target <- inputs$alpha
    if (!is.null(design_name)) {
      design_df <- design_df[design_df$design %in% design_name, , drop = FALSE]
    }
    results <- lapply(seq_len(nrow(design_df)), function(i) {
      row <- design_df[i, , drop = FALSE]
      res <- .opchar_tte(n1 = as.integer(row$n1),
                         c1 = row$c1,
                         n2 = as.integer(row$n2),
                         c2 = row$c2,
                         DA1 = row$DA1,
                         DA2 = row$DA2,
                         S0_design = as.numeric(row$S0_designed),
                         S0 = S0,
                         S1 = S1,
                         x0 = inputs$x0,
                         rate = inputs$rate,
                         alpha_target = alpha_target,
                         S0L = S0L,
                         S0U = S0U,
                         grid_points = grid_points)
      res$Design <- row$design
      res$S0_designed <- as.numeric(row$S0_designed)
      res$n1 <- as.integer(row$n1)
      res$c1 <- row$c1
      res$n2 <- as.integer(row$n2)
      res$c2 <- row$c2
      res$DA1 <- row$DA1
      res$DA2 <- row$DA2
      res
    })
    out <- dplyr::bind_rows(results)
    out <- .standardize_tte(out)
    out <- cbind(Design = design_df$design, out)
    for (col in names(out)) {
      if (!is.numeric(out[[col]])) next
      if (col %in% c("S0_designed", "alpha_target")) {
        out[[col]] <- round(out[[col]], 2)
      } else if (grepl("^EN_", col) || col == "avg_EN") {
        out[[col]] <- round(out[[col]], 2)
      } else {
        out[[col]] <- round(out[[col]], 4)
      }
    }
    return(structure(
      list(results = out,
           inputs = list(S0 = S0,
                         S1 = S1,
                         S0L = S0L,
                         S0U = S0U,
                         alpha_target = alpha_target,
                         beta_target = inputs$beta,
                         x0 = inputs$x0,
                         rate = inputs$rate,
                         source = "design",
                         endpoint = "tte"),
           design_inputs = design$inputs),
      class = "opchar"))
  }
  
  # =============================
  # CASE 2: manual input
  # =============================
  if (is.null(S0_designed)) S0_designed <- S0
  if (is.null(S0)) S0 <- S0_designed
  missing_args <- c(n1 = is.null(n1),
                    c1 = is.null(c1),
                    n2 = is.null(n2),
                    c2 = is.null(c2),
                    DA1 = is.null(DA1),
                    DA2 = is.null(DA2),
                    S0_designed = is.null(S0_designed),
                    S1 = is.null(S1),
                    x0 = is.null(x0),
                    rate = is.null(rate),
                    alpha_target = is.null(alpha_target))
  if (any(missing_args)) {
    stop("Must provide manual inputs: ",
         paste(names(missing_args)[missing_args], collapse = ", "))
  }
  res <- .opchar_tte(n1 = n1,
                     c1 = c1,
                     n2 = n2,
                     c2 = c2,
                     DA1 = DA1,
                     DA2 = DA2,
                     S0_design = S0_designed,
                     S0 = S0,
                     S1 = S1,
                     x0 = x0,
                     rate = rate,
                     alpha_target = alpha_target,
                     S0L = S0L,
                     S0U = S0U,
                     grid_points = grid_points)
  res$S0_designed <- S0_designed
  res$n1 <- n1
  res$c1 <- c1
  res$n2 <- n2
  res$c2 <- c2
  res$DA1 <- DA1
  res$DA2 <- DA2
  res <- .standardize_tte(res)
  for (col in names(res)) {
    if (!is.numeric(res[[col]])) next
    if (col %in% c("S0_designed", "alpha_target")) {
      res[[col]] <- round(res[[col]], 2)
    } else if (grepl("^EN_", col) || col == "avg_EN") {
      res[[col]] <- round(res[[col]], 2)
    } else {
      res[[col]] <- round(res[[col]], 4)
    }
  }
  
  return(structure(
    list(results = res,
         inputs = list(S0 = S0,
                       S1 = S1,
                       S0L = S0L,
                       S0U = S0U,
                       alpha_target = alpha_target,
                       x0 = x0,
                       rate = rate,
                       source = "manual",
                       endpoint = "tte")),
    class = "opchar"))
}

#' Standardize Time-to-Event Operating-Characteristic Columns
#'
#' Internal helper that orders TTE operating
#' characteristic columns before printing, tabulation, or summarization
#' under a fixed benchmark or benchmark range.
#'
#' @param df Data frame of TTE operating characteristics.
#'
#' @return A data frame with columns ordered for downstream display.
#' @keywords internal
.standardize_tte <- function(df) {
  base_cols <- c("S0_designed","n1","c1","n2","c2","DA1","DA2",
                 "EN_S0","PET_S0",
                 "alpha_target","alpha_at_S0","power_at_S1")
  interval_cols <- c("EN_S0L","EN_S0U","avg_EN",
                     "PET_S0L","PET_S0U",
                     "alpha_at_S0L","alpha_at_S0U","sup_alpha")
  cols <- base_cols
  if (any(interval_cols %in% names(df))) {
    cols <- c("S0_designed","n1","c1","n2","c2","DA1","DA2",
              "EN_S0L","EN_S0","EN_S0U","avg_EN",
              "PET_S0L","PET_S0","PET_S0U",
              "alpha_target",
              "alpha_at_S0L","alpha_at_S0","alpha_at_S0U","sup_alpha",
              "power_at_S1")
    cols <- cols[cols %in% names(df)]
  }
  df <- df[, cols, drop = FALSE]
  return(df)
}

#' Bivariate Normal Integral for r-KJ Calculations
#'
#' Internal numerical integral used in restricted Kwak-Jung rejection and power
#' calculations.
#'
#' @param upper Upper integration limit.
#' @param c1 Stage-1 transformed boundary.
#' @param rho Correlation between stage-wise test statistics.
#'
#' @return Numeric integral value, or \code{NA_real_} for invalid inputs.
#' @keywords internal
.tte_bvn_integral <- function(upper, c1, rho) {
  if (!is.finite(upper) || !is.finite(c1) || !is.finite(rho) || rho <= 0 || rho >= 1) {
    return(NA_real_)
  }
  lower <- max(upper - 8, -8)
  sr <- sqrt(1 - rho^2)
  stats::integrate(
    function(z) stats::dnorm(z) * stats::pnorm((c1 - rho * z) / sr),
    lower = lower,
    upper = upper,
    subdivisions = 200L,
    rel.tol = 1e-8
  )$value
}

#' Restricted Follow-Up Integral Component
#'
#' Internal helper for the restricted Kwak-Jung information calculation.
#'
#' @param DA1 Analysis time.
#' @param ta Accrual duration.
#' @param hz Hazard rate.
#' @param BORNE Restricted integration bound.
#'
#' @return Numeric information component.
#' @keywords internal
.tte_fct <- function(DA1, ta, hz, BORNE) {
  e <- exp(-hz * BORNE)
  (DA1 * (1 - e) + BORNE * e + e / hz - 1 / hz) / (ta * hz)
}

#' Stage-Wise Restricted Follow-Up Integral
#'
#' Internal helper that evaluates the restricted follow-up information integral
#' for one analysis time.
#'
#' @param lambda_eval Evaluation hazard rate.
#' @param DAk Analysis time for the stage.
#' @param ta Accrual duration.
#' @param x0 Restricted follow-up time.
#'
#' @return Numeric information integral.
#' @keywords internal
.tte_stage_integral <- function(lambda_eval, DAk, ta, x0) {
  if (!is.finite(lambda_eval) || lambda_eval <= 0) return(NA_real_)
  if (abs(DAk - (ta + x0)) < 1e-10) {
    return((1 - exp(-lambda_eval * x0)) / lambda_eval)
  }
  .tte_fct(DAk, ta, lambda_eval, min(DAk, x0))
}

#' Evaluate One r-KJ Candidate at One Null Value
#'
#' Internal helper that evaluates type I error, PET, and expected sample size
#' for one candidate TTE design at one null value.
#'
#' @param n1,c1,n2,c2,DA1,DA2 Candidate design parameters.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#' @param lambda_d Planning hazard rate.
#' @param lambda_eval Evaluation hazard rate.
#'
#' @return A list with \code{alpha}, \code{pet}, and \code{en}.
#' @keywords internal
.tte_eval_candidate_point <- function(n1, c1, n2, c2, DA1, DA2, x0, rate,
                                      lambda_d, lambda_eval) {
  ta <- n2 / rate
  int1 <- .tte_stage_integral(lambda_eval, DA1, ta, x0)
  int2 <- .tte_stage_integral(lambda_eval, DA2, ta, x0)
  if (!is.finite(int1) || !is.finite(int2) || int1 <= 0 || int2 <= 0) {
    return(list(alpha = NA_real_, pet = NA_real_, en = NA_real_))
  }
  
  sigma01_sq <- lambda_d * int1
  sigma02_sq <- lambda_d * int2
  sigmae1_sq <- lambda_eval * int1
  sigmae2_sq <- lambda_eval * int2
  if (min(sigma01_sq, sigma02_sq, sigmae1_sq, sigmae2_sq) <= 0) {
    return(list(alpha = NA_real_, pet = NA_real_, en = NA_real_))
  }
  
  sigma01 <- sqrt(sigma01_sq)
  sigma02 <- sqrt(sigma02_sq)
  sigmae1 <- sqrt(sigmae1_sq)
  sigmae2 <- sqrt(sigmae2_sq)
  rho <- sigmae1 / sigmae2
  omega1 <- (lambda_eval - lambda_d) * int1
  omega2 <- (lambda_eval - lambda_d) * int2
  
  cbar1 <- (sigma01 * c1 - omega1 * sqrt(n1)) / sigmae1
  cbar2 <- (sigma02 * c2 - omega2 * sqrt(n2)) / sigmae2
  
  alpha_pt <- .tte_bvn_integral(cbar2, cbar1, rho)
  pet <- 1 - stats::pnorm(cbar1)
  en <- n1 + (1 - pet) * (n2 - n1)
  list(alpha = alpha_pt, pet = pet, en = en)
}

#' Evaluate Power for One r-KJ Candidate
#'
#' Internal helper that evaluates approximate power for one candidate TTE
#' design at an alternative hazard rate.
#'
#' @param n1,c1,n2,c2,DA1,DA2 Candidate design parameters.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#' @param lambda_d Planning null hazard rate.
#' @param lambda_alt Alternative hazard rate.
#'
#' @return Numeric power value, or \code{NA_real_} for invalid inputs.
#' @keywords internal
.tte_power_candidate <- function(n1, c1, n2, c2, DA1, DA2, x0, rate,
                                 lambda_d, lambda_alt) {
  ta <- n2 / rate
  int1_alt <- .tte_stage_integral(lambda_alt, DA1, ta, x0)
  int2_alt <- .tte_stage_integral(lambda_alt, DA2, ta, x0)
  if (!is.finite(int1_alt) || !is.finite(int2_alt) || int1_alt <= 0 || int2_alt <= 0) {
    return(NA_real_)
  }
  
  sigma01_sq <- lambda_d * int1_alt
  sigma02_sq <- lambda_d * int2_alt
  omega1 <- (lambda_alt - lambda_d) * int1_alt
  omega2 <- (lambda_alt - lambda_d) * int2_alt
  
  lambda_mid <- (lambda_alt + lambda_d) / 2
  t1_mid <- .tte_stage_integral(lambda_mid, DA1, ta, x0)
  t2_mid <- .tte_stage_integral(lambda_mid, DA2, ta, x0)
  if (!is.finite(t1_mid) || !is.finite(t2_mid) || t1_mid <= 0 || t2_mid <= 0) {
    return(NA_real_)
  }
  
  sigma11 <- sqrt(lambda_mid * t1_mid)
  sigma12 <- sqrt(lambda_mid * t2_mid)
  rho1 <- sigma11 / sigma12
  cbar1 <- (sqrt(sigma01_sq) * c1 - omega1 * sqrt(n1)) / sigma11
  cbar2 <- (sqrt(sigma02_sq) * c2 - omega2 * sqrt(n2)) / sigma12
  .tte_bvn_integral(cbar2, cbar1, rho1)
}

#' Evaluate Time-to-Event Operating Characteristics
#'
#' Internal helper that computes operating characteristics for a single
#' restricted Kwak-Jung design at a single or benchmark interval.
#'
#' @param n1,c1,n2,c2,DA1,DA2 Candidate design parameters.
#' @param S0_design Planning null survival probability.
#' @param S1 Alternative survival probability.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#' @param alpha_target Target type I error rate.
#' @param S0 Optional single benchmark survival probability for reporting.
#' @param S0L,S0U Optional BUDS bounds.
#' @param grid_points Number of grid points for BUDS evaluation.
#'
#' @return A one-row data frame of operating characteristics.
#' @keywords internal
.opchar_tte <- function(
    n1, c1, n2, c2, DA1, DA2,
    S0_design, S1, x0, rate, alpha_target,
    S0 = NULL, S0L = NULL, S0U = NULL,
    grid_points = 101) {
  if (is.null(S0)) S0 <- S0_design
  vals <- c(S0_design, S0, S1, x0, rate, DA1, DA2, n1, c1, n2, c2, alpha_target)
  if (!all(is.finite(vals))) {
    return(data.frame(alpha_target = alpha_target,
                      alpha_at_S0  = NA_real_,
                      power_at_S1  = NA_real_,
                      PET_S0       = NA_real_,
                      EN_S0        = NA_real_))
  }
  pt <- list(alpha = pr_reject_tte_BUDS(n1, c1, n2, c2, DA1, DA2, S0_design, S0, x0, rate),
             pet   = pet_tte_BUDS(n1, c1, n2, c2, DA1, DA2, S0_design, S0, x0, rate),
             en    = en_tte_BUDS(n1, c1, n2, c2, DA1, DA2, S0_design, S0, x0, rate))
  power_at_S1 <- power_tte_BUDS(n1 = n1, c1 = c1, n2 = n2, c2 = c2,
                               DA1 = DA1, DA2 = DA2,
                               S0_design = S0_design, S1 = S1,
                               x0 = x0, rate = rate)
  out <- data.frame(alpha_target = alpha_target,
                    alpha_at_S0  = pmin(pmax(pt$alpha, 0), 1),
                    power_at_S1  = pmin(pmax(power_at_S1, 0), 1),
                    PET_S0       = pt$pet,
                    EN_S0        = pt$en)
  has_interval <- !is.null(S0L) && !is.null(S0U) && S0L != S0U
  if (has_interval) {
    if (S0L > S0U) stop("Require S0L <= S0U")
    eval_one <- function(s) {
      list(alpha = pr_reject_tte_BUDS(n1, c1, n2, c2, DA1, DA2, S0_design, s, x0, rate),
           pet   = pet_tte_BUDS(n1, c1, n2, c2, DA1, DA2, S0_design, s, x0, rate),
           en    = en_tte_BUDS(n1, c1, n2, c2, DA1, DA2, S0_design, s, x0, rate))
    }
    lo <- eval_one(S0L)
    hi <- eval_one(S0U)
    grid <- seq(S0L, S0U, length.out = grid_points)
    grid_eval <- lapply(grid, eval_one)
    a_grid  <- vapply(grid_eval, function(x) x$alpha, numeric(1))
    en_grid <- vapply(grid_eval, function(x) x$en, numeric(1))
    out$alpha_at_S0L <- pmin(pmax(lo$alpha, 0), 1)
    out$alpha_at_S0U <- pmin(pmax(hi$alpha, 0), 1)
    out$sup_alpha    <- pmin(pmax(max(a_grid, na.rm = TRUE), 0), 1)
    out$PET_S0L      <- lo$pet
    out$PET_S0U      <- hi$pet
    out$EN_S0L       <- lo$en
    out$EN_S0U       <- hi$en
    out$avg_EN       <- mean(en_grid, na.rm = TRUE)
  }
  out
}
