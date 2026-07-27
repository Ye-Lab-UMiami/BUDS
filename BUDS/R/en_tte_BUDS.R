#' Expected Sample Size for a Two-Stage r-KJ BUDS
#'
#' Computes the approximate expected sample size for a restricted Kwak-Jung
#' design evaluated at survival probability \code{s} at time \code{x0}.
#'
#' @param n1 Stage-1 sample size.
#' @param c1 Stage-1 futility boundary.
#' @param n2 Total sample size if the trial continues to the final analysis.
#' @param c2 Final futility boundary.
#' @param DA1 Analysis time for stage 1.
#' @param DA2 Analysis time for the final analysis.
#' @param S0_design Survival probability at \code{x0} used to construct the
#'   design.
#' @param s Survival probability at \code{x0} where EN is evaluated.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#'
#' @return Numeric expected sample size.
en_tte_BUDS <- function(n1, c1, n2, c2, DA1, DA2,
                       S0_design, s, x0, rate) {
  ev <- .tte_eval_candidate_point(
    n1 = n1, c1 = c1, n2 = n2, c2 = c2,
    DA1 = DA1, DA2 = DA2,
    x0 = x0, rate = rate,
    lambda_d = -log(S0_design) / x0,
    lambda_eval = -log(s) / x0
  )
  ev$en
}
