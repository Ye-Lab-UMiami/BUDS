#' Power for a Two-Stage r-KJ Design
#'
#' Computes the approximate power of a restricted Kwak-Jung two-stage design
#' at target survival probability \code{S1}, using the averaged-hazard
#' approximation of Kwak and Jung.
#'
#' @param n1 Stage-1 sample size.
#' @param c1 Stage-1 futility boundary.
#' @param n2 Total sample size if the trial continues to the final analysis.
#' @param c2 Final futility boundary.
#' @param DA1 Analysis time for stage 1.
#' @param DA2 Analysis time for the final analysis.
#' @param S0_design Survival probability at \code{x0} used to construct the
#'   design.
#' @param S1 Alternative survival probability at \code{x0}.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#'
#' @return Numeric power.
power_tte_BUDS <- function(n1, c1, n2, c2, DA1, DA2,
                          S0_design, S1, x0, rate) {
  .tte_power_candidate(
    n1 = n1, c1 = c1, n2 = n2, c2 = c2,
    DA1 = DA1, DA2 = DA2,
    x0 = x0, rate = rate,
    lambda_d = -log(S0_design) / x0,
    lambda_alt = -log(S1) / x0
  )
}
