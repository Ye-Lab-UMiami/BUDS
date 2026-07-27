#' Expected Sample Size for BUDS
#'
#' Computes the expected sample size (EN) at a given response rate p,
#' accounting for the probability of early termination after stage 1.
#'
#' @param n1 Stage-1 sample size.
#' @param r1 Stage-1 futility cutoff. The trial stops early for futility if
#'   \code{r1} or fewer responses are observed during stage 1.
#' @param n Total sample size if the trial continues to the final analysis.
#' @param p Response probability where EN is evaluated.
#'
#' @return Numeric expected sample size.
#' @importFrom stats pbinom
en_binary_BUDS <- function(n1, r1, n, p) {
  n1 + (1 - pbinom(r1, n1, p)) * (n - n1)
}
