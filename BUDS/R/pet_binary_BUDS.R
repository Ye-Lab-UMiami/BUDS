#' Probability of Early Termination under a Two-Stage Simon BUDS
#'
#' Computes the probability of early termination (stopping for futility
#' after stage 1) at a given response rate p.
#'
#' @param n1 Stage-1 sample size.
#' @param r1 Stage-1 futility cutoff. The trial stops early for futility if
#'   \code{r1} or fewer responses are observed during stage 1.
#' @param p Response probability where PET is evaluated.
#'
#' @return Numeric probability of early termination.
#' @importFrom stats pbinom
pet_binary_BUDS <- function(n1, r1, p) {
  pbinom(r1, n1, p)
}
