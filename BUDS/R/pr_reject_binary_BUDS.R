#' Exact Type I error for a Two-Stage Simon BUDS Design
#'
#' Computes the exact probability of rejecting the null hypothesis at a given
#' response rate p for a two-stage design.
#'
#' @param n1 Stage-1 sample size.
#' @param r1 Stage-1 futility cutoff. The trial stops early for futility if
#'   \code{r1} or fewer responses are observed during stage 1.
#' @param n Total sample size if the trial continues to the final analysis.
#' @param r Final futility cutoff. The treatment is not considered promising if
#'   \code{r} or fewer responses are observed by the end of stage 2.
#' @param p Response probability where type I error is evaluated.
#'
#' @return Numeric type I error.
#' @importFrom stats dbinom pbinom
pr_reject_binary_BUDS <- function(n1, r1, n, r, p) {
  r <- r + 1
  k <- 0:n1
  pmf1    <- dbinom(k, n1, p)
  proceed <- k > r1
  needed2 <- pmax(r - k, 0)
  tail2   <- pbinom(needed2 - 1, n - n1, p, lower.tail = FALSE)
  prob    <- sum(pmf1[proceed] * tail2[proceed])
  prob    <- min(max(prob, 0), 1)
  return(prob)
}
