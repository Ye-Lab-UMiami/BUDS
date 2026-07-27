#' BUDS: Benchmark Uncertainty Design Selection for Two-Stage Single-Arm Phase II Trials
#'
#' The BUDS package provides tools for designing and reporting two-stage
#' single-arm phase II clinical trials under point-null and interval-null
#' benchmarks. It includes binary Simon designs and restricted Kwak-Jung
#' time-to-event designs, with Interval-Null Robust (BUDS) selection criteria.
#'
#' Main user-facing functions are:
#' \describe{
#'   \item{\code{\link{binary_BUDS}}}{Search binary Simon with BUDS.}
#'   \item{\code{\link{tte_BUDS}}}{Search time-to-event restricted Kwak-Jung with BUDS.}
#'   \item{\code{\link{opchar_binary_BUDS}} and \code{\link{opchar_tte_BUDS}}}{Evaluate operating characteristics.}
#'   \item{\code{\link{table_design_opchar}}, \code{\link{table_design_rules}}, and \code{\link{table_design_inputs}}}{Create reporting tables.}
#'   \item{\code{\link{plot_rejection_rate}}}{Plot type I error across the benchmark range.}
#'   \item{\code{\link{export_design_html}}}{Export a self-contained HTML design report.}
#' }
#'
#' @keywords internal
"_PACKAGE"
