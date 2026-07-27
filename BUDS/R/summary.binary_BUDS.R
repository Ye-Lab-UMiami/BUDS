#' Summarize a Binary Simon Design under BUDS
#'
#' Produces a prose summary for one selected design from a
#' \code{"binary_BUDS"} object.
#'
#' @param object Object of class \code{"binary_BUDS"}.
#' @param design_name Name of the design to summarize.
#' @param ... Unused.
#'
#' @return An object of class \code{"binary_BUDS_summary"}.
#' @export
summary.binary_BUDS <- function(object, design_name, ...) {
  
  design_label_map <- function(name) {
    if (name == "Optimal") return("Simon (Optimal)")
    if (name == "Minimax") return("Simon (Minimax)")
    if (name == "Balanced") return("Balanced")
    name
  }
  
  if (missing(design_name) || is.null(design_name)) {
    stop("You must specify a design_name.")
  }
  if (!is.character(design_name) || length(design_name) != 1) {
    stop("design_name must be a single character string. You provided: ",
         paste(design_name, collapse = ", "))
  }
  
  design_matrix <- object$designs
  inputs <- object$inputs
  
  if (!design_name %in% colnames(design_matrix)) {
    stop("design_name must match one of: ",
         paste(colnames(design_matrix), collapse = ", "))
  }
  
  col <- design_matrix[, design_name]
  n1 <- as.integer(col["n1"])
  r1 <- as.integer(col["r1"])
  n  <- as.integer(col["n"])
  r  <- as.integer(col["r"])
  label <- design_label_map(design_name)
  is_BUDS <- grepl("^BUDS", label)
  
  p0_planned <- if (!is.null(inputs$p0)) {
    inputs$p0
  } else if (!is.null(inputs$p0U)) {
    inputs$p0U
  } else {
    inputs$p0L
  }
  
  hyp <- sprintf("H0: p <= %.2f versus H1: p > %.2f",
                 p0_planned, p0_planned)
  
  benchmark_txt <- if (!is.null(inputs$p0L) &&
                       !is.null(inputs$p0U) &&
                       inputs$p0L < inputs$p0U) {
    sprintf(" The benchmark range used for design selection is [%.2f, %.2f], and the target response probability is p1 = %.2f.",
            inputs$p0L, inputs$p0U, inputs$p1)
  } else {
    sprintf(" The target response probability is p1 = %.2f.", inputs$p1)
  }
  
  EN <- as.numeric(col["EN"])
  PET <- as.numeric(col["PET"])
  
  if (is_BUDS && !is.null(inputs$p0L) && !is.null(inputs$p0U) && inputs$p0L < inputs$p0U) {
    en_pet_txt <- sprintf(
      " At the upper bound of the benchmark range (p = %.2f), the expected sample size for this design is %.2f, with a probability of early termination of %.3f.",
      inputs$p0U, EN, PET
    )
  } else {
    en_pet_txt <- sprintf(
      " The expected sample size for this design is %.2f, with a probability of early termination of %.3f.",
      EN, PET
    )
  }
  
  intro_txt <- if (grepl("^BUDS", label)) {
    sprintf(
      "A %s-selected two-stage phase II single-arm clinical trial design is defined to evaluate whether the response probability supports continuation to the next phase of the clinical trial (%s).%s ",
      label, hyp, benchmark_txt
    )
  } else {
    sprintf(
      "A two-stage phase II single-arm clinical trial design, %s, is defined to evaluate whether the response probability supports continuation to the next phase of the clinical trial (%s).%s ",
      label, hyp, benchmark_txt
    )
  }
  
  txt <- sprintf(
    paste0(
      "%s",
      "The design is specified with a type I error rate of %.2f and power of %.0f%%. ",
      "The total number of subjects required if the trial continues to the second stage is %d, with %d subjects enrolled in the first stage and an additional %d subjects enrolled in the second stage, if necessary.%s ",
      "After the first stage, the trial will be terminated early for futility if %d or fewer responses are observed. ",
      "Otherwise, the study will continue to the second stage. ",
      "At the conclusion of the trial, the treatment will be considered promising if more than %d of the %d subjects respond."
    ),
    intro_txt,
    inputs$alpha,
    (1 - inputs$beta) * 100,
    n,
    n1,
    n - n1,
    en_pet_txt,
    r1,
    r,
    n
  )
  
  structure(txt, class = "binary_BUDS_summary")
}

#' Print a Summary of a Binary Simon BUDS Design
#'
#' @param x Object of class \code{"binary_BUDS_summary"}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.binary_BUDS_summary <- function(x, ...) {
  cat(x, "\n")
  invisible(x)
}