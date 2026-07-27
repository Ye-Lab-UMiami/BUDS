#' Summarize a Survival r-KJ under BUDS
#'
#' Produces a prose summary for one selected design from a
#' \code{"tte_BUDS"} object.
#'
#' @param object Object of class \code{"tte_BUDS"}.
#' @param design_name Name of the design to summarize.
#' @param ... Unused.
#'
#' @return An object of class \code{"tte_BUDS_summary"}.
#' @export
summary.tte_BUDS <- function(object, design_name, ...) {
  
  design_label_map <- function(name) {
    if (name == "r-KJ") return("restricted-Kwak and Jung")
    name
  }
  
  if (missing(design_name) || is.null(design_name)) {
    stop("You must specify a design_name.")
  }
  if (!is.character(design_name) || length(design_name) != 1) {
    stop("design_name must be a single character string.")
  }
  
  design_df <- .tte_design_matrix_to_df(object$designs)
  inputs <- object$inputs
  
  if (!design_name %in% design_df$design) {
    stop("design_name must match one of: ",
         paste(design_df$design, collapse = ", "))
  }
  
  row <- design_df[design_df$design == design_name, , drop = FALSE][1, ]
  
  n1 <- as.integer(row$n1)
  n  <- as.integer(row$n2)
  c1 <- as.numeric(row$c1)
  c2 <- as.numeric(row$c2)
  DA1 <- as.numeric(row$DA1)
  DA2 <- as.numeric(row$DA2)
  label <- design_label_map(row$design)
  is_BUDS <- grepl("^BUDS", label)
  
  if (identical(inputs$input_scale_alt, "lambda1") && !is.null(inputs$lambda1)) {
    S1 <- inputs$S1
    lambda1 <- inputs$lambda1
  } else if (!is.null(inputs$S1)) {
    S1 <- inputs$S1
    lambda1 <- if (!is.null(inputs$x0)) -log(inputs$S1) / inputs$x0 else NULL
  } else if (!is.null(inputs$HR)) {
    lambda0 <- -log(inputs$S0) / inputs$x0
    lambda1 <- lambda0 / inputs$HR
    S1 <- exp(-lambda1 * inputs$x0)
  } else {
    stop("Must provide either S1 or HR in inputs.")
  }
  
  S0_planned <- if ("S0_designed" %in% names(row) && !is.na(row$S0_designed)) {
    as.numeric(row$S0_designed)
  } else if (!is.null(inputs$S0)) {
    inputs$S0
  } else if (!is.null(inputs$S0U)) {
    inputs$S0U
  } else {
    inputs$S0L
  }
  
  lambda0_planned <- -log(S0_planned) / inputs$x0
  
  use_hazard_scale <- identical(inputs$input_scale_null, "hazard") &&
    identical(inputs$input_scale_alt, "lambda1")
  
  if (use_hazard_scale) {
    hyp <- sprintf("H0: lambda(x0) >= %.2f versus H1: lambda(x0) < %.2f",
                   lambda0_planned, lambda0_planned)
  } else {
    hyp <- sprintf("H0: S(x0) <= %.2f versus H1: S(x0) > %.2f",
                   S0_planned, S0_planned)
  }
  
  benchmark_txt <- if (!is.null(inputs$S0L) &&
                       !is.null(inputs$S0U) &&
                       inputs$S0L < inputs$S0U) {
    if (use_hazard_scale) {
      sprintf(
        " The benchmark range used for design selection is [%.2f, %.2f] on the hazard-rate scale, and the target alternative hazard rate is lambda1 = %.2f.",
        inputs$lambda0L_resolved, inputs$lambda0U_resolved, lambda1
      )
    } else {
      sprintf(
        " The benchmark range used for design selection is [%.2f, %.2f] on the survival-probability scale, and the target survival probability is S1 = %.2f.",
        inputs$S0L, inputs$S0U, S1
      )
    }
  } else {
    if (use_hazard_scale) {
      sprintf(" The target alternative hazard rate is lambda1 = %.2f.", lambda1)
    } else {
      sprintf(" The target survival probability is S1 = %.2f.", S1)
    }
  }
  
  EN  <- as.numeric(row$EN)
  PET <- as.numeric(row$PET)
  
  if (is_BUDS && !is.null(inputs$S0L) && !is.null(inputs$S0U) && inputs$S0L < inputs$S0U) {
    en_pet_txt <- sprintf(
      " For reporting, EN and PET are evaluated at the planning benchmark S(x0) = %.2f; the expected sample size for this design is %.2f, with a probability of early termination of %.3f.",
      S0_planned, EN, PET
    )
  } else {
    en_pet_txt <- sprintf(
      " For reporting, EN and PET are evaluated at S(x0) = %.2f; the expected sample size for this design is %.2f, with a probability of early termination of %.3f.",
      S0_planned, EN, PET
    )
  }
  
  endpoint_txt <- if (use_hazard_scale) {
    sprintf("the hazard rate at the clinically meaningful time point x0 = %.2f", inputs$x0)
  } else {
    sprintf("the survival probability at the clinically meaningful time point x0 = %.2f", inputs$x0)
  }
  
  intro_txt <- if (grepl("^BUDS", label)) {
    sprintf(
      "A %s-selected two-stage phase II single-arm clinical trial design is defined to evaluate whether %s supports continuation to the next phase of the clinical trial (%s).%s ",
      label, endpoint_txt, hyp, benchmark_txt
    )
  } else {
    sprintf(
      "A two-stage phase II single-arm clinical trial design, %s, is defined to evaluate whether %s supports continuation to the next phase of the clinical trial (%s).%s ",
      label, endpoint_txt, hyp, benchmark_txt
    )
  }
  
  txt <- sprintf(
    paste0(
      "%s",
      "The design assumes a uniform accrual rate of %.2f subjects per unit of time. ",
      "The design is specified with a type I error rate of %.2f and power of %.0f%%. ",
      "The total number of subjects required if the trial continues to the second stage is %d, with %d subjects enrolled in the first stage and an additional %d subjects enrolled in the second stage, if necessary.%s ",
      "The first interim analysis is conducted at time %.2f. ",
      "After the first stage, the trial will be terminated early for futility if the stage-1 test statistic exceeds %.4f. ",
      "Otherwise, the study will continue to the second stage. ",
      "At the conclusion of the trial, the treatment will be considered promising if the final test statistic does not exceed %.4f."
    ),
    intro_txt,
    inputs$rate,
    inputs$alpha,
    (1 - inputs$beta) * 100,
    n,
    n1,
    n - n1,
    en_pet_txt,
    DA1,
    c1,
    c2
  )
  
  structure(txt, class = "tte_BUDS_summary")
}

#' Print a Summary of a Robust r-KJ Design under BUDS
#'
#' @param x Object of class \code{"tte_BUDS_summary"}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.tte_BUDS_summary <- function(x, ...) {
  cat(x, "\n")
  invisible(x)
}
