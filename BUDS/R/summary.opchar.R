#' Summarize Operating Characteristics
#'
#' Produces a prose summary from an `opchar` object returned by
#' [opchar_binary_BUDS()] or [opchar_tte_BUDS()].
#'
#' @param object Object of class `opchar`.
#' @param design_name Optional design name when `object` contains multiple
#'   designs.
#' @param ... Unused.
#'
#' @return An object of class `opchar_summary`.
#' @export
summary.opchar <- function(object, design_name = NULL, ...) {
  
  df <- object$results
  inputs <- object$inputs
  source <- if (!is.null(inputs$source)) inputs$source else "manual"
  endpoint <- inputs$endpoint
  
  design_label_map_binary <- function(name) {
    if (name == "Optimal") return("Simon (Optimal)")
    if (name == "Minimax") return("Simon (Minimax)")
    if (name == "Balanced") return("Balanced")
    name
  }
  
  design_label_map_tte <- function(name) {
    if (name == "r-KJ") return("restricted-Kwak and Jung")
    name
  }
  
  range_txt <- function(a, b, digits = 4) {
    vals <- sort(c(as.numeric(a), as.numeric(b)))
    sprintf(paste0("%.", digits, "f to %.", digits, "f"), vals[1], vals[2])
  }
  
  if (source == "manual") {
    row <- df[1, ]
    label <- NULL
  } else {
    if (nrow(df) == 1 && is.null(design_name)) {
      row <- df[1, , drop = FALSE]
    } else {
      if (is.null(design_name) || length(design_name) != 1) {
        stop("You must specify exactly one design_name")
      }
      row <- df[df$Design == design_name, , drop = FALSE]
      if (nrow(row) == 0) stop("design_name must match exactly one design")
    }
    
    if (endpoint == "tte") {
      label <- design_label_map_tte(row$Design[1])
    } else {
      label <- design_label_map_binary(row$Design[1])
    }
  }
  
  if (endpoint == "tte") {
    n1 <- row$n1[1]
    c1 <- row$c1[1]
    DA1 <- row$DA1[1]
    n2 <- row$n2[1]
    c2 <- row$c2[1]
    DA2 <- row$DA2[1]
  } else {
    n1 <- row$n1[1]
    r1 <- row$r1[1]
    n  <- row$n[1]
    r  <- row$r[1]
  }
  
  if (source == "design") {
    design_inputs <- if (!is.null(object$design_inputs)) object$design_inputs else inputs
    
    if (endpoint == "tte") {
      input_scale_null <- if (!is.null(design_inputs$input_scale_null)) {
        tolower(as.character(design_inputs$input_scale_null))
      } else {
        "survival"
      }
      use_hazard_scale <- input_scale_null %in% c("hazard", "lambda", "lambda0")
      
      S0_planned <- if ("S0_designed" %in% names(row) && !is.na(row$S0_designed[1])) {
        as.numeric(row$S0_designed[1])
      } else if (!is.null(design_inputs$S0)) {
        design_inputs$S0
      } else if (!is.null(design_inputs$S0U)) {
        design_inputs$S0U
      } else {
        design_inputs$S0L
      }
      
      lambda0_planned <- if (!is.null(design_inputs$lambda0_resolved)) {
        design_inputs$lambda0_resolved
      } else {
        -log(S0_planned) / design_inputs$x0
      }
      
      lambda1 <- if (!is.null(design_inputs$lambda1)) {
        design_inputs$lambda1
      } else {
        -log(design_inputs$S1) / design_inputs$x0
      }
      
      if (use_hazard_scale) {
        hyp <- sprintf(
          "H0: lambda(x0) >= %.2f versus H1: lambda(x0) < %.2f",
          lambda0_planned, lambda0_planned
        )
        
        target_txt <- sprintf(
          "The target alternative hazard rate used for power evaluation is lambda1 = %.2f.",
          lambda1
        )
        
        benchmark_txt <- if (!is.null(design_inputs$lambda0L_resolved) &&
                             !is.null(design_inputs$lambda0U_resolved) &&
                             design_inputs$lambda0L_resolved != design_inputs$lambda0U_resolved) {
          sprintf(
            "The benchmark range used for design selection and operating-characteristic evaluation is [%.2f, %.2f] on the hazard-rate scale.",
            design_inputs$lambda0L_resolved, design_inputs$lambda0U_resolved
          )
        } else {
          "A single planning benchmark is used."
        }
      } else {
        hyp <- sprintf(
          "H0: S(x0) <= %.2f versus H1: S(x0) > %.2f",
          S0_planned, S0_planned
        )
        
        target_txt <- sprintf(
          "The target alternative survival probability used for power evaluation is S1 = %.2f.",
          design_inputs$S1
        )
        
        benchmark_txt <- if (!is.null(design_inputs$S0L) &&
                             !is.null(design_inputs$S0U) &&
                             design_inputs$S0L != design_inputs$S0U) {
          sprintf(
            "The benchmark range used for design selection and operating-characteristic evaluation is [%.2f, %.2f] on the survival-probability scale.",
            design_inputs$S0L, design_inputs$S0U
          )
        } else {
          "A single planning benchmark is used."
        }
      }
    } else {
      p0_planned <- if (!is.null(design_inputs$p0)) {
        design_inputs$p0
      } else if (!is.null(design_inputs$p0U)) {
        design_inputs$p0U
      } else {
        design_inputs$p0L
      }
      
      hyp <- sprintf(
        "H0: p <= %.2f versus H1: p > %.2f",
        p0_planned, p0_planned
      )
      
      target_txt <- sprintf(
        "The target alternative response probability used for power evaluation is p1 = %.2f.",
        design_inputs$p1
      )
      
      benchmark_txt <- if (!is.null(design_inputs$p0L) &&
                           !is.null(design_inputs$p0U) &&
                           design_inputs$p0L != design_inputs$p0U) {
        sprintf(
          "The benchmark range used for design selection and operating-characteristic evaluation is [%.2f, %.2f].",
          design_inputs$p0L, design_inputs$p0U
        )
      } else {
        "A single planning benchmark is used."
      }
    }
  }
  
  if (endpoint == "tte") {
    if (source == "manual") {
      design_txt <- sprintf(
        paste0(
          "A two-stage phase II single-arm design is defined by the boundaries n1 = %d, c1 = %.4f, DA1 = %.4f, n2 = %d, c2 = %.4f, and DA2 = %.4f, ",
          "with a target type I error rate of %.2f. "
        ),
        n1, c1, DA1, n2, c2, DA2, inputs$alpha_target
      )
    } else {
      intro_txt <- if (grepl("^BUDS", label)) {
        sprintf(
          "The %s-selected two-stage phase II single-arm design is summarized under the planning benchmark specification (%s), where x0 = %.2f. ",
          label, hyp, inputs$x0
        )
      } else {
        sprintf(
          "The two-stage phase II single-arm design, %s, is summarized under the planning benchmark specification (%s), where x0 = %.2f. ",
          label, hyp, inputs$x0
        )
      }
      
      design_txt <- sprintf(
        paste0(
          "%s",
          "%s %s ",
          "The design is specified with a type I error rate of %.2f and power of %.0f%%. ",
          "The design boundaries are n1 = %d, c1 = %.4f, DA1 = %.4f, n2 = %d, c2 = %.4f, DA2 = %.4f. "
        ),
        intro_txt,
        benchmark_txt,
        target_txt,
        inputs$alpha_target,
        (1 - inputs$beta_target) * 100,
        n1, c1, DA1, n2, c2, DA2
      )
    }
  } else {
    if (source == "manual") {
      design_txt <- sprintf(
        paste0(
          "A two-stage phase II single-arm design is defined by the boundaries n1 = %d, r1 = %d, n = %d, and r = %d, ",
          "with a target type I error rate of %.2f. "
        ),
        n1, r1, n, r, inputs$alpha_target
      )
    } else {
      intro_txt <- if (grepl("^BUDS", label)) {
        sprintf(
          "The %s-selected two-stage phase II single-arm design is summarized under the planning benchmark specification (%s). ",
          label, hyp
        )
      } else {
        sprintf(
          "The two-stage phase II single-arm design, %s, is summarized under the planning benchmark specification (%s). ",
          label, hyp
        )
      }
      
      design_txt <- sprintf(
        paste0(
          "%s",
          "%s %s ",
          "The design is specified with a type I error rate of %.2f and power of %.0f%%. ",
          "The design boundaries are n1 = %d, r1 = %d, n = %d, r = %d. "
        ),
        intro_txt,
        benchmark_txt,
        target_txt,
        inputs$alpha_target,
        (1 - inputs$beta_target) * 100,
        n1, r1, n, r
      )
    }
  }
  
  if (endpoint == "tte") {
    input_scale_null <- if (!is.null(inputs$input_scale_null)) {
      tolower(as.character(inputs$input_scale_null))
    } else {
      "survival"
    }
    use_hazard_scale <- input_scale_null %in% c("hazard", "lambda", "lambda0")
    has_interval <- !is.null(inputs$S0L) && !is.null(inputs$S0U) && inputs$S0L != inputs$S0U
    
    if (has_interval) {
      if (use_hazard_scale) {
        lambda0L <- if (!is.null(inputs$lambda0L_resolved)) inputs$lambda0L_resolved else -log(inputs$S0U) / inputs$x0
        lambda0U <- if (!is.null(inputs$lambda0U_resolved)) inputs$lambda0U_resolved else -log(inputs$S0L) / inputs$x0
        lambda0 <- if (!is.null(inputs$lambda0_resolved)) inputs$lambda0_resolved else -log(inputs$S0) / inputs$x0
        lambda1 <- if (!is.null(inputs$lambda1)) inputs$lambda1 else -log(inputs$S1) / inputs$x0
        
        eval_txt <- sprintf(
          "Operating characteristics are evaluated across the benchmark range [%.2f, %.2f] on the hazard-rate scale, with planning benchmark lambda0 = %.2f and target alternative lambda1 = %.2f. ",
          lambda0L, lambda0U, lambda0, lambda1
        )
      } else {
        eval_txt <- sprintf(
          "Operating characteristics are evaluated across the benchmark range [%.2f, %.2f] on the survival-probability scale, with planning benchmark S0 = %.2f and target alternative S1 = %.2f. ",
          inputs$S0L, inputs$S0U, inputs$S0, inputs$S1
        )
      }
      
      op_txt <- sprintf(
        paste0(
          eval_txt,
          "The probability of early termination ranges from %s. ",
          "The expected sample size ranges from %.2f to %.2f (average %.2f). ",
          "The type I error rate ranges from %s, ",
          "and the achieved power is %.4f."
        ),
        range_txt(row$PET_S0L, row$PET_S0U, digits = 4),
        min(row$EN_S0L, row$EN_S0U), max(row$EN_S0L, row$EN_S0U), row$avg_EN,
        range_txt(row$alpha_at_S0L, row$alpha_at_S0U, digits = 4),
        row$power_at_S1
      )
    } else {
      if (use_hazard_scale) {
        lambda0 <- if (!is.null(inputs$lambda0_resolved)) inputs$lambda0_resolved else -log(inputs$S0) / inputs$x0
        lambda1 <- if (!is.null(inputs$lambda1)) inputs$lambda1 else -log(inputs$S1) / inputs$x0
        
        eval_txt <- sprintf(
          "Operating characteristics are evaluated at planning benchmark lambda0 = %.2f and target alternative lambda1 = %.2f. ",
          lambda0, lambda1
        )
      } else {
        eval_txt <- sprintf(
          "Operating characteristics are evaluated at planning benchmark S0 = %.2f and target alternative S1 = %.2f. ",
          inputs$S0, inputs$S1
        )
      }
      
      op_txt <- sprintf(
        paste0(
          eval_txt,
          "The probability of early termination is %.4f ",
          "and the expected sample size is %.2f. ",
          "The type I error rate is %.4f, ",
          "and the achieved power is %.4f."
        ),
        row$PET_S0,
        row$EN_S0,
        row$alpha_at_S0,
        row$power_at_S1
      )
    }
  } else {
    has_interval <- !is.null(inputs$p0L) && !is.null(inputs$p0U) && inputs$p0L != inputs$p0U
    
    if (has_interval) {
      op_txt <- sprintf(
        paste0(
          "Operating characteristics are evaluated across the benchmark range [%.2f, %.2f], with the planning benchmark p0 = %.2f and target alternative p1 = %.2f. ",
          "The probability of early termination ranges from %s. ",
          "The expected sample size ranges from %.2f to %.2f (average %.2f). ",
          "The type I error rate ranges from %s, ",
          "and the achieved power is %.4f."
        ),
        inputs$p0L, inputs$p0U, inputs$p0, inputs$p1,
        range_txt(row$PET_p0L, row$PET_p0U, digits = 4),
        min(row$EN_p0L, row$EN_p0U), max(row$EN_p0L, row$EN_p0U), row$avg_EN,
        range_txt(row$alpha_at_p0L, row$alpha_at_p0U, digits = 4),
        row$power_at_p1
      )
    } else {
      op_txt <- sprintf(
        paste0(
          "Operating characteristics are evaluated at planning benchmark p0 = %.2f and target alternative p1 = %.2f. ",
          "The probability of early termination is %.4f ",
          "and the expected sample size is %.2f. ",
          "The type I error rate is %.4f, ",
          "and the achieved power is %.4f."
        ),
        inputs$p0, inputs$p1,
        row$PET_p0,
        row$EN_p0,
        row$alpha_at_p0,
        row$power_at_p1
      )
    }
  }
  
  txt <- paste0(design_txt, op_txt)
  structure(txt, class = "opchar_summary")
}


#' Print a Summary of Operating Characteristics
#'
#' @param x Object of class \code{"opchar_summary"}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.opchar_summary <- function(x, ...) {
  cat(unclass(x), "\n")
  invisible(x)
}
