#' Table of Design Inputs
#'
#' Formats the core input parameters for a binary or time-to-event two-stage
#' design object.
#'
#' @param design A design object of class \code{"binary_BUDS"},
#'   or \code{"tte_BUDS"}.
#' @param output Output type: \code{"base"}, \code{"gt"}, or
#'   \code{"flextable"}.
#' @param save_path Optional path used when saving a rendered table.
#'
#' @return A base data frame, \code{gt} table, or \code{flextable}.
#' @export
table_design_inputs <- function(design,
                                output = c("base","gt","flextable"),
                                save_path = NULL) {
  
  output <- match.arg(output)
  if (inherits(design, "binary_BUDS")) {
    inputs <- design$inputs
    endpoint <- "binary"
  } else if (inherits(design, "tte_BUDS")) {
    inputs <- design$inputs
    endpoint <- "tte"
  } else {
    stop("Input must be a binary_BUDS or tte_BUDS object")
  }
  
  # Build data
  out <- data.frame(Parameter = character(), Value = numeric(), stringsAsFactors = FALSE)
  if (endpoint == "binary") {
    ## Interval (only if non-degenerate)
    if (!is.null(inputs$p0L) && !is.null(inputs$p0U)) {
      is_degenerate <- (inputs$p0 == inputs$p0L) && (inputs$p0 == inputs$p0U)
      if (!is_degenerate) {
        out <- rbind(out,
                     data.frame(Parameter = "p0L", Value = inputs$p0L),
                     data.frame(Parameter = "p0U", Value = inputs$p0U))
      }
    }
    ## Core parameters
    out <- rbind(data.frame(Parameter = "p0", Value = inputs$p0),
                 data.frame(Parameter = "p1", Value = inputs$p1),
                 out,
                 data.frame(Parameter = "alpha", Value = inputs$alpha),
                 data.frame(Parameter = "power", Value = inputs$power))
    ## Bounds (optional)
    if (!is.null(inputs$n.ub)) {
      out <- rbind(out,
                   data.frame(Parameter = "n.ub", Value = inputs$n.ub))
    }
  } else {
    use_hazard_scale <- identical(inputs$input_scale_null, "hazard")
    alt_parameter <- if (identical(inputs$input_scale_alt, "lambda1")) {
      "lambda1"
    } else if (identical(inputs$input_scale_alt, "HR")) {
      "HR"
    } else {
      "S1"
    }

    if (use_hazard_scale) {
      lambda0_main <- if (!is.null(inputs$lambda0_resolved)) inputs$lambda0_resolved else inputs$lambda0
      lambda0L_main <- if (!is.null(inputs$lambda0L_resolved)) inputs$lambda0L_resolved else inputs$lambda0L
      lambda0U_main <- if (!is.null(inputs$lambda0U_resolved)) inputs$lambda0U_resolved else inputs$lambda0U
      if (!is.null(lambda0L_main) && !is.null(lambda0U_main) &&
          !(isTRUE(all.equal(lambda0_main, lambda0L_main)) &&
            isTRUE(all.equal(lambda0_main, lambda0U_main)))) {
        out <- rbind(out,
                     data.frame(Parameter = "lambda0L", Value = lambda0L_main),
                     data.frame(Parameter = "lambda0U", Value = lambda0U_main))
      }
      if (!is.null(lambda0_main)) {
        out <- rbind(out, data.frame(Parameter = "lambda0", Value = lambda0_main))
      }
    } else {
      if (!is.null(inputs$S0L) && !is.null(inputs$S0U) &&
          !(isTRUE(all.equal(inputs$S0, inputs$S0L)) &&
            isTRUE(all.equal(inputs$S0, inputs$S0U)))) {
        out <- rbind(out,
                     data.frame(Parameter = "S0L", Value = inputs$S0L),
                     data.frame(Parameter = "S0U", Value = inputs$S0U))
      }
      if (!is.null(inputs$S0)) {
        out <- rbind(out, data.frame(Parameter = "S0", Value = inputs$S0))
      }
    }
    if (alt_parameter == "lambda1" && !is.null(inputs$lambda1)) {
      out <- rbind(out, data.frame(Parameter = "lambda1", Value = inputs$lambda1))
    } else if (alt_parameter == "HR" && !is.null(inputs$HR)) {
      out <- rbind(out, data.frame(Parameter = "HR", Value = inputs$HR))
    } else if (!is.null(inputs$S1)) {
      out <- rbind(out, data.frame(Parameter = "S1", Value = inputs$S1))
    }
    out <- rbind(out,
                 data.frame(Parameter = "x0", Value = inputs$x0),
                 data.frame(Parameter = "rate", Value = inputs$rate),
                 data.frame(Parameter = "alpha", Value = inputs$alpha))
    if (!is.null(inputs$beta)) {
      out <- rbind(out,
                   data.frame(Parameter = "power", Value = 1 - inputs$beta))
    }
    if (!is.null(inputs$n.ub)) {
      out <- rbind(out,
                   data.frame(Parameter = "n.ub", Value = inputs$n.ub))
    }
  }
  out$Value <- round(out$Value, 2)
  
  # ================= BASE VERSION =======================
  if (output == "base") {
    return(out)
  }
  
  # ================= GT VERSION =========================
  if (output == "gt") {
    if (!requireNamespace("gt", quietly = TRUE)) {
      stop("Package 'gt' must be installed")
    }
    library(gt)
    gt_tab <- gt(out) %>%
      fmt_number(columns = "Value", decimals = 2) %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_column_labels()) %>%
      cols_align("left", "Parameter") %>%
      cols_align("center", "Value")
    if (!is.null(save_path)) {
      gtsave(gt_tab, save_path)
    }
    return(gt_tab)
  }
  
  # ================= FLEXTABLE VERSION ==================
  if (output == "flextable") {
    if (!requireNamespace("flextable", quietly = TRUE)) {
      stop("Package 'flextable' must be installed")
    }
    library(flextable)
    ft <- flextable(out)
    ft <- ft %>%
      theme_booktabs() %>%
      bold(part = "header") %>%
      align(j = "Parameter", align = "left") %>%
      align(j = "Value", align = "center") %>%
      autofit()
    if (!is.null(save_path)) {
      flextable::save_as_docx(ft, path = save_path)
    }
    return(ft)
  }
}
