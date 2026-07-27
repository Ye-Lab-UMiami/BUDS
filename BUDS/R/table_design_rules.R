#' Table of Two-Stage Decision Rules
#'
#' Formats the stage-wise sample sizes and decision rules for a binary or
#' time-to-event two-stage design.
#'
#' @param x A design object of class \code{"binary_BUDS"},
#'   or \code{"tte_BUDS"}.
#' @param design_name Name of the design to extract. Required for multi-design
#'   objects such as \code{"binary_BUDS"} and \code{"tte_BUDS"}.
#' @param output Output type: \code{"base"}, \code{"gt"}, or
#'   \code{"flextable"}.
#' @param save_path Optional path used when saving a rendered table.
#' @param ... Unused.
#'
#' @return A base data frame, \code{gt} table, or \code{flextable}.
#' @export
table_design_rules <- function(x,
                               design_name,
                               output = c("base","gt","flextable"),
                               save_path = NULL,
                               ...) {
  
  # Checks
  output <- match.arg(output)
  if (missing(design_name)) {
    stop("Argument 'design_name' must be provided.")
  }
  if (length(design_name) != 1) {
    stop("Argument 'design_name' must be a single value (length 1).")
  }
  if (!is.character(design_name)) {
    stop("Argument 'design_name' must be a character string.")
  }
  if (inherits(x, "binary_BUDS")) {
    dm <- x$designs
    if (!design_name %in% colnames(dm)) {
      stop(sprintf("Design '%s' not found. Available designs: %s",
                   design_name, paste(colnames(dm), collapse = ", ")))
    }
    design_col <- dm[, design_name]
    df <- data.frame(
      Stage = c("Stage 1", "Stage 2"),
      `Sample Size (n)` = c(design_col["n1"], design_col["n"]),
      Decision = c(sprintf("Stop for futility \u2264 %d responses", design_col["r1"]),
                   sprintf("Declare promising if > %d responses", design_col["r"])),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  } else if (inherits(x, "tte_BUDS")) {
    dm <- .tte_design_matrix_to_df(x$designs)
    if (!design_name %in% dm$design) {
      stop(sprintf("Design '%s' not found. Available designs: %s",
                   design_name, paste(dm$design, collapse = ", ")))
    }
    row <- dm[dm$design == design_name, , drop = FALSE][1, , drop = FALSE]
    df <- data.frame(
      Stage = c("Stage 1", "Stage 2"),
      `Sample Size (n)` = c(row$n1, row$n2),
      Decision = c(sprintf("Stop for futility if Z1 > %.4f", row$c1),
                   sprintf("Declare promising if Z2 \u2264 %.4f", row$c2)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  } else {
    stop("Input must be a binary_BUDS or tte_BUDS object")
  }
  
  # ================= BASE VERSION =======================
  if (output == "base") {
    df_fmt <- df
    df_fmt$Stage <- format(df_fmt$Stage, justify = "left")
    df_fmt$`Sample Size (n)` <- format(df_fmt$`Sample Size (n)`, justify = "right")
    df_fmt$Decision <- format(df_fmt$Decision, justify = "left")
    return(df_fmt) 
  }
  
  # ================= GT VERSION =========================
  if (output == "gt") {
    if (!requireNamespace("gt", quietly = TRUE)) {
      stop("Package 'gt' must be installed")
    }
    library(gt)
    gt_tab <- gt(df) %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_column_labels()) %>%
      cols_align("left", c("Stage", "Decision")) %>%
      cols_align("center", "Sample Size (n)")
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
    ft <- flextable(df)
    ft <- ft %>%
      theme_booktabs() %>%
      bold(part = "header") %>%
      align(j = "Stage", align = "left") %>%
      align(j = "Decision", align = "left") %>%
      align(j = "Sample Size (n)", align = "center") %>%
      autofit()
    if (!is.null(save_path)) {
      flextable::save_as_docx(ft, path = save_path)
    }
    return(ft)
  }
}
