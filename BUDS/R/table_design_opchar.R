#' Table of Operating Characteristics for Binary or TTE Two-Stage Designs
#'
#' Formats operating characteristics from [opchar_binary_BUDS()] or
#' [opchar_tte_BUDS()] into a compact table for reporting. Supports binary
#' `binary_BUDS` designs and time-to-event `tte_BUDS` objects.
#'
#' @param design A design object of class \code{"binary_BUDS"}
#'   or \code{"tte_BUDS"}.
#' @param design_name Optional character vector used to subset designs before
#'   tabulation.
#' @param save_path Optional file path for saving the rendered table when using
#'   \code{output = "flextable"}.
#' @param details Logical; include explanatory footnotes.
#' @param output Output type: \code{"gt"}, \code{"flextable"}, or
#'   \code{"data.frame"}.
#'
#' @return A \code{gt} table, \code{flextable}, or \code{data.frame}.
#' @export
table_design_opchar <- function(design,
                                design_name = NULL,
                                save_path = NULL,
                                details = TRUE,
                                output = c("gt", "flextable", "data.frame")) {
  
  output <- match.arg(output)
  library(dplyr)
  endpoint <- if (inherits(design, "binary_BUDS")) {
    "binary"
  } else if (inherits(design, "tte_BUDS")) {
    "tte"
  } else {
    stop("Input must be a binary_BUDS or tte_BUDS object")
  }
  
  if (!is.null(design_name)) {
    all_designs <- colnames(design$designs)
    if (any(design_name == "BUDS")) {
      BUDS_names <- grep("^BUDS", all_designs, value = TRUE)
      design_name <- unique(c(design_name[design_name != "BUDS"], BUDS_names))
    }
  }
  
  if (inherits(design, "binary_BUDS")) {
    op <- opchar_binary_BUDS(design = design, design_name = design_name)
  } else if (inherits(design, "tte_BUDS")) {
    op <- opchar_tte_BUDS(design = design, design_name = design_name)
  } else {
    stop("Input must be a binary_BUDS or tte_BUDS object")
  }
  df <- op$results
  inputs <- op$inputs
  present <- unique(df$Design)
  
  if (endpoint == "binary") {
    df <- df %>%
      mutate(Boundaries = sprintf("(%d, %d, %d, %d)", n1, r1, n, r),
             power_target = 1 - inputs$beta_target)
    cols <- c("Design", "Boundaries",
              "EN_p0L", "EN_p0", "EN_p0U", "avg_EN",
              "PET_p0L", "PET_p0", "PET_p0U",
              "alpha_target", "alpha_at_p0L", "alpha_at_p0", "alpha_at_p0U", "sup_alpha",
              "power_target", "power_at_p1")
  } else {
    df <- df %>%
      mutate(Boundaries = sprintf("(%d, %.4f, %.4f, %d, %.4f, %.4f)", n1, c1, DA1, n2, c2, DA2),
             power_target = 1 - inputs$beta_target)
    cols <- c("Design", "S0_designed", "Boundaries",
              "EN_S0L", "EN_S0", "EN_S0U", "avg_EN",
              "PET_S0L", "PET_S0", "PET_S0U",
              "alpha_target", "alpha_at_S0L", "alpha_at_S0", "alpha_at_S0U", "sup_alpha",
              "power_target", "power_at_S1")
  }
  
  df <- df[, cols[cols %in% names(df)], drop = FALSE]
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
  if ("S0_designed" %in% names(df)) {
    df$S0_designed <- round(df$S0_designed, 2)
  }
  if (output == "data.frame") return(df)
  
  if (endpoint == "binary") {
    classic_designs <- c("Optimal", "Minimax", "Balanced")
    present_classic <- intersect(classic_designs, present)
    n_classic <- length(present_classic)
    classic_label <- paste(present_classic, collapse = ", ")
    BUDS_designs <- grep("^BUDS", present, value = TRUE)
    n_BUDS <- length(BUDS_designs)
    pluralize <- function(n) ifelse(n == 1, "design", "designs")
    has_classic <- n_classic > 0
    has_BUDS <- n_BUDS > 0
    foot_txt <- ""
    if (has_classic && has_BUDS) {
      foot_txt <- sprintf(paste0("%s %s planned at p0 = %.2f and p1 = %.2f and included as reference %s for comparison. ",
                                 "BUDS-selected %s use the benchmark range p in [%.2f, %.2f] with p1 = %.2f.",
                                 "Operating characteristics evaluated at p0 = %.2f, p0L = %.2f, p0U = %.2f, and p1 = %.2f."),
                          classic_label, pluralize(n_classic),
                          inputs$p0, inputs$p1,
                          pluralize(n_classic),
                          pluralize(n_BUDS),
                          inputs$p0L, inputs$p0U, inputs$p1,
                          inputs$p0, inputs$p0L, inputs$p0U, inputs$p1)
    } else if (has_classic) {
      foot_txt <- sprintf("%s %s planned at p0 = %.2f and p1 = %.2f. Operating characteristics evaluated at p0 = %.2f, p0L = %.2f, p0U = %.2f, and p1 = %.2f.",
                          classic_label, pluralize(n_classic),
                          inputs$p0, inputs$p1,
                          inputs$p0, inputs$p0L, inputs$p0U, inputs$p1)
    } else if (has_BUDS) {
      foot_txt <- sprintf("BUDS-selected %s use the benchmark range p in [%.2f, %.2f] and p1 = %.2f. Operating characteristics evaluated at p0 = %.2f, p0L = %.2f, p0U = %.2f, and p1 = %.2f.",
                          pluralize(n_BUDS),
                          inputs$p0L, inputs$p0U, inputs$p1,
                          inputs$p0, inputs$p0L, inputs$p0U, inputs$p1)
    }
    defs <- c("Optimal" = "Optimal: Simon's (1989) design that minimizes expected sample size.",
              "Minimax" = "Minimax: Simon's (1989) design minimizes maximum sample size.",
              "Balanced" = "Balanced: Ye & Shyr (2007) design that balances stage sizes.",
              "BUDS (Least Regret)" = "BUDS (Least Regret): minimizes least-case regret.",
              "BUDS (Average EN)" = "BUDS (Average EN): minimizes average EN.",
              "BUDS (Average Regret)" = "BUDS (Average Regret): minimizes average regret.",
              "BUDS (Least EN)" = "BUDS (Least EN): minimizes least-case EN.",
              "BUDS (Min N)" = "BUDS (Min N): minimizes total sample size.")
  } else {
    classic_designs <- c("r-KJ")
    present_classic <- intersect(classic_designs, present)
    n_classic <- length(present_classic)
    classic_label <- paste(present_classic, collapse = ", ")
    BUDS_designs <- grep("^BUDS", present, value = TRUE)
    n_BUDS <- length(BUDS_designs)
    pluralize <- function(n) ifelse(n == 1, "design", "designs")
    has_classic <- n_classic > 0
    has_BUDS <- n_BUDS > 0
    S0L_show <- if (!is.null(inputs$S0L)) inputs$S0L else inputs$S0
    S0U_show <- if (!is.null(inputs$S0U)) inputs$S0U else inputs$S0
    foot_txt <- ""
    if (has_classic && has_BUDS) {
      foot_txt <- sprintf(paste0("%s %s planned at S0 = %.2f and S1 = %.2f at x0 = %.2f and included as reference %s for comparison. ",
                                 "BUDS-selected %s planned for S(x0) in [%.2f, %.2f] with S1 = %.2f. ",
                                 "Operating characteristics evaluated at S0 = %.2f, S0L = %.2f, S0U = %.2f, and S1 = %.2f."),
                          classic_label, pluralize(n_classic),
                          inputs$S0, inputs$S1, inputs$x0,
                          pluralize(n_classic),
                          pluralize(n_BUDS),
                          S0L_show, S0U_show, inputs$S1,
                          inputs$S0, S0L_show, S0U_show, inputs$S1)
    } else if (has_classic) {
      foot_txt <- sprintf("%s %s planned at S0 = %.2f and S1 = %.2f at x0 = %.2f. Operating characteristics evaluated at S0 = %.2f, S0L = %.2f, S0U = %.2f, and S1 = %.2f.",
                          classic_label, pluralize(n_classic),
                          inputs$S0, inputs$S1, inputs$x0,
                          inputs$S0, S0L_show, S0U_show, inputs$S1)
    } else if (has_BUDS) {
      foot_txt <- sprintf("BUDS-selected %s use the benchmark range S(x0) in [%.2f, %.2f] and S1 = %.2f at x0 = %.2f. Operating characteristics evaluated at S0 = %.2f, S0L = %.2f, S0U = %.2f, and S1 = %.2f.",
                          pluralize(n_BUDS),
                          S0L_show, S0U_show, inputs$S1,
                          inputs$S0, S0L_show, S0U_show, inputs$S1, inputs$x0)
    }
    defs <- c("r-KJ" = "r-KJ: restricted-Kwak and Jung (2017) design that minimizes expected sample size.",
              "BUDS (Least Regret)" = "BUDS (Lorst Regret): minimizes least-case regret.",
              "BUDS (Average EN)" = "BUDS (Average EN): minimizes average EN.",
              "BUDS (Average Regret)" = "BUDS (Average Regret): minimizes average regret.",
              "BUDS (Least EN)" = "BUDS (Least EN): minimizes least-case EN.",
              "BUDS (Min N)" = "BUDS (Min N): minimizes total sample size.")
  }
  keep_defs <- defs[names(defs) %in% present]
  
  if (output == "gt") {
    if (!requireNamespace("gt", quietly = TRUE)) {
      stop("Package 'gt' must be installed")
    }
    library(gt)
    gt_tab <- gt(df) %>% fmt_number(where(is.numeric), decimals = 4)
    if ("alpha_target" %in% names(df)) {
      gt_tab <- gt_tab %>% fmt_number(c(alpha_target), decimals = 2)
    }
    if ("power_target" %in% names(df)) {
      gt_tab <- gt_tab %>% fmt_number(c(power_target), decimals = 2)
    }
    if ("S0_designed" %in% names(df)) {
      gt_tab <- gt_tab %>% fmt_number(columns = "S0_designed", decimals = 2)
    }
    if (endpoint == "binary") {
      labels <- list(Boundaries = "(n1, r1, n, r)",
                     EN_p0L = "EN(p0L)", EN_p0 = "EN(p0)", EN_p0U = "EN(p0U)", avg_EN = "Avg. EN",
                     PET_p0L = "PET(p0L)", PET_p0 = "PET(p0)", PET_p0U = "PET(p0U)",
                     alpha_target = "Target",
                     alpha_at_p0L = "\u03B1(p0L)",
                     alpha_at_p0 = "\u03B1(p0)",
                     alpha_at_p0U = "\u03B1(p0U)",
                     sup_alpha = "Sup \u03B1",
                     power_target = "Target",
                     power_at_p1 = "Power(p1)")
    } else {
      labels <- list(S0_designed = "Planned S0",
                     Boundaries = "(n1, c1, DA1, n2, c2, DA2)",
                     EN_S0L = "EN(S0L)", EN_S0 = "EN(S0)", EN_S0U = "EN(S0U)", avg_EN = "Avg. EN",
                     PET_S0L = "PET(S0L)", PET_S0 = "PET(S0)", PET_S0U = "PET(S0U)",
                     alpha_target = "Target",
                     alpha_at_S0L = "\u03B1(S0L)",
                     alpha_at_S0 = "\u03B1(S0)",
                     alpha_at_S0U = "\u03B1(S0U)",
                     sup_alpha = "Sup \u03B1",
                     power_target = "Target",
                     power_at_S1 = "Power(S1)")
    }
    labels <- labels[names(labels) %in% names(df)]
    gt_tab <- do.call(gt::cols_label, c(list(gt_tab), labels))
    left_cols <- intersect(c("Design", "Boundaries"), names(df))
    gt_tab <- gt_tab %>%
      tab_style(cell_text(weight = "bold"), cells_column_labels()) %>%
      cols_align("left", left_cols) %>%
      cols_align("center", setdiff(names(df), left_cols))
    
    if ("Boundaries" %in% names(df)) {
      gt_tab <- gt_tab %>% tab_spanner(label = "Boundaries", columns = "Boundaries", id = "sp_boundaries")
    }
    if (endpoint == "binary") {
      en_cols  <- c("EN_p0L", "EN_p0", "EN_p0U", "avg_EN")
      pet_cols <- c("PET_p0L", "PET_p0", "PET_p0U")
    } else {
      en_cols  <- c("EN_S0L", "EN_S0", "EN_S0U", "avg_EN")
      pet_cols <- c("PET_S0L", "PET_S0", "PET_S0U")
    }
    if (any(en_cols %in% names(df))) {
      gt_tab <- gt_tab %>% tab_spanner(label = "Expected Sample Size (EN)", columns = intersect(c("EN_S0L", "EN_S0", "EN_S0U", "avg_EN"), names(df)), id = "sp_en")
    }
    if (any(pet_cols %in% names(df))) {
      gt_tab <- gt_tab %>% tab_spanner(label = "Probability of Early Termination (PET)", columns = intersect(c("PET_S0L", "PET_S0", "PET_S0U"), names(df)), id = "sp_pet")
    }
    alpha_cols <- if (endpoint == "binary") {
      c("alpha_target", "alpha_at_p0L", "alpha_at_p0", "alpha_at_p0U", "sup_alpha")
    } else {
      c("alpha_target", "alpha_at_S0L", "alpha_at_S0", "alpha_at_S0U", "sup_alpha")
    }
    if (any(alpha_cols %in% names(df))) {
      gt_tab <- gt_tab %>% tab_spanner(label = "Type I Error", columns = intersect(alpha_cols, names(df)), id = "sp_alpha")
    }
    if (any(c("power_target", "power_at_S1", "power_at_p1") %in% names(df))) {
      gt_tab <- gt_tab %>% tab_spanner(label = "Power", columns = intersect(c("power_target", "power_at_S1", "power_at_p1"), names(df)), id = "sp_power")
    }
    if (nzchar(foot_txt)) {
      gt_tab <- gt_tab %>% tab_source_note(gt::md(paste0("**", foot_txt, "**")))
    }
    if (length(keep_defs) > 0) {
      gt_tab <- gt_tab %>% tab_source_note(gt::md(paste("Designs:", paste(keep_defs, collapse = " "))))
    }
    if (details) {
      if (endpoint == "binary") {
        gt_tab <- gt_tab %>%
          tab_source_note(gt::md("n1: sample size in stage 1.")) %>%
          tab_source_note(gt::md("r1: early stopping boundary (stop if responses ≤ r1).")) %>%
          tab_source_note(gt::md("n: total sample size.")) %>%
          tab_source_note(gt::md("r: final rejection boundary (reject if responses > r).")) %>%
          tab_source_note(gt::md("Expected Sample Size (EN): expected sample size under repeated sampling.")) %>%
          tab_source_note(gt::md("EN(p0L), EN(p0), EN(p0U): evaluated at the corresponding response rates.")) %>%
          tab_source_note(gt::md("Avg. EN: average expected sample size across the interval [p0L, p0U].")) %>% 
          tab_source_note(gt::md("Probability of Early Termination (PET): probability of stopping after stage 1.")) %>%
          tab_source_note(gt::md("PET(p0L), PET(p0), PET(p0U): evaluated at the corresponding response rates.")) %>%
          tab_source_note(gt::md("Type I Error: probability of rejecting H0 when true.")) %>%
          tab_source_note(gt::md("Target: required type I error rate.")) %>%
          tab_source_note(gt::md("α(p0L), α(p0), α(p0U): achieved type I error at the corresponding response rates.")) %>%
          tab_source_note(gt::md("Power: probability of rejecting H0 when p = p1.")) %>% 
          tab_source_note(gt::md("Target: required power.")) %>% 
          tab_source_note(gt::md("Power(p1): achieved power at p1."))
      } else {
        gt_tab <- gt_tab %>%
          tab_source_note(gt::md("Planned S0: value of S0 at which the selected design was constructed.")) %>% 
          tab_source_note(gt::md("n1: sample size in stage 1.")) %>%
          tab_source_note(gt::md("c1: stage-1 futility boundary; stop early if Z1 > c1.")) %>%
          tab_source_note(gt::md("n2: additional sample size in stage 2.")) %>%
          tab_source_note(gt::md("c2: final boundary; declare treatment promising if Z2 ≤ c2.")) %>%
          tab_source_note(gt::md("DA1, DA2: analysis time points for stage 1 and stage 2.")) %>%
          tab_source_note(gt::md("Expected Sample Size (EN): expected sample size under repeated sampling.")) %>%
          tab_source_note(gt::md("EN(S0L), EN(S0), EN(S0U): evaluated at the corresponding survival probabilities.")) %>%
          tab_source_note(gt::md("Avg. EN: average expected sample size across the interval [S0L, S0U].")) %>% 
          tab_source_note(gt::md("Probability of Early Termination (PET): probability of stopping after stage 1.")) %>%
          tab_source_note(gt::md("PET(S0L), PET(S0), PET(S0U): evaluated at the corresponding survival probabilities.")) %>%
          tab_source_note(gt::md("Type I Error: probability of rejecting H0 when true.")) %>%
          tab_source_note(gt::md("Target: required type I error rate.")) %>%
          tab_source_note(gt::md("α(S0L), α(S0), α(S0U): achieved type I error at the corresponding survival probabilities.")) %>%
          tab_source_note(gt::md("Power: probability of rejecting H0 when S(x0) = S1.")) %>% 
          tab_source_note(gt::md("Target: required power.")) %>% 
          tab_source_note(gt::md("Power(S1): achieved power at S1."))
      }
    }
    
    return(gt_tab)
  }
  
  if (output == "flextable") {
    if (!requireNamespace("flextable", quietly = TRUE)) {
      stop("Package 'flextable' must be installed")
    }
    if (!requireNamespace("officer", quietly = TRUE)) {
      stop("Package 'officer' must be installed")
    }
    library(flextable)
    ft <- flextable(df) %>%
      theme_booktabs() %>%
      bold(part = "header") %>%
      align(j = intersect(c("Design", "Boundaries"), names(df)), align = "left") %>%
      align(j = setdiff(names(df), c("Design", "Boundaries")), align = "center") %>%
      autofit()
    if (nzchar(foot_txt)) {
      ft <- add_footer_lines(ft, foot_txt)
    }
    if (length(keep_defs) > 0) {
      ft <- add_footer_lines(ft, paste("Designs:", paste(keep_defs, collapse = " ")))
    }
    if (details) {
      if (endpoint == "binary") {
        ft <- add_footer_lines(ft, c("n1: sample size in stage 1.",
                                     "r1: early stopping boundary (stop if responses ≤ r1).",
                                     "n: total sample size.",
                                     "r: final rejection boundary (reject if responses > r).",
                                     "Expected Sample Size (EN): expected sample size under repeated sampling.",
                                     "EN(p0L), EN(p0), EN(p0U): evaluated at the corresponding response rates.",
                                     "Avg. EN: average expected sample size across the interval [p0L, p0U].",
                                     "Probability of Early Termination (PET): probability of stopping after stage 1.",
                                     "PET(p0L), PET(p0), PET(p0U): evaluated at the corresponding response rates.",
                                     "Type I Error: probability of rejecting H0 when true.",
                                     "Target: required type I error rate.",
                                     "α(p0L), α(p0), α(p0U): achieved type I error at the corresponding response rates.",
                                     "Power: probability of rejecting H0 when p = p1.",
                                     "Target: required power.",
                                     "Power(p1): achieved power at p1."))
      } else {
        ft <- add_footer_lines(ft, c("Planned S0: value of S0 at which the selected design was constructed.",
                                     "n1: sample size in stage 1.",
                                     "c1: stage-1 futility boundary; stop early if Z1 > c1.",
                                     "n2: additional sample size in stage 2.",
                                     "c2: final boundary; declare treatment promising if Z2 ≤ c2.",
                                     "DA1, DA2: analysis time points for stage 1 and stage 2.",
                                     "Expected Sample Size (EN): expected sample size under repeated sampling.",
                                     "EN(S0L), EN(S0), EN(S0U): evaluated at the corresponding survival probabilities.",
                                     "Avg. EN: average expected sample size across the interval [S0L, S0U].",
                                     "Probability of Early Termination (PET): probability of stopping after stage 1.",
                                     "PET(S0L), PET(S0), PET(S0U): evaluated at the corresponding survival probabilities.",
                                     "Type I Error: probability of rejecting H0 when true.",
                                     "Target: required type I error rate.",
                                     "α(S0L), α(S0), α(S0U): achieved type I error at the corresponding survival probabilities.",
                                     "Power: probability of rejecting H0 when S(x0) = S1.",
                                     "Target: required power.",
                                     "Power(S1): achieved power at S1."))
      }
    }
    if (!is.null(save_path)) {
      flextable::save_as_docx(ft, path = save_path)
    }
    return(ft)
  }
}
