#' Binary Type I Error Rates for Plotting
#'
#' Internal helper used by [plot_type_I_error()] to evaluate exact binary
#' type I error rates over the supplied benchmark range.
#'
#' @param row One-row data frame containing binary design parameters.
#' @param grid_n Number of grid points used for the response-rate grid.
#'
#' @return A tibble with grid values and type I error rates.
#' @keywords internal
.binary_type_I_error <- function(row, grid_n = 101) {
  n1  <- as.integer(row[["n1"]])
  r1  <- as.integer(row[["r1"]])
  n   <- as.integer(row[["n"]])
  r   <- as.integer(row[["r"]])
  p0L <- as.numeric(row[["p0L"]])
  p0U <- as.numeric(row[["p0U"]])
  
  vals <- c(n1, r1, n, r, p0L, p0U)
  if (!all(is.finite(vals))) {
    return(tibble::tibble(Design_all = row$Design_all,
                          p = NA_real_,
                          p0L = p0L,
                          p0U = p0U,
                          reject_prob = NA_real_,
                          skipped = TRUE))
  }
  grid <- seq(p0L, p0U, length.out = grid_n)
  reject_prob <- vapply(grid, function(pp) pr_reject_binary_BUDS(n1, r1, n, r, pp), numeric(1))
  tibble::tibble(Design_all = row$Design_all,
                 p = grid,
                 p0L = p0L,
                 p0U = p0U,
                 reject_prob = reject_prob,
                 skipped = FALSE)
}

#' Time-to-Event Type I Error Rates for Plotting
#'
#' Internal helper used by [plot_type_I_error()] to evaluate approximate
#' type I errors for restricted Kwak-Jung designs over the supplied
#' survival-probability interval.
#'
#' @param row One-row data frame containing TTE design parameters.
#' @param S0L,S0U Null survival-probability bounds.
#' @param x0 Restricted follow-up time.
#' @param rate Accrual rate per time unit.
#' @param grid_n Number of grid points used for the survival grid.
#'
#' @return A tibble with grid values and type I errors.
#' @keywords internal
.tte_type_I_error <- function(row, S0L, S0U, x0, rate, grid_n = 101) {
  vals <- c(row$n1, row$c1, row$n2, row$c2, row$DA1, row$DA2, row$S0_designed, S0L, S0U, x0, rate)
  if (!all(is.finite(vals))) {
    return(tibble::tibble(Design_all = row$Design_all,
                          s = NA_real_,
                          lambda = NA_real_,
                          S0L = S0L,
                          S0U = S0U,
                          reject_prob = NA_real_,
                          skipped = TRUE))
  }
  grid <- seq(S0L, S0U, length.out = grid_n)
  reject_prob <- vapply(grid, function(ss) {
    pr_reject_tte_BUDS(n1 = as.integer(row$n1), c1 = row$c1,
                      n2 = as.integer(row$n2), c2 = row$c2,
                      DA1 = row$DA1, DA2 = row$DA2,
                      S0_design = as.numeric(row$S0_designed),
                      s = ss, x0 = x0, rate = rate)
  }, numeric(1))
  tibble::tibble(Design_all = row$Design_all,
                 s = grid,
                 lambda = -log(grid) / x0,
                 S0L = S0L,
                 S0U = S0U,
                 reject_prob = reject_prob,
                 skipped = FALSE)
}

#' Plot Type I Error Rates Across the Benchmark Range
#'
#' Draws the type I error curve across the benchmark range for binary
#' two-stage designs or restricted Kwak-Jung time-to-event designs.
#'
#' @param design A design object of class \code{"binary_BUDS"}
#'   or \code{"tte_BUDS"}.
#' @param designs Optional character vector used to subset designs to plot.
#' @param grid_n Number of grid points used for the x-axis evaluation grid.
#' @param p0L,p0U Optional benchmark bounds overriding the values stored in a
#'   \code{"binary_BUDS"} design object.
#' @param S0L,S0U Optional benchmark bounds overriding the values stored in an
#'   \code{"tte_BUDS"} design object.
#'
#' @return A \code{ggplot} object.
#' @export
plot_type_I_error <- function(design,
                                designs = NULL,
                                grid_n = 101,
                                p0L = NULL,
                                p0U = NULL,
                                S0L = NULL,
                                S0U = NULL) {
  
  library(dplyr)
  library(purrr)
  library(ggplot2)
  alpha_ref <- NULL
  
  # =============================
  # SIMON (BINARY)
  # =============================
  if (inherits(design, "binary_BUDS")) {
    inputs <- design$inputs
    alpha_ref <- inputs$alpha
    if (!is.null(S0L) || !is.null(S0U)) {
      stop("Use p0L/p0U with binary_BUDS objects.")
    }
    p0L <- if (is.null(p0L)) inputs$p0L else p0L
    p0U <- if (is.null(p0U)) inputs$p0U else p0U
    df <- as.data.frame(t(design$designs))
    df$Design_all <- rownames(df)
    rownames(df) <- NULL
    df$p0L <- p0L
    df$p0U <- p0U
    if (!is.null(designs)) {
      design_names_all <- df$Design_all
      matched <- unique(unlist(lapply(designs, function(p) {
        p_lc <- tolower(p)
        d_lc <- tolower(design_names_all)
        exact <- design_names_all[d_lc == p_lc]
        if (length(exact) > 0) return(exact)
        design_names_all[grepl(p_lc, d_lc, fixed = TRUE)]
      })))
      if (length(matched) == 0) stop("No matching designs found.")
      df <- df %>% filter(Design_all %in% matched)
    }
    reject_df <- df %>%
      group_split(Design_all) %>%
      map_df(.binary_type_I_error, grid_n = grid_n) %>%
      ungroup()
    xlab <- "True Response Rate"
    xvar <- "p"
    # =============================
    # RKJ BUDS
    # =============================
  } else if (inherits(design, "tte_BUDS")) {
    inputs <- design$inputs
    alpha_ref <- inputs$alpha
    if (!is.null(p0L) || !is.null(p0U)) {
      stop("Use S0L/S0U with tte_BUDS objects.")
    }
    S0L <- if (is.null(S0L)) inputs$S0L else S0L
    S0U <- if (is.null(S0U)) inputs$S0U else S0U
    df <- .tte_design_matrix_to_df(design$designs)
    df$Design_all <- df$design
    if (!is.null(designs)) {
      design_names_all <- df$Design_all
      matched <- unique(unlist(lapply(designs, function(p) {
        p_lc <- tolower(p)
        d_lc <- tolower(design_names_all)
        exact <- design_names_all[d_lc == p_lc]
        if (length(exact) > 0) return(exact)
        design_names_all[grepl(p_lc, d_lc, fixed = TRUE)]
      })))
      if (length(matched) == 0) stop("No matching designs found.")
      df <- df %>% filter(Design_all %in% matched)
    }
    reject_df <- df %>%
      group_split(Design_all) %>%
      map_df(.tte_type_I_error,
             S0L = S0L,
             S0U = S0U,
             x0 = inputs$x0,
             rate = inputs$rate,
             grid_n = grid_n) %>%
      ungroup()
    use_hazard_scale <- any(!vapply(inputs[c("lambda0", "lambda0L", "lambda0U")], is.null, logical(1)))
    if (use_hazard_scale) {
      xlab <- "True Hazard Rate"
      xvar <- "lambda"
      x_left <- min(-log(c(S0L, S0U)) / inputs$x0)
      x_right <- max(-log(c(S0L, S0U)) / inputs$x0)
      reject_df <- reject_df %>% arrange(Design_all, .data[[xvar]])
    } else {
      xlab <- "True Survival Probability"
      xvar <- "s"
      x_left <- S0L
      x_right <- S0U
      reject_df <- reject_df %>% arrange(Design_all, .data[[xvar]])
    }
  } else {
    stop("design must be a binary_BUDS or tte_BUDS object")
  }
  
  present <- unique(reject_df$Design_all)
  binary_order <- c("Optimal", "Minimax", "Balanced")
  BUDS_order <- c("BUDS (Worst Regret)", "BUDS (Average EN)", "BUDS (Worst EN)", "BUDS (Average Regret)", "BUDS (Min N)")
  tte_order <- c("r-KJ", BUDS_order)
  if (any(present %in% tte_order)) {
    final_order <- c(tte_order[tte_order %in% present], setdiff(present, tte_order))
  } else {
    final_order <- c(binary_order[binary_order %in% present],
                     BUDS_order[BUDS_order %in% present],
                     setdiff(present, c(binary_order, BUDS_order)))
  }
  reject_df$Design_all <- factor(reject_df$Design_all, levels = final_order)
  if (!exists("x_left", inherits = FALSE) || !exists("x_right", inherits = FALSE)) {
    x_left <- p0L
    x_right <- p0U
  }
  x_limits <- c(x_left, x_right)
  if (isTRUE(all.equal(x_left, x_right))) {
    eps <- max(1e-4, abs(x_left) * 0.01)
    x_limits <- c(x_left - eps, x_right + eps)
  }

  x_scale <- scale_x_continuous(
    limits = x_limits,
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 0.0001)
  )
  
  # ---------------------------
  # Plot
  # ---------------------------
  ggplot(reject_df, aes(x = .data[[xvar]], y = reject_prob, color = Design_all)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = alpha_ref, linetype = "dashed",
               color = "#444444", linewidth = 0.7) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 0.0001)) +
    x_scale +
    labs(x = xlab, y = "Type I Error Rate",
         color = "Design") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")
}
