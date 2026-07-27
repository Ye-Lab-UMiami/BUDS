#' Operating Characteristics for Simon and with BUDS criterion
#'
#' Calculates operating characteristics either from a `binary_BUDS` object or
#' from manually supplied two-stage decision boundaries.
#'
#' @param design Optional object of class `binary_BUDS`.
#' @param n1,r1,n,r Manual binary decision boundaries used when `design` is
#'   `NULL`.
#' @param p0,p1 Null and alternative response rates.
#' @param alpha_target Target one-sided type I error rate.
#' @param p0L,p0U Optional benchmark bounds. When omitted, a fixed benchmark is
#'   used.
#' @param grid_points Number of grid points used to evaluate operating characteristics
#'   under BUDS
#' @param design_name Optional character vector of design names to keep when a
#'   `design` object is supplied.
#'
#' @return An object of class `opchar`.
#' @export
opchar_binary_BUDS <- function(design = NULL,
                             n1 = NULL, r1 = NULL, n = NULL, r = NULL,
                             p0 = NULL, p1 = NULL, alpha_target = NULL,
                             p0L = NULL, p0U = NULL,
                             grid_points = 101, design_name = NULL) {
  
  # =============================
  # CASE 1: design object
  # =============================
  if (!is.null(design)) {
    if (!inherits(design, "binary_BUDS")) {
      stop("design must be a binary_BUDS object")
    }
    design_matrix <- design$designs
    inputs <- design$inputs
    if (is.null(p0))           p0 <- inputs$p0
    if (is.null(p1))           p1 <- inputs$p1
    if (is.null(alpha_target)) alpha_target <- inputs$alpha
    if (is.null(p0L))          p0L <- inputs$p0L
    if (is.null(p0U))          p0U <- inputs$p0U
    all_designs <- colnames(design_matrix)
    if (is.null(design_name)) {
      design_name <- all_designs
    } else {
      design_name <- intersect(design_name, all_designs)
      if (length(design_name) == 0) {
        warning("No valid design names provided; using all designs.")
        design_name <- all_designs
      }
    }
    
    results <- lapply(design_name, function(dname) {
      col <- design_matrix[, dname]
      res <- .opchar_binary(n1 = as.integer(col["n1"]),
                            r1 = as.integer(col["r1"]),
                            n  = as.integer(col["n"]),
                            r  = as.integer(col["r"]),
                            p0 = p0,
                            p1 = p1,
                            alpha_target = alpha_target,
                            p0L = p0L,
                            p0U = p0U,
                            grid_points = grid_points)
      res$Design <- dname
      res$n1 <- as.integer(col["n1"])
      res$r1 <- as.integer(col["r1"])
      res$n  <- as.integer(col["n"])
      res$r  <- as.integer(col["r"])
      res
    })
    out <- dplyr::bind_rows(results)
    out <- .standardize_binary(out)
    for (col in names(out)) {
      if (!is.numeric(out[[col]])) next
      if (col %in% c("alpha_target")) {
        out[[col]] <- round(out[[col]], 2)
      } else if (grepl("^EN_", col) || col == "avg_EN") {
        out[[col]] <- round(out[[col]], 2)
      } else if (col == "power_at_p1") {
        out[[col]] <- signif(out[[col]], 3)
      } else {
        out[[col]] <- round(out[[col]], 4)
      }
    }
    out$Design <- design_name
    out <- out[, c("Design", setdiff(names(out), "Design"))]
    return(structure(
      list(results = out,
           inputs = list(
             p0 = p0,
             p1 = p1,
             p0L = p0L,
             p0U = p0U,
             alpha_target = alpha_target,
             beta_target = inputs$beta,
             source = "design",
             endpoint = "binary"),
           design_inputs = design$inputs),
      class = "opchar"))
  }
  
  # =============================
  # CASE 2: manual input
  # =============================
  if (any(is.null(c(n1, r1, n, r, p0, p1, alpha_target)))) {
    stop("Must provide either a design object OR full inputs")
  }
  res <- .opchar_binary(n1, r1, n, r,
                        p0, p1,
                        alpha_target,
                        p0L, p0U,
                        grid_points)
  res$n1 <- n1
  res$r1 <- r1
  res$n  <- n
  res$r  <- r
  res <- .standardize_binary(res)
  for (col in names(res)) {
    if (!is.numeric(res[[col]])) next
    if (col %in% c("alpha_target")) {
      res[[col]] <- round(res[[col]], 2)
    } else if (grepl("^EN_", col) || col == "avg_EN") {
      res[[col]] <- round(res[[col]], 2)
    } else if (col == "power_at_p1") {
      res[[col]] <- signif(res[[col]], 3)
    } else {
      res[[col]] <- round(res[[col]], 4)
    }
  }
  
  return(structure(
    list(results = res,
         inputs = list(
           p0 = p0,
           p1 = p1,
           p0L = p0L,
           p0U = p0U,
           alpha_target = alpha_target,
           source = "manual",
           endpoint = "binary")),
    class = "opchar"))
}


#' Standardize Binary Operating-Characteristic Columns
#'
#' Internal helper that orders binary single benchmark and uncertainty-aware operating
#' characteristic columns before printing, tabulation, or summarization.
#'
#' @param df Data frame of binary operating characteristics.
#'
#' @return A data frame with columns ordered for downstream display.
#' @keywords internal
.standardize_binary <- function(df) {
  base_cols <- c("n1","r1","n","r",
                 "EN_p0","PET_p0",
                 "alpha_target","alpha_at_p0","power_at_p1")
  interval_cols <- c("EN_p0L","EN_p0U","avg_EN",
                     "PET_p0L","PET_p0U",
                     "alpha_at_p0L","alpha_at_p0U","sup_alpha")
  cols <- base_cols
  if (any(interval_cols %in% names(df))) {
    cols <- c("n1","r1","n","r",
              "EN_p0L","EN_p0","EN_p0U","avg_EN",
              "PET_p0L","PET_p0","PET_p0U",
              "alpha_target",
              "alpha_at_p0L","alpha_at_p0","alpha_at_p0U","sup_alpha",
              "power_at_p1")
    cols <- cols[cols %in% names(df)]
  }
  df <- df[, cols, drop = FALSE]
  return(df)
}

#' Evaluate Binary Operating Characteristics
#'
#' Internal helper that computes exact operating characteristics for a single
#' binary two-stage design at a single benchmark or on a benchmark range.
#'
#' @param n1,r1 Stage-1 sample size and futility cutoff.
#' @param n,r Total sample size and final futility cutoff.
#' @param p0,p1 Null and alternative response probabilities.
#' @param alpha_target Target type I error rate.
#' @param p0L,p0U Optional benchmark bounds.
#' @param grid_points Number of grid points for benchmark range evaluation.
#'
#' @return A one-row data frame of operating characteristics.
#' @keywords internal
.opchar_binary <- function(n1, r1, n, r, p0, p1, alpha_target,
                           p0L = NULL, p0U = NULL, grid_points = 101) {
  has_interval <- !is.null(p0L) && !is.null(p0U) && p0L != p0U
  # Check for non-finite inputs
  base_vals  <- c(n1, r1, n, r, p0, p1, alpha_target)
  extra_vals <- if (has_interval) c(p0L, p0U) else NULL
  all_vals   <- c(base_vals, extra_vals)
  if (!all(is.finite(all_vals))) {
    message("Non-finite inputs detected. Skipping calculation.")
    out <- data.frame(alpha_target = alpha_target,
                      alpha_at_p0  = NA_real_,
                      power_at_p1  = NA_real_,
                      PET_p0       = NA_real_,
                      EN_p0        = NA_real_,
                      skipped      = TRUE)
    if (has_interval) {
      out$alpha_at_p0L <- NA_real_
      out$alpha_at_p0U <- NA_real_
      out$sup_alpha    <- NA_real_
      out$PET_p0L      <- NA_real_
      out$PET_p0U      <- NA_real_
      out$EN_p0L       <- NA_real_
      out$EN_p0U       <- NA_real_
      out$avg_EN       <- NA_real_
    }
    return(out)
  }
  
  # Point null calculations
  alpha_at_p0 <- pr_reject_binary_BUDS(n1, r1, n, r, p0)
  power_at_p1 <- pr_reject_binary_BUDS(n1, r1, n, r, p1)
  PET_p0      <- pet_binary_BUDS(n1, r1, p0)
  EN_p0       <- en_binary_BUDS(n1, r1, n, p0)
  out <- data.frame(alpha_target = alpha_target,
                    alpha_at_p0  = pmin(pmax(alpha_at_p0, 0), 1),
                    power_at_p1  = pmin(pmax(power_at_p1, 0), 1),
                    PET_p0       = PET_p0,
                    EN_p0        = EN_p0)
  
  # Interval null calculations
  if (has_interval) {
    if (p0L > p0U) stop("Require p0L <= p0U")
    alpha_at_p0L <- pr_reject_binary_BUDS(n1, r1, n, r, p0L)
    alpha_at_p0U <- pr_reject_binary_BUDS(n1, r1, n, r, p0U)
    PET_p0L      <- pet_binary_BUDS(n1, r1, p0L)
    PET_p0U      <- pet_binary_BUDS(n1, r1, p0U)
    EN_p0L       <- en_binary_BUDS(n1, r1, n, p0L)
    EN_p0U       <- en_binary_BUDS(n1, r1, n, p0U)
    grid      <- seq(p0L, p0U, length.out = grid_points)
    a_grid    <- vapply(grid, function(pp) pr_reject_binary_BUDS(n1, r1, n, r, pp), numeric(1))
    en_grid   <- vapply(grid, function(pp) en_binary_BUDS(n1, r1, n, pp), numeric(1))
    sup_alpha <- max(a_grid, na.rm = TRUE)
    avg_EN    <- mean(en_grid, na.rm = TRUE)
    out$alpha_at_p0L <- pmin(pmax(alpha_at_p0L, 0), 1)
    out$alpha_at_p0U <- pmin(pmax(alpha_at_p0U, 0), 1)
    out$sup_alpha    <- pmin(pmax(sup_alpha, 0), 1)
    out$PET_p0L      <- PET_p0L
    out$PET_p0U      <- PET_p0U
    out$EN_p0L       <- EN_p0L
    out$EN_p0U       <- EN_p0U
    out$avg_EN       <- avg_EN
  }
  return(out)
}

#' Print Operating Characteristics
#'
#' Prints operating-characteristics results produced by
#' [opchar_binary_BUDS()] or [opchar_tte_BUDS()] for binary or time-to-event
#' two-stage designs.
#'
#' @param x Object of class \code{"opchar"}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.opchar <- function(x, ...) {
  
  inputs <- x$inputs
  endpoint <- if (!is.null(inputs$endpoint)) inputs$endpoint else "binary"
  cat("\nOperating Characteristics evaluated at:\n\n")
  if (endpoint == "tte" && !is.null(inputs$S0) && length(inputs$S0) > 0 && is.finite(inputs$S0)) {
    cat(sprintf(" Null survival at x0 (S0): %.2f\n", inputs$S0))
  } else if (!is.null(inputs$p0) && length(inputs$p0) > 0 && is.finite(inputs$p0)) {
    if (endpoint == "tte") {
      cat(sprintf(" Null survival at x0 (p0): %.2f\n", inputs$p0))
    } else {
      cat(sprintf(" Null response rate (p0): %.2f\n", inputs$p0))
    }
  }
  if (endpoint == "tte" && !is.null(inputs$S0L) && !is.null(inputs$S0U) && inputs$S0L != inputs$S0U) {
    cat(sprintf(" Benchmark range [S0L, S0U]: [%.2f, %.2f]\n",
                inputs$S0L, inputs$S0U))
  } else if (!is.null(inputs$p0L) && !is.null(inputs$p0U) && inputs$p0L != inputs$p0U) {
    if (endpoint == "tte") {
      cat(sprintf(" Benchmark range [p0L, p0U]: [%.2f, %.2f]\n",
                  inputs$p0L, inputs$p0U))
    } else {
      cat(sprintf(" Benchmark range [p0L, p0U]: [%.2f, %.2f]\n",
                  inputs$p0L, inputs$p0U))
    }
  }
  if (endpoint == "tte") {
    cat(sprintf(" Target survival at x0 (S1): %.2f\n", if (!is.null(inputs$S1)) inputs$S1 else inputs$p1))
  } else {
    cat(sprintf(" Target response rate (p1): %.2f\n", inputs$p1))
  }
  
  # Print design inputs if manual entry
  if (!is.null(inputs$n1)) {
    cat("\n Design parameters:\n")
    cat(sprintf("  n1 = %d, r1 = %d\n", inputs$n1, inputs$r1))
    cat(sprintf("  n = %d, r = %d\n", inputs$n, inputs$r))
  } else if (endpoint == "tte" && !is.null(inputs$x0)) {
    cat(sprintf(" x0 = %.2f, accrual rate = %.2f\n", inputs$x0, inputs$rate))
  }
  cat("\n")
  df <- x$results
  
  # Rounding
  alpha_target_col <- intersect("alpha_target", names(df))
  prob_cols <- c(
    "alpha_at_p0", "power_at_p1", "alpha_at_p0L", "alpha_at_p0U",
    "alpha_at_S0", "alpha_at_S0L", "alpha_at_S0U", "sup_alpha",
    "power_at_S1"
  )
  pet_cols <- grep("PET", names(df), value = TRUE)
  prob_cols <- unique(c(prob_cols, pet_cols))
  en_cols <- grep("EN", names(df), value = TRUE)
  ## Keep only existing columns
  prob_cols <- intersect(prob_cols, names(df))
  en_cols   <- intersect(en_cols, names(df))
  
  # Formatting
  df[alpha_target_col] <- lapply(df[alpha_target_col], function(col) {
    formatC(col, format = "f", digits = 2)
  })
  df[prob_cols] <- lapply(df[prob_cols], function(col) {
    formatC(col, format = "f", digits = 4)
  })
  df[en_cols] <- lapply(df[en_cols], function(col) {
    formatC(col, format = "f", digits = 2)
  })
  
  # Print results
  df <- as.data.frame(df)
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}
