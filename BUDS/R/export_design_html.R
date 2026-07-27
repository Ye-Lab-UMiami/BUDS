#' Export a Design Report to HTML
#'
#' Generates a self-contained HTML report for a binary or TTE two-stage design,
#' including operating-characteristics tables, type I error plots, design
#' rules, and summary text.
#'
#' @param design A design object of class \code{"binary_BUDS"}
#'   or \code{"tte_BUDS"}.
#' @param design_name Name of the focal design to summarize in the report.
#' @param designs_keep Optional character vector indicating which designs should
#'   appear in the operating-characteristics table and type I error rate figure.
#' @param file Output HTML file path.
#'
#' @return Invisibly returns the output file path after writing the HTML file.
#' @export
export_design_html <- function(design,
                               design_name = NULL,
                               designs_keep = NULL,
                               file) {
  
  library(htmltools)
  library(gt)
  if (missing(file) || is.null(file)) {
    stop("You must specify a file path.")
  }
  endpoint <- if (inherits(design, "binary_BUDS")) {
    "binary"
  } else if (inherits(design, "tte_BUDS")) {
    "tte"
  } else {
    stop("design must be a binary_BUDS or tte_BUDS object")
  }
  design_names_all <- colnames(design$designs)
  
  .match_designs <- function(patterns) {
    if (is.null(patterns)) return(design_names_all)
    unique(unlist(lapply(patterns, function(p) {
      p_lc <- tolower(p)
      d_lc <- tolower(design_names_all)
      # exact match
      exact <- design_names_all[d_lc == p_lc]
      if (length(exact) > 0) return(exact)
      # starts-with match
      starts <- design_names_all[grepl(paste0("^", p_lc), d_lc)]
      if (length(starts) > 0) return(starts)
      # fallback contains
      design_names_all[grepl(p_lc, d_lc, fixed = TRUE)]
    })))
  }
  if (!is.null(designs_keep) && any(designs_keep == "BUDS")) {
    BUDS_names <- grep("^BUDS", design_names_all, value = TRUE)
    designs_keep <- unique(c(designs_keep[designs_keep != "BUDS"], BUDS_names))
  }
  if (!is.null(design_name) && design_name == "BUDS") {
    BUDS_names <- grep("^BUDS", design_names_all, value = TRUE)
    design_name <- BUDS_names[1]
  }
  if (is.null(design_name)) {
    design_name <- design_names_all[1]
  } else if (!design_name %in% design_names_all) {
    matched <- .match_designs(design_name)
    if (length(matched) == 0) {
      stop("design_name must match one of: ",
           paste(design_names_all, collapse = ", "))
    }
    design_name <- matched[1]
  }
  if (!is.null(designs_keep)) {
    designs_keep <- .match_designs(designs_keep)
  }
  if (is.null(designs_keep)) {
    designs_keep <- design_names_all
  } else if (!(design_name %in% designs_keep)) {
    designs_keep <- unique(c(designs_keep, design_name))
  }
  
  # ---------------------------
  # Operating characteristics
  # ---------------------------
  if (inherits(design, "binary_BUDS")) {
    op <- opchar_binary_BUDS(design = design, design_name = designs_keep)
  } else if (inherits(design, "tte_BUDS")) {
    op <- opchar_tte_BUDS(design = design, design_name = designs_keep)
  } else {
    stop("design must be a binary_BUDS or tte_BUDS object")
  }
  
  # ---------------------------
  # Summaries
  # ---------------------------
  design_summary <- paste(capture.output(summary(design, design_name)), collapse = "\n")
  
  op_summary <- paste(capture.output(summary(op, design_name)), collapse = "\n")
  
  # ---------------------------
  # Tables
  # ---------------------------
  opchar_html <- HTML(gt::as_raw_html(
    table_design_opchar(design, design_name = designs_keep, output = "gt", details = TRUE)))
  
  rules_html <- HTML(gt::as_raw_html(
    table_design_rules(design, design_name, output = "gt")))
  
  inputs_html <- HTML(gt::as_raw_html(
    table_design_inputs(design, output = "gt")))
  
  # ---------------------------
  # Plot
  # ---------------------------
  plot_file <- tempfile(fileext = ".png")
  plot_obj <- plot_type_I_error(design, designs = designs_keep)
  ggplot2::ggsave(plot_file, plot_obj, width = 6.5, height = 4, dpi = 300)
  plot_uri <- base64enc::dataURI(file = plot_file, mime = "image/png")
  
  # ---------------------------
  # References
  # ---------------------------
  refs <- tags$ol(
    tags$li("Simon, R. (1989). Optimal two-stage designs for phase II clinical trials. Controlled Clinical Trials, 10(1), 1-10. https://doi.org/10.1016/0197-2456(89)90015-9"),
    tags$li("Belin, L., De Rycke, Y., & Broet, P. (2017). A two-stage design for phase II trials with time-to-event endpoint using restricted follow-up. Contemporary Clinical Trials Communications, 8, 127-134. https://doi.org/10.1016/j.conctc.2017.09.010"),
    tags$li("Irlmeier, R., Jin, Z., & Ye, F. BUDS: Benchmark Uncertainty Design Selection for Two-Stage Single-Arm Phase II Trials.")
  )
  
  # ---------------------------
  # Titles
  # ---------------------------
  title <- if (endpoint == "tte") {
    "Two-Stage Designs for Time-to-Event Endpoints"
  } else {
    "Two-Stage Designs for Tests of One Proportion"
  }
  
  subtitle <- if (endpoint == "tte") {
    "restricted-Kwak and Jung and BUDS Criterion (BUDS)"
  } else {
    "Simon and BUDS Criterion (BUDS)"
  }
  
  figure_title <- if (endpoint == "tte") {
    "Figure: Type I Error Rate Across the Benchmark Range"
  } else {
    "Figure: Type I Error Rate Across the Benchmark Range"
  }
  
  # ---------------------------
  # HTML output
  # ---------------------------
  html <- tags$html(
    tags$head(
      tags$title(title),
      tags$style(HTML("
        body {
          font-family: Arial;
          margin: 40px;
          max-width: 1400px;
          line-height: 1.5;
        }
        h1 { font-size: 28px; margin-bottom: 5px; }
        .subtitle {
          font-size: 20px;
          font-weight: 600;
          color: #555;
          margin-top: -5px;
          margin-bottom: 25px;
        }
        h3 { margin-top: 30px; }
        h4 { margin-top: 20px; }
        .gt_table, table {
          margin-left: 0 !important;
          margin-right: auto !important;
        }
      "))
    ),
    tags$body(
      h1(title),
      tags$div(subtitle, class = "subtitle"),
      tags$hr(style = "margin-top: 5px; margin-bottom: 25px;"),
      div(style = "width: 100%; overflow-x: auto;", opchar_html),
      h4(figure_title),
      div(
        style = "text-align: left;",
        tags$img(src = plot_uri,
                 style = "width: 800px; margin-top: 10px;")
      ),
      h3(paste0(design_name, " Design")),
      h4("Design Summary"),
      tags$p(style = "white-space: pre-line;", design_summary),
      h4("Decision Rules"),
      rules_html,
      h4("Operating Characteristics Summary"),
      tags$p(style = "white-space: pre-line;", op_summary),
      h3("References"),
      refs,
      h3("Appendix: Design Inputs"),
      inputs_html
    )
  )
  
  save_html(html, file)
  message("Saved to: ", file)
}
