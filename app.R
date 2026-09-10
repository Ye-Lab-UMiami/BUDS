library(shiny)
library(bslib)
library(DT)
library(plotly)
library(waiter)
library(ggplot2)

app_pkg_root <- file.path(getwd(), "BUDS")
app_desc <- file.path(app_pkg_root, "DESCRIPTION")
if (file.exists(app_desc)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required when launching the app from a local BUDS source folder.")
  }
  pkgload::load_all(app_pkg_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
}

design_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1B2A4A",
  base_font = font_google("Source Sans Pro"),
  heading_font = font_google("DM Sans"),
  code_font = font_google("JetBrains Mono")
)

css_rules <- "
  .navbar { background-color: #1B2A4A !important; }
  .navbar .navbar-brand, .navbar .nav-link { color: #fff !important; }
  .navbar .nav-link.active { border-bottom: 3px solid #fff; font-weight: 600; }
  .app-sidebar { background: #f7f9fc; border-radius: 14px; padding: 18px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
  .card-shell { background: #fff; border-radius: 12px; padding: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); margin-bottom: 16px; }
  .card-title { font-weight: 700; color: #1B2A4A; margin-bottom: 12px; font-size: 14px; text-transform: uppercase; letter-spacing: 1px; }
  .param-group { background: #eef3f8; padding: 14px; border-radius: 10px; margin-bottom: 14px; }
  .param-group h6 { color: #1B2A4A; font-weight: 700; margin-bottom: 10px; font-size: 11px; text-transform: uppercase; letter-spacing: 1.4px; }
  .hypothesis-block-label { color: #1B2A4A; font-size: 18px; font-weight: 600; margin-bottom: 10px; }
  .input-block-label { color: #1B2A4A; margin-bottom: 10px; min-height: 72px; display: flex; flex-direction: column; justify-content: flex-end; }
  .input-block-label .input-label-main { font-size: 18px; font-weight: 600; line-height: 1.1; }
  .input-block-label .input-label-sub { font-size: 14px; line-height: 1.2; margin-top: 4px; color: #4f5d75; }
  .tte-time-input .input-block-label { min-height: 0; margin-bottom: 8px; justify-content: flex-start; }
  .tte-time-input { margin-bottom: 18px; }
  .hypothesis-vs { color: #1B2A4A; font-size: 18px; font-weight: 600; margin: 8px 0 10px 0; }
  .summary-text { font-size: 14px; line-height: 1.8; color: #333; background: #f8f9fa; padding: 18px; border-radius: 8px; border-left: 4px solid #1B2A4A; }
  .summary-text > .shiny-html-output { margin: 0; padding: 0; min-height: 0 !important; }
  .summary-text .html-fill-item, .summary-text .html-fill-container { min-height: 0 !important; }
  .summary-text p { margin: 0; }
  .result-footnote { font-size: 13px; line-height: 1.7; color: #4f5d75; margin-top: 14px; }
  .run-btn { font-size: 15px; font-weight: 600; letter-spacing: 0.4px; padding: 12px; border-radius: 8px; width: 100%; }
  .timing-badge { display: inline-block; background: #e8f5e9; color: #2e7d32; padding: 4px 12px; border-radius: 20px; font-size: 12px; font-weight: 600; }
  .vbox { background: #fff; border-radius: 10px; padding: 18px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); text-align: center; }
  .vbox-icon { font-size: 24px; margin-bottom: 8px; }
  .vbox-value { font-size: 28px; font-weight: 700; color: #1B2A4A; }
  .vbox-title { font-size: 11px; text-transform: uppercase; letter-spacing: 1px; color: #888; margin-top: 4px; }
  .vbox-primary .vbox-icon { color: #1B2A4A; }
  .vbox-info .vbox-icon { color: #2980b9; }
  .vbox-success .vbox-icon { color: #27ae60; }
  .vbox-secondary .vbox-icon { color: #6c757d; }
  .nav-menu .radio label { display: block; width: 100%; margin-bottom: 8px; padding: 10px 12px; border-radius: 8px; background: #fff; border: 1px solid #dce4ef; font-weight: 600; color: #1B2A4A; }
  .nav-menu .radio input[type='radio'] { margin-right: 8px; }
"

make_vbox <- function(icon_name, value, title, cls) {
  div(
    class = paste("vbox", cls),
    div(class = "vbox-icon", icon(icon_name)),
    div(class = "vbox-value", value),
    div(class = "vbox-title", title)
  )
}

matrix_to_df <- function(design_matrix) {
  df <- as.data.frame(t(design_matrix), check.names = FALSE)
  df$Design <- rownames(df)
  rownames(df) <- NULL
  df[, c("Design", setdiff(names(df), "Design")), drop = FALSE]
}

is_single_point_design <- function(design) {
  if (inherits(design, "binary_BUDS")) {
    return(isTRUE(all.equal(design$inputs$p0L, design$inputs$p0U)))
  }
  if (inherits(design, "tte_BUDS")) {
    return(isTRUE(all.equal(design$inputs$S0L, design$inputs$S0U)))
  }
  FALSE
}

display_design_labels <- function(design) {
  design_names <- colnames(design$designs)
  labels <- if (inherits(design, "binary_BUDS") && is_single_point_design(design)) {
    mapped <- c(
      "Optimal" = "Simon's Optimal",
      "Minimax" = "Simon's Minimax",
      "Balanced" = "Simon's Balanced"
    )
    ifelse(design_names %in% names(mapped), unname(mapped[design_names]), design_names)
  } else if (is_single_point_design(design)) {
    ifelse(grepl("^BUDS \\(", design_names), design_names, paste0("BUDS (", design_names, ")"))
  } else {
    design_names
  }
  stats::setNames(design_names, labels)
}

featured_report_design_choices <- function(design) {
  choices <- display_design_labels(design)
  design_names <- unname(choices)
  keep <- grepl("^BUDS", design_names)
  filtered <- choices[keep]
  if (length(filtered) == 0) choices else filtered
}

apply_design_labels <- function(df, design) {
  if (!("Design" %in% names(df))) return(df)
  label_map <- display_design_labels(design)
  reverse_map <- stats::setNames(names(label_map), unname(label_map))
  df$Design <- unname(reverse_map[df$Design])
  df
}

augment_binary_single_point_results <- function(df, design) {
  if (!inherits(design, "binary_BUDS") || !is_single_point_design(design) || !("Design" %in% names(df))) {
    return(df)
  }
  if (!("Simon's Optimal" %in% df$Design) || any(grepl("^BUDS ", df$Design))) {
    return(df)
  }
  BUDS_row <- df[df$Design == "Simon's Optimal", , drop = FALSE]
  objectives <- if (is.null(design$inputs$app_robust_objective)) {
    design$inputs$robust_objective
  } else {
    design$inputs$app_robust_objective
  }
  BUDS_label <- if ("least_regret" %in% objectives) {
    "BUDS Least Regret"
  } else if ("avg_en" %in% objectives) {
    "BUDS Average EN"
  } else {
    "BUDS"
  }
  BUDS_row$Design <- BUDS_label
  rbind(df, BUDS_row)
}

num_or_null <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return(NULL)
  as.numeric(x)
}

subscript_label <- function(base, sub) {
  paste0(base, "<sub>", sub, "</sub>")
}

term_with_subscript <- function(term) {
  m <- regexec("^([A-Za-z\u03b1-\u03c9\u0391-\u03a9\u03bb]+)([0-9A-Za-z]+)$", term, perl = TRUE)
  parts <- regmatches(term, m)[[1]]
  if (length(parts) == 3) {
    paste0(parts[2], "<sub>", parts[3], "</sub>")
  } else {
    term
  }
}

metric_label <- function(base, term) {
  paste0(base, "<sub>", term_with_subscript(term), "</sub>")
}

bsub <- function(base, sub) {
  tags$b(HTML(metric_label(base, sub)))
}

var_label <- function(base, sub = NULL) {
  if (is.null(sub)) return(base)
  HTML(subscript_label(base, sub))
}

bvar <- function(base, sub = NULL) {
  if (is.null(sub)) return(tags$b(base))
  tags$b(HTML(subscript_label(base, sub)))
}

summary_to_html <- function(txt) {
  if (length(txt) == 0 || is.null(txt)) return(HTML(""))
  x <- as.character(txt)
  replacements <- list(
    "H0" = "H<sub>0</sub>",
    "H1" = "H<sub>1</sub>",
    "x0" = "x<sub>0</sub>",
    "lambda0L" = "&lambda;<sub>0L</sub>",
    "lambda0U" = "&lambda;<sub>0U</sub>",
    "lambda1" = "&lambda;<sub>1</sub>",
    "lambda0" = "&lambda;<sub>0</sub>",
    "S0L" = "S<sub>0L</sub>",
    "S0U" = "S<sub>0U</sub>",
    "S1" = "S<sub>1</sub>",
    "S0" = "S<sub>0</sub>",
    "p0L" = "p<sub>0L</sub>",
    "p0U" = "p<sub>0U</sub>",
    "p1" = "p<sub>1</sub>",
    "p0" = "p<sub>0</sub>",
    "n1" = "n<sub>1</sub>",
    "n2" = "n<sub>2</sub>",
    "c1" = "c<sub>1</sub>",
    "c2" = "c<sub>2</sub>",
    "r1" = "r<sub>1</sub>"
  )
  for (pat in names(replacements)) {
    x <- gsub(paste0("\\b", pat, "\\b"), replacements[[pat]], x, perl = TRUE)
  }
  HTML(x)
}

select_existing <- function(df, cols) {
  df[, intersect(cols, names(df)), drop = FALSE]
}

default_binary_nmax <- function(alpha, power, p0u, p1) {
  base_nmax <- binary_nmax_reference(alpha = alpha, power = power, p0u = p0u, p1 = p1)
  if (is.null(base_nmax)) return(NULL)
  as.integer(base_nmax + 10L)
}

binary_nmax_reference <- function(alpha, power, p0u, p1) {
  vals <- c(alpha, power, p0u, p1)
  if (any(!is.finite(vals))) return(NULL)
  if (alpha <= 0 || alpha >= 1 || power <= 0 || power >= 1) return(NULL)
  if (p0u <= 0 || p0u >= 1 || p1 <= 0 || p1 >= 1) return(NULL)
  if (p1 <= p0u) return(NULL)
  za <- qnorm(1 - alpha)
  zb <- qnorm(power)
  num <- za * sqrt(p0u * (1 - p0u)) + zb * sqrt(p1 * (1 - p1))
  as.integer(ceiling((num * num) / ((p1 - p0u) * (p1 - p0u))))
}

binary_opchar_colnames <- function(df) {
  cols <- names(df)
  labels <- cols
  map <- c(
    "n1" = subscript_label("n", "1"),
    "r1" = subscript_label("r", "1"),
    "EN_p0L" = metric_label("EN", "p0L"),
    "EN_p0" = metric_label("EN", "p0"),
    "EN_p0U" = metric_label("EN", "p0U"),
    "avg_EN" = "Average EN",
    "PET_p0L" = metric_label("PET", "p0L"),
    "PET_p0" = metric_label("PET", "p0"),
    "PET_p0U" = metric_label("PET", "p0U"),
    "sup_alpha" = "max(alpha)",
    "power_at_p1" = metric_label("power", "p1")
  )
  labels[cols %in% names(map)] <- unname(map[cols[cols %in% names(map)]])
  labels
}

tte_opchar_colnames <- function(df) {
  cols <- names(df)
  labels <- cols
  map <- c(
    "n1" = subscript_label("n", "1"),
    "n2" = subscript_label("n", "2"),
    "c1" = subscript_label("c", "1"),
    "c2" = subscript_label("c", "2"),
    "DA1" = "DA1",
    "DA2" = "DA2",
    "EN_S0L" = metric_label("EN", "S0L"),
    "EN_S0" = metric_label("EN", "S0"),
    "EN_S0U" = metric_label("EN", "S0U"),
    "avg_EN" = "Average EN",
    "PET_S0L" = metric_label("PET", "S0L"),
    "PET_S0" = metric_label("PET", "S0"),
    "PET_S0U" = metric_label("PET", "S0U"),
    "sup_alpha" = "max(alpha)",
    "power_at_S1" = metric_label("power", "S1")
  )
  labels[cols %in% names(map)] <- unname(map[cols[cols %in% names(map)]])
  labels
}

binary_design_colnames <- function(df) {
  cols <- names(df)
  labels <- cols
  map <- c(
    "r1" = subscript_label("r", "1"),
    "n1" = subscript_label("n", "1"),
    "r" = "r",
    "n" = "n",
    "power" = "power",
    "ratio" = paste0("ratio (", subscript_label("n", "1"), ":", subscript_label("n", "2"), ")")
  )
  labels[cols %in% names(map)] <- unname(map[cols[cols %in% names(map)]])
  labels
}

tte_design_colnames <- function(df) {
  cols <- names(df)
  labels <- cols
  map <- c(
    "c1" = subscript_label("c", "1"),
    "n1" = subscript_label("n", "1"),
    "c2" = subscript_label("c", "2"),
    "n2" = subscript_label("n", "2"),
    "DA1" = "DA1",
    "DA2" = "DA2",
    "PET" = "PET",
    "EN" = "EN",
    "alpha" = "alpha",
    "power" = "power",
    "ratio" = paste0("ratio (", subscript_label("n", "1"), ":", subscript_label("n", "2"), ")")
  )
  labels[cols %in% names(map)] <- unname(map[cols[cols %in% names(map)]])
  labels
}

default_tte_grid_points <- function(null_scale, x0, s0l = NULL, s0u = NULL,
                                    lambda0l = NULL, lambda0u = NULL,
                                    step = 0.01) {
  if (identical(null_scale, "lambda")) {
    if (is.null(x0) || !is.finite(x0) || x0 <= 0 ||
        is.null(lambda0l) || is.null(lambda0u)) {
      return(21L)
    }
    s0l <- exp(-lambda0u * x0)
    s0u <- exp(-lambda0l * x0)
  }
  if (is.null(s0l) || is.null(s0u) || !is.finite(s0l) || !is.finite(s0u) || s0l > s0u) {
    return(21L)
  }
  if (isTRUE(all.equal(s0l, s0u))) {
    return(1L)
  }
  as.integer(max(2, floor((s0u - s0l) / step) + 1))
}

tte_hypothesis_inputs_ui <- function(scale = c("S", "lambda")) {
  scale <- match.arg(scale)
  if (identical(scale, "S")) {
    tagList(
      div(class = "hypothesis-block-label", HTML("(S<sub>0L</sub>, S<sub>0U</sub>)")),
      fluidRow(
        column(6, numericInput("tte_S0L", NULL, NA, min = 0.01, max = 0.99, step = 0.01)),
        column(6, numericInput("tte_S0U", NULL, NA, min = 0.01, max = 0.99, step = 0.01))
      ),
      div(class = "hypothesis-vs", "versus"),
      selectInput(
        "tte_alt_input", "Alternative input",
        choices = c("S1 (alternative survival prob)" = "S1", "HR (hazard ratio)" = "HR"),
        selected = "S1",
        width = "100%"
      ),
      conditionalPanel(
        condition = "input.tte_alt_input == 'S1'",
        div(class = "hypothesis-block-label", HTML("S<sub>1</sub>")),
        fluidRow(
          column(12, numericInput("tte_S1", NULL, NA, min = 0.01, max = 0.99, step = 0.001))
        )
      ),
      conditionalPanel(
        condition = "input.tte_alt_input == 'HR'",
        div(class = "hypothesis-block-label", HTML("HR (&lambda;<sub>0</sub> / &lambda;<sub>1</sub>)")),
        fluidRow(
          column(12, numericInput("tte_HR", NULL, NA, min = 1.001, step = 0.001))
        )
      )
    )
  } else {
    tagList(
      div(class = "hypothesis-block-label", HTML("(\u03bb<sub>0L</sub>, \u03bb<sub>0U</sub>)")),
      fluidRow(
        column(6, numericInput("tte_lambda0L", NULL, NA, min = 0.001, step = 0.001)),
        column(6, numericInput("tte_lambda0U", NULL, NA, min = 0.001, step = 0.001))
      ),
      div(class = "hypothesis-vs", "versus"),
      div(class = "hypothesis-block-label", HTML("\u03bb<sub>1</sub>")),
      fluidRow(
        column(12, numericInput("tte_lambda1", NULL, NA, min = 0.001, step = 0.001))
      )
    )
  }
}

section_menu <- function(id) {
  div(
    class = "param-group nav-menu",
    h6("Results Navigation"),
    radioButtons(
      id, NULL,
      choiceNames = list(
        tagList(icon("table-list"), " Search Results"),
        tagList(icon("chart-line"), " Operating Characteristics"),
        tagList(icon("file-lines"), " Design Details"),
        tagList(icon("download"), " Report")
      ),
      choiceValues = c("results", "opchar", "details", "report"),
      selected = "results"
    )
  )
}

build_results_tab <- function(prefix, title_input, inputs_ui) {
  tabPanel(
    title_input,
    icon = if (prefix == "binary") icon("flask") else icon("clock"),
    br(),
    fluidRow(
      column(
        3,
        div(
          class = "app-sidebar",
          inputs_ui,
          actionButton(paste0("run_", prefix), "Run",
                       class = "btn-primary run-btn", icon = icon("play")),
          br(), br(),
          section_menu(paste0(prefix, "_section"))
        )
      ),
      column(
        9,
        tabsetPanel(
          id = paste0(prefix, "_tabs"),
          type = "hidden",
          tabPanel(
            "results",
            div(class = "card-shell",
                div(class = "card-title", "Selected Designs"),
                DTOutput(paste0(prefix, "_design_table")),
                uiOutput(paste0(prefix, "_results_footnote")))
          ),
          tabPanel(
            "opchar",
            div(class = "card-shell",
                div(class = "card-title", "Operating Characteristics"),
                DTOutput(paste0(prefix, "_opchar_table")),
                uiOutput(paste0(prefix, "_opchar_footnote"))),
            div(class = "card-shell",
                div(class = "card-title", "Type I Error Rate Across the Benchmark Range"),
                div(style = "max-width: 860px; margin: 0 auto;",
                    plotlyOutput(paste0(prefix, "_type_I_error_plotly"), height = "560px")))
          ),
          tabPanel(
            "details",
            div(class = "card-shell",
                div(class = "card-title", "Choose a Design"),
                selectInput(paste0(prefix, "_detail_design"), "Design", choices = NULL, width = "100%")),
            fluidRow(
              column(6,
                     div(class = "card-shell",
                         div(class = "card-title", "Design Summary"),
                        div(class = "summary-text", uiOutput(paste0(prefix, "_summary"))))),
              column(6,
                     div(class = "card-shell",
                         div(class = "card-title", "Decision Rules"),
                         tableOutput(paste0(prefix, "_rules"))))
            ),
            div(class = "card-shell",
                div(class = "card-title", "Operating Characteristics Summary"),
                div(class = "summary-text", uiOutput(paste0(prefix, "_opchar_summary"))))
          ),
          tabPanel(
            "report",
            div(class = "card-shell",
                div(class = "card-title", "Export Local HTML Report"),
                selectInput(paste0(prefix, "_report_design"), "Featured design", choices = NULL, width = "100%"),
                checkboxGroupInput(paste0(prefix, "_report_keep"), "Include designs", choices = NULL),
                downloadButton(paste0(prefix, "_report_download"), "Download HTML Report", class = "btn-primary"),
                br(), br(),
                div(class = "summary-text",
                    "This saves a self-contained HTML report to your local machine."))
          )
        )
      )
    )
  )
}

ui <- navbarPage(
  title = span("Benchmark Uncertainty Design Selection (BUDS)", style = "font-weight:700; letter-spacing:0.5px;"),
  theme = design_theme,
  windowTitle = "BUDS",
  header = tagList(
    useWaiter(),
    tags$head(
      tags$title("BUDS"),
      tags$style(HTML(css_rules)),
      tags$script(HTML("
        document.addEventListener('DOMContentLoaded', function() {
          function setPlaceholders() {
            const placeholders = {
              binary_p0L: 'p0L',
              binary_p0U: 'p0U',
              binary_p1: 'p1',
              tte_S0L: 'S0L',
              tte_S0U: 'S0U',
              tte_S1: 'S1',
              tte_HR: 'HR',
              tte_lambda0L: '\u03bb0L',
              tte_lambda0U: '\u03bb0U',
              tte_lambda1: '\u03bb1'
            };
            Object.entries(placeholders).forEach(([id, value]) => {
              const el = document.getElementById(id);
              if (el) el.setAttribute('placeholder', value);
            });
          }
          setPlaceholders();
          document.addEventListener('shiny:connected', setPlaceholders);
          document.addEventListener('shiny:value', setPlaceholders, true);
        });
      "))
    )
  ),
  build_results_tab(
    "binary",
    "Binary",
    tagList(
      div(
        class = "param-group",
        h6("Hypothesis"),
        div(class = "hypothesis-block-label", HTML("(p<sub>0L</sub>, p<sub>0U</sub>)")),
        fluidRow(
          column(6, numericInput("binary_p0L", NULL, NA, min = 0, max = 0.99, step = 0.01)),
          column(6, numericInput("binary_p0U", NULL, NA, min = 0, max = 0.99, step = 0.01))
        ),
        div(class = "hypothesis-vs", "versus"),
        div(class = "hypothesis-block-label", HTML("p<sub>1</sub>")),
        fluidRow(
          column(12, numericInput("binary_p1", NULL, NA, min = 0.01, max = 0.99, step = 0.01))
        )
      ),
      div(
        class = "param-group",
        h6("Error Rate Constraints"),
        fluidRow(
          column(6, numericInput("binary_alpha", "Alpha", 0.05, min = 0.001, max = 0.25, step = 0.005)),
          column(6, numericInput("binary_power", "Power", 0.80, min = 0.50, max = 0.99, step = 0.01))
        )
      ),
      div(
        class = "param-group",
        h6("Search Space"),
        numericInput("binary_nub", HTML("n<sub>max</sub> (the upper bound on the total sample size)"), 150, min = 10, max = 1000, step = 10)
      ),
      div(
        class = "param-group",
        h6("BUDS Objectives"),
        checkboxGroupInput(
          "binary_objective", NULL,
          choices = c(
            "Least Regret" = "least_regret",
            "Average EN" = "avg_en"
          ),
          selected = c("least_regret", "avg_en")
        )
      )
    )
  ),
  build_results_tab(
    "tte",
    "TTE",
    tagList(
      div(
        class = "param-group",
        h6("Hypothesis Estimand"),
        selectInput("tte_null_scale", NULL,
                    choices = c("Survival probability" = "S",
                                "Hazard rate" = "lambda"),
                    selected = "S", width = "100%"),
        uiOutput("tte_hypothesis_inputs")
      ),
      div(
        class = "param-group",
        h6("Time-To-Event"),
        div(
          class = "tte-time-input",
          div(
            class = "input-block-label",
            div(class = "input-label-main", HTML("x<sub>0</sub>")),
            div(class = "input-label-sub", "restricted follow-up time")
          ),
          numericInput("tte_x0", NULL, 1, min = 0.1, step = 0.1)
        ),
        div(
          class = "tte-time-input",
          div(
            class = "input-block-label",
            div(class = "input-label-main", "rate"),
            div(class = "input-label-sub", "accrual rate per time unit")
          ),
          numericInput("tte_rate", NULL, 15, min = 0.01, step = 0.1)
        )
      ),
      div(
        class = "param-group",
        h6("Error Rate Constraints"),
        fluidRow(
          column(6, numericInput("tte_alpha", "Alpha", 0.05, min = 0.001, max = 0.25, step = 0.005)),
          column(6, numericInput("tte_power", "Power", 0.90, min = 0.50, max = 0.99, step = 0.01))
        )
      ),
      div(
        class = "param-group",
        h6("Search Space"),
        numericInput("tte_nub", HTML("n<sub>max</sub> (the upper bound on the total sample size)"), 100, min = 2, step = 1)
      ),
      div(
        class = "param-group",
        h6("BUDS Objectives"),
        checkboxGroupInput(
          "tte_objective", NULL,
          choices = c(
            "Average EN" = "avg_en",
            "Least Regret" = "least_regret"
          ),
          selected = c("avg_en", "least_regret")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  binary_waiter <- Waiter$new(
    html = tagList(spin_fading_circles(), h4("Searching Binary designs...", style = "color:white; margin-top:20px;")),
    color = "rgba(27,42,74,0.85)"
  )
  tte_waiter <- Waiter$new(
    html = tagList(spin_fading_circles(), h4("Searching TTE designs...", style = "color:white; margin-top:20px;")),
    color = "rgba(27,42,74,0.85)"
  )

  binary_result <- reactiveVal(NULL)
  tte_result <- reactiveVal(NULL)

  observeEvent(input$binary_section, {
    updateTabsetPanel(session, "binary_tabs", selected = input$binary_section)
  }, ignoreInit = TRUE)

  observeEvent(input$tte_section, {
    updateTabsetPanel(session, "tte_tabs", selected = input$tte_section)
  }, ignoreInit = TRUE)

  observeEvent(input$tte_null_scale, {
    tte_result(NULL)
    removeNotification("tte_status")
  }, ignoreInit = TRUE)

  observeEvent(input$tte_alt_input, {
    tte_result(NULL)
    removeNotification("tte_status")
  }, ignoreInit = TRUE)

  observeEvent(
    list(input$binary_p0L, input$binary_p0U, input$binary_p1, input$binary_alpha, input$binary_power, input$binary_objective),
    {
      p0L <- num_or_null(input$binary_p0L)
      p0U <- num_or_null(input$binary_p0U)
      p1 <- num_or_null(input$binary_p1)
      alpha <- num_or_null(input$binary_alpha)
      power <- num_or_null(input$binary_power)
      if (is.null(p0L) || is.null(p0U) || is.null(p1) || is.null(alpha) || is.null(power)) return()
      if (p0L > p0U) return()
      suggested_nmax <- default_binary_nmax(alpha = alpha, power = power, p0u = p0U, p1 = p1)
      if (is.null(suggested_nmax)) return()
      updateNumericInput(session, "binary_nub", value = suggested_nmax)
    },
    ignoreInit = TRUE
  )

  output$tte_hypothesis_inputs <- renderUI({
    null_scale <- if (is.null(input$tte_null_scale)) "S" else input$tte_null_scale
    tte_hypothesis_inputs_ui(null_scale)
  })

  observeEvent(input$run_binary, {
    req(length(input$binary_objective) > 0)
    p0L <- num_or_null(input$binary_p0L)
    p0U <- num_or_null(input$binary_p0U)
    p1 <- num_or_null(input$binary_p1)
    if (is.null(p0L) || is.null(p0U) || is.null(p1)) {
      showNotification("Enter p0L, p0U, and p1 before running the Binary search.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    if (p0L > p0U) {
      showNotification("Require p0L <= p0U for the Binary hypothesis.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    p0 <- (p0L + p0U) / 2
    binary_waiter$show()
    tryCatch({
      d <- BUDS::binary_BUDS(
        alpha = input$binary_alpha,
        power = input$binary_power,
        p0 = p0,
        p1 = p1,
        p0L = p0L,
        p0U = p0U,
        n.ub = as.integer(input$binary_nub),
        robust_objective = input$binary_objective
      )
      d$inputs$app_robust_objective <- input$binary_objective
      binary_result(d)
      design_choices <- display_design_labels(d)
      report_choices <- featured_report_design_choices(d)
      BUDS_cols <- grep("^BUDS", colnames(d$designs), value = TRUE)
      updateSelectInput(session, "binary_detail_design", choices = design_choices, selected = colnames(d$designs)[1])
      updateSelectInput(session, "binary_report_design", choices = report_choices, selected = unname(report_choices)[1])
      updateCheckboxGroupInput(session, "binary_report_keep", choices = design_choices, selected = colnames(d$designs))
      updateRadioButtons(session, "binary_section", selected = "results")
      binary_waiter$hide()
      if (p0L < p0U && length(BUDS_cols) == 0) {
        suggested_nmax <- binary_nmax_reference(alpha = input$binary_alpha, power = input$binary_power, p0u = p0U, p1 = p1)
        showNotification(
          HTML(
            paste0(
              "No BUDS framework found. Hint: increase the gap between [p<sub>0L</sub>, p<sub>0U</sub>] and p<sub>1</sub>",
              if (!is.null(suggested_nmax)) {
                paste0(", or try starting with n<sub>max</sub> = ", suggested_nmax, ".")
              } else {
                ", or use a larger n<sub>max</sub>."
              }
            )
          ),
          type = "error",
          duration = 8
        )
      }
      showNotification(paste("Binary search done in", d$meta$execution_time), type = "message", duration = 4)
    }, error = function(e) {
      binary_waiter$hide()
      showNotification(paste("Binary error:", e$message), type = "error", duration = 8)
    })
  })

  observeEvent(input$run_tte, {
    req(length(input$tte_objective) > 0)
    removeNotification("tte_status")
    null_scale <- if (is.null(input$tte_null_scale)) "S" else input$tte_null_scale
    alt_input <- if (identical(null_scale, "S") && !is.null(input$tte_alt_input)) input$tte_alt_input else "S1"
    x0_val <- num_or_null(input$tte_x0)
    s0l_val <- if (identical(null_scale, "S")) num_or_null(input$tte_S0L) else NULL
    s0u_val <- if (identical(null_scale, "S")) num_or_null(input$tte_S0U) else NULL
    lambda0l_val <- if (identical(null_scale, "lambda")) num_or_null(input$tte_lambda0L) else NULL
    lambda0u_val <- if (identical(null_scale, "lambda")) num_or_null(input$tte_lambda0U) else NULL
    s1_val <- if (identical(null_scale, "S") && identical(alt_input, "S1")) num_or_null(input$tte_S1) else NULL
    hr_val <- if (identical(null_scale, "S") && identical(alt_input, "HR")) num_or_null(input$tte_HR) else NULL
    lambda1_val <- if (identical(null_scale, "lambda")) num_or_null(input$tte_lambda1) else NULL
    if (identical(null_scale, "S") && (is.null(s0l_val) || is.null(s0u_val) || s0l_val > s0u_val)) {
      showNotification("Require S0L and S0U with S0L <= S0U for the TTE hypothesis.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    if (identical(null_scale, "S") && identical(alt_input, "S1") && is.null(s1_val)) {
      showNotification("Enter S1 before running the TTE search.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    if (identical(null_scale, "S") && identical(alt_input, "HR") && (is.null(hr_val) || hr_val <= 1)) {
      showNotification("Enter HR > 1 before running the TTE search.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    if (identical(null_scale, "lambda") && (is.null(lambda0l_val) || is.null(lambda0u_val) || lambda0l_val > lambda0u_val)) {
      showNotification("Require lambda0L and lambda0U with lambda0L <= lambda0U for the TTE hypothesis.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    if (identical(null_scale, "lambda") && is.null(lambda1_val)) {
      showNotification("Enter lambda1 before running the TTE search.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    s0_val <- if (identical(null_scale, "S")) (s0l_val + s0u_val) / 2 else NULL
    lambda0_val <- if (identical(null_scale, "lambda")) (lambda0l_val + lambda0u_val) / 2 else NULL
    grid_points_val <- default_tte_grid_points(
      null_scale = null_scale,
      x0 = x0_val,
      s0l = s0l_val,
      s0u = s0u_val,
      lambda0l = lambda0l_val,
      lambda0u = lambda0u_val,
      step = 0.01
    )
    if (identical(null_scale, "lambda")) {
      if (is.null(x0_val) || !is.finite(x0_val) || x0_val <= 0) {
        showNotification("x0 must be positive when using hazard-rate inputs.", type = "error", duration = 8)
        return(invisible(NULL))
      }
    }
    tte_waiter$show()
    tryCatch({
      d <- withCallingHandlers(
        BUDS::tte_BUDS(
          S0 = s0_val,
          S0L = s0l_val,
          S0U = s0u_val,
          lambda0 = lambda0_val,
          lambda0L = lambda0l_val,
          lambda0U = lambda0u_val,
          S1 = s1_val,
          HR = hr_val,
          lambda1 = lambda1_val,
          x0 = x0_val,
          rate = num_or_null(input$tte_rate),
          alpha = input$tte_alpha,
          beta = 1 - input$tte_power,
          n.ub = as.integer(input$tte_nub),
          n_design_points = grid_points_val,
          n_eval_grid = grid_points_val,
          BUDS_objective = input$tte_objective,
          verbose = FALSE
        ),
        warning = function(w) invokeRestart("muffleWarning")
      )
      tte_result(d)
      design_choices <- display_design_labels(d)
      report_choices <- featured_report_design_choices(d)
      updateSelectInput(session, "tte_detail_design", choices = design_choices, selected = colnames(d$designs)[1])
      updateSelectInput(session, "tte_report_design", choices = report_choices, selected = unname(report_choices)[1])
      updateCheckboxGroupInput(session, "tte_report_keep", choices = design_choices, selected = colnames(d$designs))
      updateRadioButtons(session, "tte_section", selected = "results")
      tte_waiter$hide()
      showNotification(paste("TTE search done in", d$meta$execution_time), type = "message", duration = 4, id = "tte_status")
    }, error = function(e) {
      tte_waiter$hide()
      tte_result(NULL)
      if (identical(e$message, "No feasible designs found within n.ub.") ||
          identical(e$message, "No BUDS-feasible designs. Try narrower [S0L,S0U] or larger n.ub.")) {
        hint_html <- if (identical(null_scale, "lambda")) {
          paste0(
            "TTE error: No feasible designs found within n<sub>max</sub>. ",
            "Hint: narrow [\u03bb<sub>0L</sub>, \u03bb<sub>0U</sub>], move it farther from \u03bb<sub>1</sub>, ",
            "or use a larger n<sub>max</sub>."
          )
        } else {
          paste0(
            "TTE error: No feasible designs found within n<sub>max</sub>. ",
            "Hint: narrow [S<sub>0L</sub>, S<sub>0U</sub>], move it farther from S<sub>1</sub>, ",
            "or use a larger n<sub>max</sub>."
          )
        }
        showNotification(HTML(hint_html), type = "error", duration = 8, id = "tte_status")
      } else {
        showNotification(paste("TTE error:", e$message), type = "error", duration = 8, id = "tte_status")
      }
    })
  })

  binary_opchar <- reactive({
    d <- binary_result()
    req(d)
    BUDS::opchar_binary_BUDS(design = d)
  })

  tte_opchar <- reactive({
    d <- tte_result()
    req(d)
    BUDS::opchar_tte_BUDS(
      design = d,
      S0 = d$inputs$S0,
      S0L = d$inputs$S0L,
      S0U = d$inputs$S0U
    )
  })

  output$binary_value_boxes <- renderUI({
    d <- binary_result()
    if (is.null(d)) {
      return(div(style = "text-align:center; padding:60px; color:#999;",
                 icon("flask", style = "font-size:48px; opacity:0.3;"), br(), br(),
                 "Configure parameters and click ", strong("Run Binary Search")))
    }
    dm <- d$designs
    opt_n  <- if ("Optimal" %in% colnames(dm)) as.integer(dm["n", "Optimal"]) else "--"
    mm_n   <- if ("Minimax" %in% colnames(dm)) as.integer(dm["n", "Minimax"]) else "--"
    opt_en <- if ("Optimal" %in% colnames(dm)) round(dm["EN", "Optimal"], 1) else "--"
    BUDS_cols <- grep("^BUDS", colnames(dm), value = TRUE)
    BUDS_n <- if (length(BUDS_cols) > 0) as.integer(dm["n", BUDS_cols[1]]) else "--"
    fluidRow(
      style = "margin-bottom: 16px;",
      column(3, make_vbox("bullseye", opt_n, "Optimal N", "vbox-primary")),
      column(3, make_vbox("compress", mm_n, "Minimax N", "vbox-info")),
      column(3, make_vbox("chart-line", opt_en, "Optimal EN(midpoint)", "vbox-success")),
      column(3, make_vbox("shield-halved", BUDS_n, "BUDS N", "vbox-secondary"))
    )
  })

  output$binary_timing_badge <- renderUI({
    d <- binary_result()
    req(d)
    tags$span(class = "timing-badge", icon("bolt"), d$meta$execution_time)
  })

  output$binary_design_table <- renderDT({
    d <- binary_result()
    req(d)
    df <- matrix_to_df(d$designs)
    df <- apply_design_labels(df, d)
    df <- augment_binary_single_point_results(df, d)
    op_df <- binary_opchar()$results
    op_df <- apply_design_labels(op_df, d)
    op_df <- augment_binary_single_point_results(op_df, d)
    if ("sup_alpha" %in% names(op_df) && "alpha" %in% names(df)) {
      df$alpha <- unname(op_df$sup_alpha[match(df$Design, op_df$Design)])
    }
    for (col in c("r1", "n1", "r", "n")) if (col %in% names(df)) df[[col]] <- as.integer(df[[col]])
    if ("beta" %in% names(df)) {
      df$power <- 1 - df$beta
    }
    if (all(c("n1", "n") %in% names(df))) {
      stage2_n <- df$n - df$n1
      df$ratio <- ifelse(stage2_n <= 0, NA_real_, df$n1 / stage2_n)
    }
    drop_cols <- intersect(c("alpha1", "alpha2", "beta1", "beta2", "beta"), names(df))
    if (length(drop_cols) > 0) {
      df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
    }
    df <- select_existing(df, c("Design", "r1", "n1", "PET", "r", "n", "EN", "alpha", "power", "ratio"))
    df$.alpha_bad <- ""
    df$.power_bad <- ""
    if (all(c("alpha", "power") %in% names(df))) {
      df$.alpha_bad[round(df$alpha, 2) > round(input$binary_alpha, 2)] <- "violate"
      df$.power_bad[round(df$power, 2) < round(input$binary_power, 2)] <- "violate"
    }
    hidden_targets <- which(names(df) %in% c(".alpha_bad", ".power_bad")) - 1L
    DT::datatable(
      df,
      colnames = binary_design_colnames(df),
      escape = FALSE,
      options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE,
                     columnDefs = list(
                       list(className = "dt-center", targets = "_all"),
                       list(visible = FALSE, targets = hidden_targets)
                     )),
      rownames = FALSE,
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("alpha", valueColumns = ".alpha_bad",
                      backgroundColor = DT::styleEqual("violate", "#fdeaea"),
                      color = DT::styleEqual("violate", "#b00020")) |>
      DT::formatStyle("power", valueColumns = ".power_bad",
                      backgroundColor = DT::styleEqual("violate", "#fdeaea"),
                      color = DT::styleEqual("violate", "#b00020")) |>
      DT::formatRound(intersect(c("PET", "alpha", "power", "ratio"), names(df)), 2) |>
      DT::formatRound(intersect("EN", names(df)), 2)
  })

  output$binary_results_footnote <- renderUI({
    d <- binary_result()
    req(d)
    midpoint_txt <- if (isTRUE(all.equal(d$inputs$p0L, d$inputs$p0U))) {
      objectives <- if (is.null(d$inputs$app_robust_objective)) {
        d$inputs$robust_objective
      } else {
        d$inputs$app_robust_objective
      }
      BUDS_label <- if ("least_regret" %in% objectives) {
        "BUDS Least Regret"
      } else if ("avg_en" %in% objectives) {
        "BUDS Average EN"
      } else {
        "BUDS"
      }
      tagList(
        "Single-point benchmark input is used with ", bvar("p", "0"), " = ", sprintf("%.2f", d$inputs$p0), ". ",
        "In this case, the ", BUDS_label, " design coincides with Simon's Optimal design."
      )
    } else {
      tagList(
        "When the benchmark range is provided, the classical Simon comparator uses the midpoint ",
        bvar("p", "0"), " = (", bvar("p", "0L"), " + ", bvar("p", "0U"), ") / 2 = ",
        sprintf("%.2f", d$inputs$p0), "."
      )
    }
    div(
      class = "result-footnote",
      div(tags$b("n"), " is the total number of subjects required."),
      div(bsub("n", "1"), " is the number of subjects to accrue during stage 1."),
      div(bsub("r", "1"), ", if ", bsub("r", "1"), " or fewer responses are observed during stage 1, the trial is stopped early for futility."),
      div(tags$b("r"), ", if ", tags$b("r"), " or fewer responses are observed by the end of stage two, then no further investigation of the drug is warranted."),
      div(tags$b("alpha"), " is the maximum type I error over the benchmark range."),
      div(tags$b("power"), " is shown as 1 - beta."),
      div(tags$b("EN"), " is the expected sample size for the trial when response rate is ", bvar("p", "0"), "."),
      div(tags$b("PET"), " is the probability of early termination."),
      div("* alpha value highlighted exceeds type I error constraint."),
      div("* ", midpoint_txt)
    )
  })

  output$binary_opchar_table <- renderDT({
    d <- binary_result()
    req(d)
    df <- binary_opchar()$results
    df <- apply_design_labels(df, d)
    df <- augment_binary_single_point_results(df, d)
    if ("alpha_target" %in% names(df)) df$alpha_target <- NULL
    if (!("sup_alpha" %in% names(df)) && "alpha_at_p0" %in% names(df)) df$sup_alpha <- df$alpha_at_p0
    drop_alpha_cols <- intersect(c("alpha_at_p0L", "alpha_at_p0", "alpha_at_p0U"), names(df))
    if (length(drop_alpha_cols) > 0) df <- df[, setdiff(names(df), drop_alpha_cols), drop = FALSE]
    df <- select_existing(df, c("Design", "n1", "r1", "n", "r",
                                "EN_p0L", "EN_p0", "EN_p0U", "avg_EN",
                                "PET_p0L", "PET_p0", "PET_p0U",
                                "sup_alpha", "power_at_p1"))
    df$.alpha_bad <- ""
    if ("sup_alpha" %in% names(df)) {
      df$.alpha_bad[df$sup_alpha > input$binary_alpha] <- "violate"
    }
    hidden_targets <- which(names(df) %in% c(".alpha_bad")) - 1L
    DT::datatable(
      df,
      colnames = binary_opchar_colnames(df),
      escape = FALSE,
      options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE,
                     columnDefs = list(
                       list(className = "dt-center", targets = "_all"),
                       list(visible = FALSE, targets = hidden_targets)
                     )),
      rownames = FALSE,
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("sup_alpha", valueColumns = ".alpha_bad",
                      backgroundColor = DT::styleEqual("violate", "#fdeaea"),
                      color = DT::styleEqual("violate", "#b00020")) |>
      DT::formatRound(setdiff(names(df)[vapply(df, is.numeric, logical(1))], c("n1", "r1", "n", "r")), 4)
  })

  output$binary_opchar_footnote <- renderUI({
    d <- binary_result()
    req(d)
    df <- binary_opchar()$results
    req(df)
    cols <- names(df)
    midpoint_txt <- if (isTRUE(all.equal(d$inputs$p0L, d$inputs$p0U))) {
      tagList("Single-point null input is used with ", bvar("p", "0"), " = ", sprintf("%.2f", d$inputs$p0), ".")
    } else {
      tagList(
        "When the benchmark range is provided, the classical Simon comparator uses the midpoint ",
        bvar("p", "0"), " = (", bvar("p", "0L"), " + ", bvar("p", "0U"), ") / 2 = ",
        sprintf("%.2f", d$inputs$p0), "."
      )
    }
    lines <- list(
      div(tags$b("n"), " is the total number of subjects required."),
      div(bsub("n", "1"), " is the number of subjects to accrue during stage 1."),
      div(bsub("r", "1"), ", if ", bsub("r", "1"), " or fewer responses are observed during stage 1, the trial is stopped early for futility."),
      div(tags$b("r"), ", if ", tags$b("r"), " or fewer responses are observed by the end of stage two, then no further investigation of the drug is warranted.")
    )
    if (any(c("EN_p0L", "EN_p0", "EN_p0U", "avg_EN") %in% cols)) {
      lines <- c(lines,
                 list(div(bsub("EN", "p0L"), ", ", bsub("EN", "p0"), ", and ", bsub("EN", "p0U"),
                         " are the expected sample sizes at ", bvar("p", "0L"), ", ", bvar("p", "0"), ", and ", bvar("p", "0U"), ", respectively."),
                      div(tags$b("Average EN"), " is the average expected sample size over the benchmark range.")))
    } else if ("EN_p0" %in% cols) {
      lines <- c(lines, list(div(bsub("EN", "p0"), " is the expected sample size at ", bvar("p", "0"), ".")))
    }
    if (any(c("PET_p0L", "PET_p0", "PET_p0U") %in% cols)) {
      lines <- c(lines,
                 list(div(bsub("PET", "p0L"), ", ", bsub("PET", "p0"), ", and ", bsub("PET", "p0U"),
                         " are the probabilities of early termination at ", bvar("p", "0L"), ", ", bvar("p", "0"), ", and ", bvar("p", "0U"), ", respectively.")))
    } else if ("PET_p0" %in% cols) {
      lines <- c(lines, list(div(bsub("PET", "p0"), " is the probability of early termination at ", bvar("p", "0"), ".")))
    }
    if ("sup_alpha" %in% cols) {
      lines <- c(lines, list(div(tags$b("max(alpha)"), " is the maximum type I error over the benchmark range.")))
    }
    if ("power_at_p1" %in% cols) {
      lines <- c(lines, list(div(bsub("power", "p1"), " is the power evaluated at ", bvar("p", "1"), ".")))
    }
    lines <- c(lines, list(
      div("* alpha value highlighted exceeds type I error constraint."),
      div("* ", midpoint_txt)
    ))
    do.call(div, c(list(class = "result-footnote"), lines))
  })

  output$binary_type_I_error_plotly <- renderPlotly({
    d <- binary_result()
    req(d)
    p <- BUDS::plot_type_I_error(d)
    gg <- plotly::ggplotly(p, tooltip = c("x", "y", "colour"))
    x_label <- "True response rate"
    line_traces <- which(vapply(gg$x$data, function(tr) identical(tr[["mode"]], "lines"), logical(1)))
    if (length(line_traces) > 0) {
      gg <- plotly::style(
        gg,
        hovertemplate = paste0(
          x_label, ": %{x:.4f}<br>",
          "Type I error: %{y:.4f}<br>",
          "Design: %{fullData.name}<extra></extra>"
        ),
        traces = line_traces
      )
    }
    gg |>
      plotly::layout(
        legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center"),
        hovermode = "x unified",
        margin = list(b = 80),
        xaxis = list(tickformat = ".4f"),
        yaxis = list(tickformat = ".4%")
      )
  })

  output$binary_summary <- renderUI({
    d <- binary_result()
    req(d, input$binary_detail_design)
    summary_to_html(summary(d, input$binary_detail_design))
  })

  output$binary_rules <- renderTable({
    d <- binary_result()
    req(d, input$binary_detail_design)
    BUDS::table_design_rules(d, input$binary_detail_design)
  }, rownames = FALSE)

  output$binary_opchar_summary <- renderUI({
    req(input$binary_detail_design)
    summary_to_html(summary(binary_opchar(), input$binary_detail_design))
  })

  output$binary_report_download <- downloadHandler(
    filename = function() paste0("binary_report_", Sys.Date(), ".html"),
    content = function(file) {
      d <- binary_result()
      req(d, input$binary_report_design)
      BUDS::export_design_html(
        design = d,
        design_name = input$binary_report_design,
        designs_keep = input$binary_report_keep,
        file = file
      )
    }
  )

  output$tte_value_boxes <- renderUI({
    d <- tte_result()
    if (is.null(d)) {
      return(div(style = "text-align:center; padding:60px; color:#999;",
                 icon("clock", style = "font-size:48px; opacity:0.3;"), br(), br(),
                 "Configure parameters and click ", strong("Run TTE Search")))
    }
    design_df <- matrix_to_df(d$designs)
    first_n2 <- if (nrow(design_df) > 0 && "n2" %in% names(design_df)) as.integer(round(design_df$n2[1])) else "--"
    first_en <- if (nrow(design_df) > 0 && "EN" %in% names(design_df)) round(design_df$EN[1], 1) else "--"
    total_n <- if (!is.null(d$meta$n_total)) d$meta$n_total else "--"
    robust_n <- if (!is.null(d$meta$n_robust)) d$meta$n_robust else "--"
    pooled_n <- if (!is.null(d$meta$n_pooled)) d$meta$n_pooled else "--"
    fluidRow(
      style = "margin-bottom: 16px;",
      column(3, make_vbox("layer-group", total_n, "Generated", "vbox-primary")),
      column(3, make_vbox("shield-halved", robust_n, "Robust", "vbox-info")),
      column(3, make_vbox("diagram-project", pooled_n, "Pooled", "vbox-success")),
      column(3, make_vbox("chart-line", paste0(first_n2, " / ", first_en), "Top n2 / EN", "vbox-secondary"))
    )
  })

  output$tte_timing_badge <- renderUI({
    d <- tte_result()
    req(d)
    tags$span(class = "timing-badge", icon("bolt"), d$meta$execution_time)
  })

  output$tte_design_table <- renderDT({
    d <- tte_result()
    req(d)
    df <- matrix_to_df(d$designs)
    df <- apply_design_labels(df, d)
    op_df <- tte_opchar()$results
    op_df <- apply_design_labels(op_df, d)
    if ("sup_alpha" %in% names(df)) names(df)[names(df) == "sup_alpha"] <- "alpha"
    if ("S0_designed" %in% names(df)) df$S0_designed <- NULL
    for (col in c("n1", "n2")) if (col %in% names(df)) df[[col]] <- as.integer(round(df[[col]]))
    if (all(c("n1", "n2") %in% names(df))) df$ratio <- ifelse(df$n2 == 0, NA_real_, df$n1 / df$n2)
    if ("power_at_S1" %in% names(op_df)) df$power <- unname(op_df$power_at_S1[match(df$Design, op_df$Design)])
    df$.power_check <- if ("power_at_S1" %in% names(op_df)) unname(op_df$power_at_S1[match(df$Design, op_df$Design)]) else NA_real_
    df$.alpha_bad <- ""
    if ("alpha" %in% names(df)) {
      bad <- (round(df$alpha, 2) > round(input$tte_alpha, 2)) |
        (!is.na(df$.power_check) & round(df$.power_check, 2) < round(input$tte_power, 2))
      df$.alpha_bad[bad] <- "violate"
    }
    df <- select_existing(df, c("Design", "c1", "n1", "DA1", "PET", "c2", "n2", "DA2", "EN", "alpha", "power", "ratio",
                                ".power_check", ".alpha_bad"))
    hidden_targets <- which(names(df) %in% c(".power_check", ".alpha_bad")) - 1L
    DT::datatable(
      df,
      colnames = tte_design_colnames(df),
      escape = FALSE,
      options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE,
                     columnDefs = list(
                       list(className = "dt-center", targets = "_all"),
                       list(visible = FALSE, targets = hidden_targets)
                     )),
      rownames = FALSE,
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("alpha", valueColumns = ".alpha_bad",
                      backgroundColor = DT::styleEqual("violate", "#fdeaea"),
                      color = DT::styleEqual("violate", "#b00020")) |>
      DT::formatRound(intersect(c("c1", "c2"), names(df)), 4) |>
      DT::formatRound(intersect(c("DA1", "DA2", "EN", "PET", "alpha", "power", "ratio"), names(df)), 2)
  })

  output$tte_results_footnote <- renderUI({
    d <- tte_result()
    req(d)
    designed_null_txt <- if (identical(d$inputs$input_scale_null, "hazard")) {
      if (isTRUE(all.equal(d$inputs$lambda0L_resolved, d$inputs$lambda0U_resolved))) {
        tagList("The designed benchmark is ", var_label("\u03bb", "0"), " = ", sprintf("%.3f", d$inputs$lambda0_resolved), " for all displayed rows.")
      } else {
        tagList(
          "When the hazard rate benchmark range is supplied, the classical r-KJ comparator uses the midpoint ",
          var_label("\u03bb", "0"), " = (", var_label("\u03bb", "0L"), " + ", var_label("\u03bb", "0U"), ") / 2 = ",
          sprintf("%.3f", d$inputs$lambda0_resolved),
          ", whereas BUDS rows use ", var_label("\u03bb", "0"), " = ", var_label("\u03bb", "0L"), " = ",
          sprintf("%.3f", d$inputs$lambda0L_resolved), "."
        )
      }
    } else {
      if (isTRUE(all.equal(d$inputs$S0L, d$inputs$S0U))) {
        tagList("The designed benchmark is ", var_label("S", "0"), " = ", sprintf("%.2f", d$inputs$S0), " for all displayed rows.")
      } else {
        tagList(
          "When the survival probability benchmark range is supplied, the classical r-KJ comparator uses the midpoint ",
          var_label("S", "0"), " = (", var_label("S", "0L"), " + ", var_label("S", "0U"), ") / 2 = ",
          sprintf("%.2f", d$inputs$S0),
          ", whereas BUDS rows use ", var_label("S", "0"), " = ", var_label("S", "0U"), " = ",
          sprintf("%.2f", d$inputs$S0U), "."
        )
      }
    }
    div(
      class = "result-footnote",
      div(bsub("n", "1"), " and ", bsub("n", "2"), " are the numbers of subjects required for stage 1 and for the full trial, respectively."),
      div(bsub("c", "1"), " is the stage-1 futility boundary. Stop for futility if the stage-1 test statistic exceeds ", bsub("c", "1"), "."),
      div(bsub("c", "2"), " is the final futility boundary. The treatment is declared promising if the final test statistic does not exceed ", bsub("c", "2"), "."),
      div(tags$b("DA1"), " and ", tags$b("DA2"), " are the time of analysis for stage 1 and stage 2, respectively."),
      div(tags$b("alpha"), " is the maximum type I error over the benchmark range."),
      div(
        tags$b("power"), " is the achieved power at ",
        if (identical(d$inputs$input_scale_null, "hazard")) var_label("\u03bb", "1") else var_label("S", "1"), "."
      ),
      div(
        tags$b("EN"), " is the expected sample size for the trial when ",
        if (identical(d$inputs$input_scale_null, "hazard")) "hazard rate is the reported null setting." else "survival probability is the reported null setting."
      ),
      div(tags$b("PET"), " is the probability of early termination."),
      div("* alpha value highlighted exceeds type I error constraint."),
      div("* ", designed_null_txt)
    )
  })

  output$tte_opchar_table <- renderDT({
    d <- tte_result()
    req(d)
    df <- tte_opchar()$results
    df <- apply_design_labels(df, d)
    if ("S0_designed" %in% names(df)) df$S0_designed <- NULL
    if ("alpha_target" %in% names(df)) df$alpha_target <- NULL
    if (!("sup_alpha" %in% names(df)) && "alpha_at_S0" %in% names(df)) df$sup_alpha <- df$alpha_at_S0
    drop_alpha_cols <- intersect(c("alpha_at_S0L", "alpha_at_S0", "alpha_at_S0U"), names(df))
    if (length(drop_alpha_cols) > 0) df <- df[, setdiff(names(df), drop_alpha_cols), drop = FALSE]
    df <- select_existing(df, c("Design", "n1", "c1", "n2", "c2", "DA1", "DA2",
                                "EN_S0L", "EN_S0", "EN_S0U", "avg_EN",
                                "PET_S0L", "PET_S0", "PET_S0U",
                                "sup_alpha", "power_at_S1"))
    df$.alpha_bad <- ""
    if ("sup_alpha" %in% names(df)) df$.alpha_bad[df$sup_alpha > input$tte_alpha] <- "violate"
    hidden_targets <- which(names(df) %in% c(".alpha_bad")) - 1L
    DT::datatable(
      df,
      colnames = tte_opchar_colnames(df),
      escape = FALSE,
      options = list(dom = "t", paging = FALSE, ordering = FALSE, scrollX = TRUE,
                     columnDefs = list(
                       list(className = "dt-center", targets = "_all"),
                       list(visible = FALSE, targets = hidden_targets)
                     )),
      rownames = FALSE,
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("sup_alpha", valueColumns = ".alpha_bad",
                      backgroundColor = DT::styleEqual("violate", "#fdeaea"),
                      color = DT::styleEqual("violate", "#b00020")) |>
      DT::formatRound(intersect(c("c1", "c2"), names(df)), 4) |>
      DT::formatRound(setdiff(names(df)[vapply(df, is.numeric, logical(1))], c("n1", "n2", "c1", "c2")), 2)
  })

  output$tte_opchar_footnote <- renderUI({
    d <- tte_result()
    req(d)
    df <- tte_opchar()$results
    req(df)
    cols <- names(df)
    designed_null_txt <- if (identical(d$inputs$input_scale_null, "hazard")) {
      if (isTRUE(all.equal(d$inputs$lambda0L_resolved, d$inputs$lambda0U_resolved))) {
        tagList("The designed null setting is ", var_label("\u03bb", "0"), " = ", sprintf("%.3f", d$inputs$lambda0_resolved), " for all displayed rows.")
      } else {
        tagList(
          "When the hazard-rate interval is supplied, the classical r-KJ comparator uses the midpoint ",
          var_label("\u03bb", "0"), " = (", var_label("\u03bb", "0L"), " + ", var_label("\u03bb", "0U"), ") / 2 = ",
          sprintf("%.3f", d$inputs$lambda0_resolved),
          ", whereas BUDS rows use ", var_label("\u03bb", "0"), " = ", var_label("\u03bb", "0L"), " = ",
          sprintf("%.3f", d$inputs$lambda0L_resolved), "."
        )
      }
    } else {
      if (isTRUE(all.equal(d$inputs$S0L, d$inputs$S0U))) {
        tagList("The designed null setting is ", var_label("S", "0"), " = ", sprintf("%.2f", d$inputs$S0), " for all displayed rows.")
      } else {
        tagList(
          "When the survival-probability interval is supplied, the classical r-KJ comparator uses the midpoint ",
          var_label("S", "0"), " = (", var_label("S", "0L"), " + ", var_label("S", "0U"), ") / 2 = ",
          sprintf("%.2f", d$inputs$S0),
          ", whereas BUDS rows use ", var_label("S", "0"), " = ", var_label("S", "0U"), " = ",
          sprintf("%.2f", d$inputs$S0U), "."
        )
      }
    }
    lines <- list(
      div(bsub("n", "1"), " and ", bsub("n", "2"), " are the numbers of subjects required for stage 1 and for the full trial, respectively."),
      div(bsub("c", "1"), " is the stage-1 futility boundary. Stop for futility if the stage-1 test statistic exceeds ", bsub("c", "1"), "."),
      div(bsub("c", "2"), " is the final futility boundary. The treatment is declared promising if the final test statistic does not exceed ", bsub("c", "2"), "."),
      div(tags$b("DA1"), " and ", tags$b("DA2"), " are the time of analysis for stage 1 and stage 2, respectively.")
    )
    if (any(c("EN_S0L", "EN_S0", "EN_S0U", "avg_EN") %in% cols)) {
      lines <- c(lines,
                 list(div(bsub("EN", "S0L"), ", ", bsub("EN", "S0"), ", and ", bsub("EN", "S0U"),
                         " are the expected sample sizes at ", bvar("S", "0L"), ", ", bvar("S", "0"), ", and ", bvar("S", "0U"), ", respectively."),
                      div(tags$b("Average EN"), " is the average expected sample size over the benchmark range.")))
    } else if ("EN_S0" %in% cols) {
      lines <- c(lines, list(div(bsub("EN", "S0"), " is the expected sample size at ", bvar("S", "0"), ".")))
    }
    if (any(c("PET_S0L", "PET_S0", "PET_S0U") %in% cols)) {
      lines <- c(lines,
                 list(div(bsub("PET", "S0L"), ", ", bsub("PET", "S0"), ", and ", bsub("PET", "S0U"),
                         " are the probabilities of early termination at ", bvar("S", "0L"), ", ", bvar("S", "0"), ", and ", bvar("S", "0U"), ", respectively.")))
    } else if ("PET_S0" %in% cols) {
      lines <- c(lines, list(div(bsub("PET", "S0"), " is the probability of early termination at ", bvar("S", "0"), ".")))
    }
    if ("sup_alpha" %in% cols) {
      lines <- c(lines, list(div(tags$b("max(alpha)"), " is the maximum type I error over the benchmark range.")))
    }
    if ("power_at_S1" %in% cols) {
      lines <- c(lines, list(div(bsub("power", "S1"), " is the power evaluated at ", bvar("S", "1"), ".")))
    }
    lines <- c(lines, list(
      div("* alpha value highlighted exceeds type I error constraint."),
      div("* ", designed_null_txt)
    ))
    do.call(div, c(list(class = "result-footnote"), lines))
  })

  output$tte_type_I_error_plotly <- renderPlotly({
    d <- tte_result()
    req(d)
    p <- BUDS::plot_type_I_error(d)
    gg <- plotly::ggplotly(p, tooltip = c("x", "y", "colour"))
    x_label <- if (identical(d$inputs$input_scale_null, "hazard")) {
      "True hazard rate"
    } else {
      "True survival probability"
    }
    line_traces <- which(vapply(gg$x$data, function(tr) identical(tr[["mode"]], "lines"), logical(1)))
    if (length(line_traces) > 0) {
      gg <- plotly::style(
        gg,
        hovertemplate = paste0(
          x_label, ": %{x:.4f}<br>",
          "Type I error: %{y:.4f}<br>",
          "Design: %{fullData.name}<extra></extra>"
        ),
        traces = line_traces
      )
    }
    gg |>
      plotly::layout(
        legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center"),
        hovermode = "x unified",
        margin = list(b = 80),
        xaxis = list(tickformat = ".4f"),
        yaxis = list(tickformat = ".4%")
      )
  })

  output$tte_summary <- renderUI({
    d <- tte_result()
    req(d, input$tte_detail_design)
    summary_to_html(summary(d, input$tte_detail_design))
  })

  output$tte_rules <- renderTable({
    d <- tte_result()
    req(d, input$tte_detail_design)
    BUDS::table_design_rules(d, input$tte_detail_design)
  }, rownames = FALSE)

  output$tte_opchar_summary <- renderUI({
    req(input$tte_detail_design)
    summary_to_html(summary(tte_opchar(), input$tte_detail_design))
  })

  output$tte_report_download <- downloadHandler(
    filename = function() paste0("tte_report_", Sys.Date(), ".html"),
    content = function(file) {
      d <- tte_result()
      req(d, input$tte_report_design)
      BUDS::export_design_html(
        design = d,
        design_name = input$tte_report_design,
        designs_keep = input$tte_report_keep,
        file = file
      )
    }
  )
}

shinyApp(ui = ui, server = server)
