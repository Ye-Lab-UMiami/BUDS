# BUDS

Benchmark Uncertainty Design Selection (BUDS) for two-stage single-arm phase II clinical trial designs with binary and time-to-event endpoints.

This package accompanies the manuscript **BUDS: Benchmark Uncertainty Design Selection for Two-Stage Single-Arm Phase II Trials**. It implements classical two-stage phase II designs and BUDS design-selection criteria for settings where the historical benchmark used in trial planning is uncertain.

## Features

- Binary endpoint designs based on Simon two-stage designs.
- Time-to-event endpoint designs based on restricted Kwak-Jung designs with
  restricted follow-up.
- Interval-null robust design selection using least-regret and average-EN
  criteria.
- Operating-characteristic evaluation across single and interval benchmark
  settings.
- Reporting helpers for tables, type I error plots, and self-contained HTML
  reports.
- Optional Shiny app interface for interactive design search.

## R Package

The R package is stored in the `BUDS/` subfolder of this repository. Installing the package gives access to the design functions, but it does not launch the Shiny app.

Install the R package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("Ye-Lab-UMiami/BUDS", subdir = "BUDS")
```

Alternatively, download the GitHub repository as a ZIP file, unzip it, and install from the local package folder. GitHub ZIP downloads may create a folder named `INRDesign-main`:

```r
install.packages(
  "/path/to/BUDS-main/BUDS",
  repos = NULL,
  type = "source"
)
```

Then load the package:

```r
library(BUDS)
```

## Package Documentation

After installation, help pages are available with:

```r
?binary_BUDS
?tte_BUDS
?opchar_binary_BUDS
?opchar_tte_BUDS
?table_design_opchar
?plot_type_I_error
?export_design_html
help(package = "BUDS")
```

## Package Examples

### Binary Endpoint Example

```r
bin_design <- binary_BUDS(
  alpha = 0.05,
  power = 0.80,
  p0L = 0.08,
  p0U = 0.12,
  p1 = 0.25,
  n.ub = 150,
  robust_objective = c("least_regret", "avg_en")
)

bin_design
opchar_binary_BUDS(bin_design)
plot_type_I_error(bin_design)
```

For a single benchmark setting, set `p0L = p0U` or supply `p0` directly.

### Time-to-Event Endpoint Example

Use survival probabilities:

```r
tte_design <- tte_BUDS(
  alpha = 0.05,
  beta = 0.10,
  S0L = 0.48,
  S0U = 0.52,
  S1 = 0.70,
  x0 = 1,
  rate = 15,
  n.ub = 150,
  inr_objective = c("avg_en", "least_regret")
)

tte_design
opchar_tte_BUDS(tte_design)
plot_type_I_error(tte_design)
```

Alternatively, specify the alternative using a hazard ratio:

```r
tte_design_hr <- tte_BUDS(
  alpha = 0.05,
  beta = 0.10,
  S0L = 0.48,
  S0U = 0.52,
  HR = 2,
  x0 = 1,
  rate = 15,
  n.ub = 150
)
```

Hazard-rate inputs are also supported through `lambda0`, `lambda0L`,
`lambda0U`, and `lambda1`.

### HTML Report

```r
export_design_html(
  design = tte_design,
  design_name = "BUDS (Least Regret)",
  file = "tte_report.html"
)
```

## Shiny App

The Shiny app is separate from the R package. It is stored at the repository root as `app.R`, so you need the full repository, not only the installed package, to run the app locally.

After downloading the repository as a ZIP file from GitHub, the folder may be named `BUDS-main`. In R, set the working directory to that repository folder, then run the app:

```r
setwd("/path/to/BUDS-main")
shiny::runApp("app.R")
```

If you clone the repository with Git, run:

```bash
git clone https://github.com/Ye-Lab-UMiami/BUDS.git
cd BUDS
R -q -e 'shiny::runApp("app.R")'
```

## Package Authors

The R package authors are:

- Rebecca Irlmeier
- Zhuoli Jin

The package maintainer is Zhuoli Jin.

## Manuscript Authors and Affiliations

Manuscript title: **BUDS: Benchmark Uncertainty Design Selection for Two-Stage Single-Arm Phase II Trials**

Running head: **BUDS for Two-Stage Single-Arm Phase II Trials**

Authors: Rebecca Irlmeier, Zhuoli Jin, and Fei Ye.

Affiliations:

- Biostatistics and Bioinformatics Shared Resource, Sylvester Comprehensive
  Cancer Center, Miami, Florida.
- Sylvester Comprehensive Cancer Center, Miami, Florida.
- Division of Biostatistics and Bioinformatics, Department of Public Health
  Sciences, University of Miami, Miami, Florida.

Corresponding author:

Fei Ye, PhD  
Department of Public Health Sciences  
Miller School of Medicine, University of Miami  
Sylvester Comprehensive Cancer Center  
Don Soffer Clinical Research Center

## References

Simon, R. (1989). Optimal two-stage designs for phase II clinical trials.
*Controlled Clinical Trials*, 10(1), 1-10.
https://doi.org/10.1016/0197-2456(89)90015-9

Belin, L., De Rycke, Y., & Broet, P. (2017). A two-stage design for phase II
trials with time-to-event endpoint using restricted follow-up. *Contemporary
Clinical Trials Communications*, 8, 127-134.
https://doi.org/10.1016/j.conctc.2017.09.010

Irlmeier, R., Jin, Z., & Ye, F. BUDS: Benchmark Uncertainty Design Selection for Two-Stage Single-Arm Phase II Trials.
