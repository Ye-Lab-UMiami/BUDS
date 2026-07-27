#!/bin/bash
# ============================================================
# install.sh -- Install the BUDS R package
#
# Run from the directory containing the BUDS/ folder:
#   bash install.sh
# ============================================================

set -e
cd "$(dirname "$0")"

echo "============================================"
echo " Installing BUDS R package"
echo "============================================"
echo ""

# Check R is available
if ! command -v R &> /dev/null; then
    echo "ERROR: R not found. Please install R first."
    exit 1
fi

# Check Rcpp is installed
Rscript -e 'if (!requireNamespace("Rcpp", quietly=TRUE)) install.packages("Rcpp", repos="https://cloud.r-project.org")'

# Install the package
R CMD INSTALL BUDS

echo ""
echo "[OK] BUDS installed successfully."
echo ""
echo "Test it in R:"
echo '  library(BUDS)'
echo '  design <- binary_iBUDS(alpha=0.05, power=0.80,'
echo '              p0=0.20, p1=0.35,'
echo '              p0L=0.18, p0U=0.22, n.ub=150)'
echo '  design'
