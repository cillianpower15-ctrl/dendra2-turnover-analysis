#!/usr/bin/env bash
# Dendra2 turnover analysis launcher for macOS and Linux.
# Usage: ./run_analysis.sh [path/to/Dendra2_Input_Workbook.xlsx]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "R is not installed, or Rscript is not on your PATH."
  echo "Install R from https://cran.r-project.org and try again."
  exit 1
fi

workbook="${1:-}"
if [ -z "$workbook" ]; then
  workbook="$here/Dendra2_Input_Workbook.xlsx"
  echo "No workbook given, using: $workbook"
fi
if [ ! -f "$workbook" ]; then
  echo "Workbook not found: $workbook"
  exit 1
fi

Rscript --vanilla "$here/scripts/run_analysis.R" "$workbook"
