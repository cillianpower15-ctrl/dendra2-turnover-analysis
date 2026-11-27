## One-click orchestrator for the Dendra2 analysis package.

## Where the caller was standing. A workbook given as a relative path is
## relative to this, not to the script folder, so it has to be captured before
## the setwd() below moves us.
.invoked_from <- normalizePath(getwd(), winslash = "/")

.file_arg <- grep("^--file=", commandArgs(), value = TRUE)
.script_dir <- if (length(.file_arg)) {
  dirname(normalizePath(sub("^--file=", "", .file_arg[1]), winslash = "/"))
} else .invoked_from
setwd(.script_dir)

.args <- commandArgs(trailingOnly = TRUE)
if (!length(.args) || !nzchar(.args[1])) {
  if (.Platform$OS.type == "windows") {
    .args <- tryCatch(c(utils::choose.files(
      caption = "Select the Dendra2 Excel workbook",
      filters = matrix(c("Excel workbooks (*.xlsx)", "*.xlsx"), ncol = 2, byrow = TRUE),
      multi = FALSE)), error = function(e) character())
  }
}
if (!length(.args) || !nzchar(.args[1]))
  stop("No workbook was selected. Run RUN_DENDRA2_ANALYSIS.bat and choose an .xlsx file.")

workbook <- .args[1]
if (!grepl("^([A-Za-z]:|[/\\\\])", workbook))
  workbook <- file.path(.invoked_from, workbook)
workbook <- normalizePath(workbook, winslash = "/", mustWork = FALSE)
if (!file.exists(workbook)) stop("Workbook not found: ", workbook)
if (tolower(tools::file_ext(workbook)) != "xlsx")
  stop("The selected input must be an .xlsx workbook: ", workbook)

cat("\n============================================================\n")
cat(" Dendra2 turnover analysis\n")
cat("============================================================\n")
cat("Workbook: ", workbook, "\n\n", sep = "")

required <- c("readxl", "dplyr", "tidyr", "ggplot2", "lme4", "lmerTest",
              "purrr", "patchwork", "emmeans", "nlme")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("First run: installing required R packages. This may take several minutes.\n")
  cat("Packages: ", paste(missing, collapse = ", "), "\n\n", sep = "")
  try(install.packages(missing, repos = "https://cloud.r-project.org",
                       dependencies = c("Depends", "Imports", "LinkingTo")),
      silent = FALSE)
  still_missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing))
    stop("These required packages could not be installed: ",
         paste(still_missing, collapse = ", "),
         "\nConnect to the internet, open R, and run:\ninstall.packages(c(",
         paste(sprintf("\"%s\"", still_missing), collapse = ", "), "))")
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results <- file.path(dirname(workbook), paste0("Dendra2_Results_", stamp))
suffix <- 1
while (dir.exists(results)) {
  results <- file.path(dirname(workbook), paste0("Dendra2_Results_", stamp, "_", suffix))
  suffix <- suffix + 1
}
dir.create(results, recursive = TRUE)
results <- normalizePath(results, winslash = "/")

snapshot <- file.path(results, "INPUT_WORKBOOK_SNAPSHOT.xlsx")
if (!file.copy(workbook, snapshot, overwrite = TRUE))
  stop("Could not create the immutable input-workbook snapshot.")
writeLines(c(
  paste("Selected workbook:", workbook),
  paste("Selected workbook MD5:", unname(tools::md5sum(workbook))),
  paste("Snapshot MD5:", unname(tools::md5sum(snapshot))),
  paste("Workbook modified:", format(file.info(workbook)$mtime, "%Y-%m-%d %H:%M:%S")),
  paste("Analysis started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("R version:", R.version.string),
  paste("Results folder:", results)), file.path(results, "RUN_INFORMATION.txt"))
pkg_versions <- data.frame(
  package = required,
  version = vapply(required, function(p) as.character(packageVersion(p)), character(1)),
  stringsAsFactors = FALSE)
write.csv(pkg_versions, file.path(results, "PACKAGE_VERSIONS.csv"), row.names = FALSE)
writeLines(results, file.path(.script_dir, "LAST_RESULTS_PATH.txt"))

## A failed run leaves a marker; a successful one does not need one.
completed <- FALSE
on.exit({
  if (!completed) {
    writeLines(c("FAILED", paste("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                 "Read the console message above and RUN_INFORMATION.txt."),
               file.path(results, "RUN_FAILED.txt"))
  }
}, add = TRUE)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
## 06 runs before 05 so the report is written last and can reference every
## figure that exists on disk.
steps <- c("01_statistics.R", "02_figures.R", "03_sensitivity_checks.R",
           "04_nonlinear_mixed_model.R", "06_publication_figure.R",
           "07_factorial_figure.R", "08_curate_figures.R", "05_build_report.R")
labels <- c("Statistics and input QC", "Publication figures", "Sensitivity checks",
            "Nonlinear turnover-rate model", "Composite main figure",
            "Per-background figure (crossed designs)",
            "Curating figures for the paper", "Plain-English report")

for (i in seq_along(steps)) {
  cat(sprintf("\n[%d of %d] %s...\n", i, length(steps), labels[i]))
  status <- system2(rscript,
    args = c("--vanilla", shQuote(file.path(.script_dir, steps[i])),
             shQuote(workbook), shQuote(results)))
  if (!identical(status, 0L))
    stop("Step failed: ", labels[i], ". The completed outputs and logs remain in:\n", results)
}

## Only the outputs that survive curation (step 08 deletes the intermediates).
required_outputs <- c(
  "Dendra2_Analysis_Report.html", "ANALYSIS_SUMMARY.txt",
  "INPUT_WORKBOOK_SNAPSHOT.xlsx", "RUN_INFORMATION.txt",
  "PACKAGE_VERSIONS.csv",
  file.path("tables", "00_publication_readiness.csv"),
  file.path("tables", "14_nonlinear_model_comparisons.csv"),
  file.path("tables", "12_nonlinear_group_rates.csv"),
  file.path("figures", "Figure1_MAIN.pdf"))
required_paths <- file.path(results, required_outputs)
bad <- !file.exists(required_paths) | is.na(file.info(required_paths)$size) |
  file.info(required_paths)$size <= 0
if (any(bad))
  stop("Output verification failed. Missing or empty files:\n",
       paste(required_outputs[bad], collapse = "\n"))

## ---- tidy the results folder ---------------------------------------------
## Verification is done, so nothing downstream depends on these paths.
## Provenance and console logs matter for the methods section but are not what
## anyone opens first: they move one level down, leaving the top of the folder
## holding the report, the summary, the figures and the tables and nothing else.
detail <- file.path(results, "run details")
dir.create(detail, showWarnings = FALSE)
for (f in c("RUN_INFORMATION.txt", "PACKAGE_VERSIONS.csv",
            "INPUT_WORKBOOK_SNAPSHOT.xlsx", "METHODS_REPORTING_TEMPLATE.txt",
            "analysis_log.txt", "sensitivity_log.txt", "nonlinear_model_log.txt")) {
  src <- file.path(results, f)
  if (file.exists(src)) file.rename(src, file.path(detail, f))
}
## figure legends belong beside the figures they describe
for (f in c("Figure1_MAIN_LEGEND.txt", "Figure1_MAIN_LEGEND_SHORT.txt")) {
  src <- file.path(results, f)
  if (file.exists(src)) file.rename(src, file.path(results, "figures", f))
}

## Completion is reported on screen and by RUN_INFORMATION.txt; a status file
## and a per-file checksum manifest only add clutter to the results folder.
completed <- TRUE
cat("\n============================================================\n")
cat(" ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Results: ", results, "\n", sep = "")
cat("Start with: Dendra2_Analysis_Report.html\n")
