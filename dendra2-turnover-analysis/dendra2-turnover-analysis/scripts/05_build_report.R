## Build a concise, readable HTML and text report from the analysis outputs.

source("00_setup.R")

read_if <- function(path) {
  if (file.exists(path)) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) else NULL
}
table_path <- function(name) file.path(TABLES, name)
qc <- read_if(table_path("00_input_QC.csv"))
photo_qc <- read_if(table_path("00_photoconversion_QC.csv"))
readiness <- read_if(table_path("00_publication_readiness.csv"))
exclusions <- read_if(table_path("00_excluded_records.csv"))
manifest <- read_if(table_path("00_statistics_performed.csv"))
headline <- read_if(table_path("09_headline_summary.csv"))
nonlinear_cmp <- read_if(table_path("14_nonlinear_model_comparisons.csv"))
nonlinear_rates <- read_if(table_path("12_nonlinear_group_rates.csv"))

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", formatC(p, digits = 3, format = "f")))
}
## The threshold is stated once in the plan; repeating it after every p-value
## made the report read like a compliance printout.
sig_text <- function(p) {
  ifelse(is.na(p), "not available",
         ifelse(p < 0.05, "significant", "not significant"))
}
esc <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
html_table <- function(df, digits = 4) {
  if (is.null(df) || !nrow(df)) return("<p><em>Not available.</em></p>")
  d <- df
  d[] <- lapply(d, function(x) if (is.numeric(x)) format(signif(x, digits), trim = TRUE) else as.character(x))
  header <- paste0("<th>", esc(names(d)), "</th>", collapse = "")
  rows <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

n_tbl <- dendra %>% distinct(Group, Slice) %>% count(Group, name = "n_slices")
ok_status <- c("Recorded", "Yes", "Passed")
warnings <- if (is.null(readiness)) "The pre-submission checklist table was not produced." else
  paste(readiness$item[!readiness$status %in% ok_status],
        readiness$status[!readiness$status %in% ok_status], sep = " - ")
warning_html <- if (!length(warnings)) {
  "<div class='good'><strong>Input QC:</strong> all slices passed the automated checks.</div>"
} else {
  paste0("<div class='warn'><strong>Still to record before submission:</strong><ul>",
         paste0("<li>", esc(warnings), "</li>", collapse = ""), "</ul></div>")
}

headline_text <- character()
if (!is.null(headline) && nrow(headline)) {
  headline_text <- paste0("- ", headline$test, ": p = ", fmt_p(headline$p),
                          " (", sig_text(headline$p), ")")
}
if (!is.null(nonlinear_cmp) && nrow(nonlinear_cmp)) {
  headline_text <- c(headline_text,
    paste0("- ", nonlinear_cmp$comparison, ": p = ", fmt_p(nonlinear_cmp$p),
           " (", sig_text(nonlinear_cmp$p), ")"))
}

stats_text <- c(
  "DENDRA2 ANALYSIS: STATISTICS PERFORMED",
  "=======================================",
  paste("Workbook:", WORKBOOK),
  paste("Reference group:", GROUPS[1]),
  paste("Primary analysis:", PRIMARY_ANALYSIS),
  paste("Preparation handling:", PREP_STATUS),
  paste("Biological replicates:", paste(paste0(n_tbl$Group, " n=", n_tbl$n_slices), collapse = "; ")),
  "",
  if (!is.null(manifest)) paste0(manifest$order, ". ", manifest$analysis, " -- ", manifest$method) else
    "See analysis_log.txt for the full model specification.")
## stats_text is summarised in the HTML report; no separate file needed.

summary_text <- c(
  "DENDRA2 ANALYSIS SUMMARY",
  "=========================",
  paste("Workbook:", WORKBOOK),
  paste("Analysed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Groups:", paste(paste0(n_tbl$Group, " n=", n_tbl$n_slices), collapse = "; ")),
  paste("Reference group:", GROUPS[1]),
  paste("Primary analysis:", PRIMARY_ANALYSIS),
  paste("Preparation handling:", PREP_STATUS),
  paste0("Baselines: GFP day ", BASELINE_GFP,
         "; RFP day ", BASELINE_RFP_DISPLAY, " post-conversion"),
  "",
  "HEADLINE TESTS",
  if (length(headline_text)) headline_text else "No headline table was produced.",
  "",
  ## Two facts about how the model was fitted, which a reader needs. The
  ## general warnings about not over-reading a p-value that used to sit here
  ## applied to any analysis ever run and told nobody anything about this one.
  "MODEL NOTES",
  paste0("The day-", BASELINE_RFP_DISPLAY,
         " post-conversion RFP anchor is excluded from the nonlinear likelihood."),
  "K is constrained above zero through log(K); plateau P is constrained to 0-1 through logit(P).")
writeLines(summary_text, file.path(RESULTS, "ANALYSIS_SUMMARY.txt"))

methods_text <- c(
  "METHODS REPORTING TEMPLATE — EDIT BRACKETED DETAILS BEFORE USE",
  "=============================================================",
  "",
  paste0("Longitudinal Dendra2 measurements were analysed in R (", R.version.string,
         "). The prespecified primary endpoint was the between-group difference in ",
         "the post-photoconversion RFP turnover-rate constant. Fraction retained was modelled as ",
         "(1-P)exp(-Kt)+P using nonlinear mixed effects, with log(K) and logit(P) parameterisations ",
         "to constrain K>0 and 0<P<1, and a slice-level random effect on log(K). ",
         "The normalised post-conversion baseline defined fraction(0)=1 structurally but was excluded ",
         "from the likelihood. Nested fixed-effect models were fitted by maximum likelihood and compared ",
         "with a likelihood-ratio test; pairwise rate ratios were Holm adjusted. ",
         if (ADJUST_FOR_PREP) "Preparation ID was included as a fixed blocking factor. "
         else "[Preparation/animal blocking was unavailable or non-identifiable; explain this limitation.] ",
         "Secondary linear mixed models used repeated measurements within slice and retained a random ",
         "Day slope whenever estimable, with random-intercept fallback only for singular/non-convergent fits. ",
         "All exclusions were recorded before analysis. [Add microscope, registration, retained-Z-plane, ",
         "signal/background and photoconversion QC details from the QC tables.]"),
  "",
  "Report the primary likelihood-ratio statistic, degrees of freedom and exact p-value,",
  "adjusted group K estimates with 95% confidence intervals, Holm-adjusted rate ratios,",
  "slice counts, preparation/animal counts, exclusions and leave-one-slice-out results.")
writeLines(methods_text, file.path(RESULTS, "METHODS_REPORTING_TEMPLATE.txt"))

headline_for_html <- headline
if (!is.null(headline_for_html) && "p" %in% names(headline_for_html))
  headline_for_html$p <- fmt_p(headline_for_html$p)
if (!is.null(nonlinear_cmp) && "p" %in% names(nonlinear_cmp))
  nonlinear_cmp$p <- fmt_p(nonlinear_cmp$p)

html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>Dendra2 Analysis Report</title>",
"<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:32px auto;color:#17202a;line-height:1.45;padding:0 18px}",
"h1{color:#1f4e78;margin-bottom:4px}h2{color:#c55a11;border-bottom:2px solid #f4d7c6;padding-bottom:5px;margin-top:30px}",
".meta{background:#f2f6f9;padding:14px 18px;border-left:5px solid #1f4e78}.good{background:#eaf5e4;padding:12px;border-left:5px solid #4f8a3b}",
".warn{background:#fff4dd;padding:12px;border-left:5px solid #d88b00}table{border-collapse:collapse;width:100%;font-size:14px;margin:10px 0 22px}",
"th{background:#1f4e78;color:white;text-align:left;padding:7px}td{border-bottom:1px solid #d9e1e7;padding:7px;vertical-align:top}",
"tr:nth-child(even){background:#f7f9fa}.fig{max-width:100%;border:1px solid #d9e1e7;margin:8px 0 18px}",
"code{background:#eef2f4;padding:2px 4px}.small{font-size:13px;color:#59636d}</style></head><body>",
## The old subtitle announced the report as machine-generated, which is the
## first thing any reader saw. The run date and workbook below already say
## everything that line was there to say.
"<h1>Dendra2 turnover analysis</h1>",
"<div class='meta'><strong>Workbook:</strong> ", esc(basename(WORKBOOK)), "<br><strong>Analysed:</strong> ",
esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "<br><strong>Groups:</strong> ",
esc(paste(paste0(n_tbl$Group, " n=", n_tbl$n_slices), collapse = "; ")), "<br><strong>Reference:</strong> ",
esc(GROUPS[1]), "<br><strong>Primary analysis:</strong> ", esc(PRIMARY_ANALYSIS),
"<br><strong>Preparation handling:</strong> ", esc(PREP_STATUS),
"<br><strong>Baselines:</strong> GFP day ", BASELINE_GFP, "; RFP day ",
BASELINE_RFP_DISPLAY, " post-conversion</div>",
warning_html,
"<h2>Before reporting these results</h2>", html_table(readiness),
"<h2>Headline statistical results</h2>", html_table(headline_for_html),
"<h2>Primary nonlinear RFP turnover-rate model</h2><p>M1 vs M2 tests whether the exponential rate constant K differs across groups. Models use maximum likelihood, log(K) constrains K above zero, logit(P) constrains the plateau to 0–1, and the structural t=0 anchor is not fitted as an observed datum.</p>",
html_table(nonlinear_cmp), html_table(nonlinear_rates),
"<h2>What statistics were performed?</h2>", html_table(manifest),
"<h2>Input quality control</h2>", html_table(qc),
"<h3>Photoconversion QC</h3>", html_table(photo_qc),
"<h3>Excluded records</h3>", html_table(exclusions),
"<h2>Figures</h2>",
"<p class='small'>Curated down to what a paper actually carries: one main figure plus a ",
"small number of supplementary figures. Every other analysis view is still available in ",
"<code>figures/other_analyses/</code> &mdash; nothing was deleted, just kept out of the way.</p>",
if (file.exists(file.path(FIGURES, "Figure1_MAIN.png"))) {
  leg <- if (file.exists(file.path(RESULTS, "Figure1_MAIN_LEGEND.txt")))
    "<code>Figure1_MAIN_LEGEND.txt</code>" else "the analysis log"
  paste0("<h3>Figure 1 (main)</h3>",
         "<p class='small'>This is the figure intended for the manuscript. A ready-to-paste ",
         "legend with the current numbers is in ", leg, "; the vector version for the ",
         "journal is <code>figures/Figure1_MAIN.pdf</code>. Individual panels are saved ",
         "separately (<code>Figure1_MAIN_panelA.pdf</code> etc.) for building a custom layout.</p>",
         "<img class='fig' src='figures/Figure1_MAIN.png'>")
} else "",
if (file.exists(file.path(FIGURES, "FigureS1_per_replicate_metrics.png")))
  "<h3>Figure S1 &mdash; per-replicate half-life and AUC</h3><p class='small'>Model-free, per-slice values behind the main figure's model-based estimates.</p><img class='fig' src='figures/FigureS1_per_replicate_metrics.png'>" else "",
if (file.exists(file.path(FIGURES, "FigureS2_GFP_synthesis.png")))
  "<h3>Figure S2 &mdash; unconverted (green) channel</h3><p class='small'>Residual unconverted signal plus new synthesis; see the axis label for what this is and is not.</p><img class='fig' src='figures/FigureS2_GFP_synthesis.png'>" else "",
if (file.exists(file.path(FIGURES, "FigureS3_model_diagnostics.png")))
  "<h3>Figure S3 &mdash; model diagnostics</h3><img class='fig' src='figures/FigureS3_model_diagnostics.png'>" else "",
if (file.exists(file.path(FIGURES, "FigureS4_decay_model_sensitivity.png")))
  "<h3>Figure S4 &mdash; sensitivity to the decay model</h3><img class='fig' src='figures/FigureS4_decay_model_sensitivity.png'>" else "",
"<h2>How to report this responsibly</h2><ul><li>Treat each slice as the repeated-measures unit; document whether slices are nested within the same animal/preparation.</li>",
"<li>Report estimates, confidence intervals and exact adjusted p-values alongside individual-slice data.</li>",
"<li>Check model diagnostics and leave-one-slice-out results before making a strong claim.</li>",
"<li>The nonlinear M1 vs M2 likelihood-ratio test is the prespecified primary test of a group difference in RFP turnover rate.</li>",
"<li>Resolve or explicitly justify every non-PASS item in the publication-readiness checklist.</li></ul>",
"<p class='small'>Full numerical output remains in <code>tables/</code>; detailed transcripts are in the log files.</p></body></html>")
writeLines(html, file.path(RESULTS, "Dendra2_Analysis_Report.html"), useBytes = TRUE)
message("Report written to ", file.path(RESULTS, "Dendra2_Analysis_Report.html"))
