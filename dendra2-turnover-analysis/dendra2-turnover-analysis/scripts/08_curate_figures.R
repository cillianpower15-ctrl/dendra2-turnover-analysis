## 08_curate_figures.R -- cut the figure set down to what a paper actually uses.

source("00_setup.R")

other_dir <- file.path(FIGURES, "other_analyses")
dir.create(other_dir, showWarnings = FALSE)

has <- function(name) file.exists(file.path(FIGURES, paste0(name, ".png")))
mv <- function(from, to) {
  for (ext in c("pdf", "png")) {
    src <- file.path(FIGURES, paste0(from, ".", ext))
    if (file.exists(src)) file.rename(src, file.path(FIGURES, paste0(to, ".", ext)))
  }
}
archive <- function(name) {
  for (ext in c("pdf", "png")) {
    src <- file.path(FIGURES, paste0(name, ".", ext))
    if (file.exists(src)) file.rename(src, file.path(other_dir, paste0(name, ".", ext)))
  }
}

## ---- MAIN FIGURE -------------------------------------------------------
if (has("Figure5_by_background")) {
  ## filenames follow the C/D panel tags set in 07
  mv("Figure5_by_background",              "Figure1_MAIN")
  mv("Figure5A_decay_by_background",       "Figure1_MAIN_panelC")
  mv("Figure5B_HSPB1_effect_by_background", "Figure1_MAIN_panelD")
  main_used <- "background"
  ## 06's legend describes a composite deleted below, so it goes too
  unlink(file.path(RESULTS, "Figure4_LEGEND.txt"))
} else if (has("Figure4_MAIN_publication")) {
  mv("Figure4_MAIN_publication", "Figure1_MAIN")
  mv("Figure4A_turnover_curves", "Figure1_MAIN_panelA")
  mv("Figure4B_halflife",        "Figure1_MAIN_panelB")
  mv("Figure4C_new_synthesis",   "Figure1_MAIN_panelC")
  mv("Figure4D_effect_size",     "Figure1_MAIN_panelD")
  main_used <- "reference"
  f4l <- file.path(RESULTS, "Figure4_LEGEND.txt")
  if (file.exists(f4l))
    invisible(file.rename(f4l, file.path(RESULTS, "Figure1_MAIN_LEGEND.txt")))
} else {
  main_used <- "none"
}

## ---- SUPPLEMENTARY -------------------------------------------------------
## per-slice half-life/AUC: the model-free spread behind the main figure's fit
mv("Figure2_per_replicate_metrics", "FigureS1_per_replicate_metrics")

## GFP channel, only if the main figure is the by-background layout (no GFP there)
if (main_used == "background" && has("Figure4C_new_synthesis"))
  mv("Figure4C_new_synthesis", "FigureS2_GFP_synthesis")

## ---- everything else is deleted, not archived ----------------------------
## All reproducible by re-running, so copies would only clutter the folder.
drop <- function(name) {
  for (ext in c("pdf", "png")) {
    f <- file.path(FIGURES, paste0(name, ".", ext))
    if (file.exists(f)) unlink(f)
  }
}
for (nm in c("Figure1_decay_trajectories", "Figure3_nonlinear_mixed_model",
             "Figure4_MAIN_publication", "Figure4A_turnover_curves",
             "Figure4B_halflife", "Figure4C_new_synthesis", "Figure4D_effect_size",
             "Figure5_by_background", "Figure5A_decay_by_background",
             "Figure5B_HSPB1_effect_by_background",
             "FigureS1_model_diagnostics", "FigureS2_halflife_model_sensitivity"))
  drop(nm)
unlink(other_dir, recursive = TRUE)

## keep only the tables the figures and methods text actually cite
keep_tables <- c("00_input_QC.csv", "00_publication_readiness.csv",
                 "00_excluded_records.csv", "04_pairwise_contrasts.csv",
                 "06_per_replicate_halflife.csv", "07_per_replicate_AUC.csv",
                 "08_group_comparisons.csv", "12_nonlinear_group_rates.csv",
                 "14_nonlinear_model_comparisons.csv",
                 "17_within_background_HSPB1_vs_Control.csv",
                 "18_perday_RFP_HSPB1_vs_Control.csv",
                 "19_background_x_construct_interaction.csv")
for (f in setdiff(list.files(TABLES), keep_tables))
  unlink(file.path(TABLES, f))
for (f in c("OUTPUT_MANIFEST.csv", "RUN_STATUS.txt", "SESSION_INFO.txt",
            "STATISTICS_PERFORMED.txt", "fitted_objects.rds"))
  unlink(file.path(RESULTS, f))

kept <- list.files(FIGURES, pattern = "\\.png$")
message("\nCurated figure set (", length(kept), " figures) in ", normalizePath(FIGURES), ":")
for (k in sort(kept)) message("  ", sub("\\.png$", "", k))
message(length(list.files(TABLES)), " tables kept.")
