## 00_setup.R -- shared configuration, validated data loading and figure theme.
## User-editable options live in the workbook; this file is not edited by users.

## ---- 1. INPUT AND OUTPUT LOCATIONS -----------------------------------------
SHEET <- "Master Data"
.args <- commandArgs(trailingOnly = TRUE)
MASTER_WORKBOOK <- "Dendra2_Input_Workbook.xlsx"
WORKBOOK <- if (length(.args) >= 1 && nzchar(.args[1])) .args[1] else MASTER_WORKBOOK
WORKBOOK <- normalizePath(WORKBOOK, winslash = "/", mustWork = FALSE)
if (!file.exists(WORKBOOK))
  stop("Workbook not found:\n  ", WORKBOOK,
       "\nUse RUN_DENDRA2_ANALYSIS.bat and select the workbook when asked.")

RESULTS <- if (length(.args) >= 2 && nzchar(.args[2])) .args[2] else
  file.path(dirname(WORKBOOK), paste0("Dendra2_Results_", format(Sys.time(), "%Y%m%d_%H%M%S")))
RESULTS <- normalizePath(RESULTS, winslash = "/", mustWork = FALSE)
TABLES <- file.path(RESULTS, "tables")
FIGURES <- file.path(RESULTS, "figures")
for (d in c(RESULTS, TABLES, FIGURES)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- 2. PACKAGES ------------------------------------------------------------
.pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "lme4", "lmerTest", "purrr",
           "patchwork", "emmeans", "nlme")
.missing <- .pkgs[!vapply(.pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing)) stop("Missing R packages: ", paste(.missing, collapse = ", "))
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(lme4); library(lmerTest); library(purrr); library(patchwork); library(emmeans)
})

## ---- 3. WORKBOOK SETTINGS ---------------------------------------------------
REFERENCE_GROUP <- "Control"
BASELINE_RFP <- 1.1
BASELINE_GFP <- 1
MIN_DECAY_R2 <- 0.8
MIN_GFP_LOSS <- 0.20
MIN_RFP_FOLD <- 2.0
MIN_SBR <- 2.0
PRIMARY_ANALYSIS <- "Constrained nonlinear RFP rate model"
RANDOM_SLOPE_STRATEGY <- "Use if estimable"

## ---- Backgrounds excluded from the analysis -------------------------------
## Default: exclude nothing. Set "Exclude backgrounds" in the Analysis Settings
## sheet to drop groups whose name contains a given text, comma separated.
## This is deliberately NOT hardcoded: a hardcoded value silently removes groups
## whose names happen to contain that text, for example a promoter name.
EXCLUDE_BACKGROUNDS <- character()

.as_yes <- function(x, default = TRUE) {
  y <- tolower(trimws(as.character(x)))
  blank <- is.na(x) | !nzchar(y)
  out <- y %in% c("yes", "y", "true", "1", "include", "included")
  out[blank] <- default
  out
}

.sheets <- excel_sheets(WORKBOOK)
if (!SHEET %in% .sheets)
  stop("The workbook has no sheet named '", SHEET, "'. Select the supplied input workbook.")

if ("Analysis Settings" %in% .sheets) {
  .set <- suppressMessages(read_excel(WORKBOOK, sheet = "Analysis Settings", skip = 2))
  if (ncol(.set) >= 2) {
    .set <- .set[, 1:2]; names(.set) <- c("Setting", "Value")
    .setting <- function(name, default = NA_character_) {
      hit <- .set$Value[tolower(trimws(.set$Setting)) == tolower(name)]
      if (!length(hit) || is.na(hit[1]) || !nzchar(trimws(as.character(hit[1])))) default else hit[1]
    }
    REFERENCE_GROUP <- as.character(.setting("Reference group", REFERENCE_GROUP))
    BASELINE_GFP <- suppressWarnings(as.numeric(.setting("GFP baseline day", BASELINE_GFP)))
    BASELINE_RFP <- suppressWarnings(as.numeric(.setting("RFP baseline day", BASELINE_RFP)))
    MIN_DECAY_R2 <- suppressWarnings(as.numeric(.setting("Minimum acceptable decay R2", MIN_DECAY_R2)))
    MIN_GFP_LOSS <- suppressWarnings(as.numeric(.setting("Minimum GFP photoconversion loss fraction", MIN_GFP_LOSS)))
    MIN_RFP_FOLD <- suppressWarnings(as.numeric(.setting("Minimum RFP photoconversion fold increase", MIN_RFP_FOLD)))
    MIN_SBR <- suppressWarnings(as.numeric(.setting("Minimum signal/background ratio", MIN_SBR)))
    PRIMARY_ANALYSIS <- as.character(.setting("Primary analysis", PRIMARY_ANALYSIS))
    RANDOM_SLOPE_STRATEGY <- as.character(.setting("Random slope strategy", RANDOM_SLOPE_STRATEGY))
    .excl <- .setting("Exclude backgrounds", "")
    .excl <- trimws(unlist(strsplit(as.character(.excl), ",", fixed = TRUE)))
    EXCLUDE_BACKGROUNDS <- .excl[nzchar(.excl)]
  }
}
if (!is.finite(BASELINE_GFP)) BASELINE_GFP <- 1
if (!is.finite(BASELINE_RFP)) BASELINE_RFP <- 1.1
if (!is.finite(MIN_DECAY_R2)) MIN_DECAY_R2 <- 0.8
if (!is.finite(MIN_GFP_LOSS)) MIN_GFP_LOSS <- 0.20
if (!is.finite(MIN_RFP_FOLD)) MIN_RFP_FOLD <- 2.0
if (!is.finite(MIN_SBR)) MIN_SBR <- 2.0

## ---- 4. LOAD, NORMALISE AND VALIDATE ---------------------------------------
EXAMPLE_PREFIXES <- c("e.g.", "eg ", "example")
.is_example <- function(x) {
  x <- tolower(trimws(x))
  Reduce(`|`, lapply(EXAMPLE_PREFIXES, function(p) startsWith(x, tolower(p))))
}
.norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))
.optional_name <- function(df, target, aliases, default) {
  if (target %in% names(df)) return(df)
  hit <- which(.norm_name(names(df)) %in% .norm_name(aliases))
  if (length(hit)) names(df)[hit[1]] <- target else df[[target]] <- default
  df
}

.raw <- suppressMessages(read_excel(WORKBOOK, sheet = SHEET, skip = 1))
if (ncol(.raw) < 7)
  stop("'Master Data' needs at least seven columns: Group, Slice ID, RecordingDay, ",
       "GFP_RawSignal, GFP_dFF0, RFP_RawSignal and RFP_dFF0.")

.raw <- .raw %>% rename(Group = 1, Slice = 2, Day = 3,
                        GFP_Raw = 4, GFP_dFF0 = 5, RFP_Raw = 6, RFP_dFF0 = 7)
.raw <- .optional_name(.raw, "PreparationID",
                       c("PreparationID", "Preparation ID", "AnimalID", "Animal ID", "PrepID"), NA_character_)
.raw <- .optional_name(.raw, "Include", c("Include", "Included", "Use"), "Yes")
.raw <- .optional_name(.raw, "ExclusionReason",
                       c("ExclusionReason", "Exclusion Reason", "Reason excluded"), "")
.raw <- .optional_name(.raw, "RetainedZPlanes",
                       c("RetainedZPlanes", "Retained Z Planes", "Z planes retained"), NA_real_)
.raw <- .optional_name(.raw, "RegistrationQC",
                       c("RegistrationQC", "Registration QC", "Registration status"), "Not recorded")
.raw <- .optional_name(.raw, "SignalBackgroundRatio",
                       c("SignalBackgroundRatio", "Signal Background Ratio", "SBR", "SNR"), NA_real_)

dendra_all <- .raw %>%
  filter(!is.na(Group), nzchar(trimws(Group)), !.is_example(Group), !is.na(Day)) %>%
  mutate(
    Group = trimws(as.character(Group)), Slice = trimws(as.character(Slice)),
    PreparationID = trimws(ifelse(is.na(PreparationID), "", as.character(PreparationID))),
    Include = .as_yes(Include),
    ExclusionReason = trimws(ifelse(is.na(ExclusionReason), "", as.character(ExclusionReason))),
    RegistrationQC = trimws(ifelse(is.na(RegistrationQC) | !nzchar(trimws(as.character(RegistrationQC))),
                                    "Not recorded", as.character(RegistrationQC))),
    across(c(Day, GFP_Raw, GFP_dFF0, RFP_Raw, RFP_dFF0, RetainedZPlanes,
             SignalBackgroundRatio), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(!is.na(Day))

if (!nrow(dendra_all)) stop("No usable rows were found in 'Master Data'.")
if (any(is.na(dendra_all$Slice) | !nzchar(dendra_all$Slice)))
  stop("Every data row needs a Slice ID.")

include_check <- dendra_all %>% group_by(Group, Slice) %>%
  summarise(n_include_values = n_distinct(Include), .groups = "drop") %>%
  filter(n_include_values > 1)
if (nrow(include_check))
  stop("Include must be identical across every day of a slice. Fix: ",
       paste(paste(include_check$Group, include_check$Slice, sep = " | "), collapse = ", "))
if (any(!dendra_all$Include & !nzchar(dendra_all$ExclusionReason)))
  stop("Every excluded slice needs a prespecified ExclusionReason.")

excluded_records <- dendra_all %>% filter(!Include)
write.csv(excluded_records, file.path(TABLES, "00_excluded_records.csv"), row.names = FALSE)
dendra <- dendra_all %>% filter(Include)
if (!nrow(dendra)) stop("All rows are excluded; no analysis can be run.")

## Drop excluded backgrounds (see EXCLUDE_BACKGROUNDS at the top of this file).
## Done here, before any modelling, so an excluded background cannot reach a
## figure, table or p-value.
if (length(EXCLUDE_BACKGROUNDS)) {
  .pat <- paste(EXCLUDE_BACKGROUNDS, collapse = "|")
  .drop <- grepl(.pat, dendra$Group, ignore.case = TRUE)
  if (any(.drop)) {
    .dropped_groups <- sort(unique(as.character(dendra$Group[.drop])))
    warning("Excluding background(s) [", paste(EXCLUDE_BACKGROUNDS, collapse = ", "),
            "]: dropped ", length(.dropped_groups), " group(s) -- ",
            paste(.dropped_groups, collapse = ", "),
            ". This comes from 'Exclude backgrounds' in the Analysis Settings sheet. ",
            "Clear that setting if the drop was not intended.", call. = FALSE,
            immediate. = TRUE)
    dendra <- dendra[!.drop, , drop = FALSE]
    if (!nrow(dendra)) stop("EXCLUDE_BACKGROUNDS removed every group.")
  }
}
EXCLUDED_BACKGROUND_NOTE <- if (length(EXCLUDE_BACKGROUNDS))
  paste0("Excluded from this analysis: ", paste(EXCLUDE_BACKGROUNDS, collapse = ", "),
         " (archived separately)") else "None"

if (any(duplicated(dendra[c("Group", "Slice", "Day")]))) {
  bad <- dendra[duplicated(dendra[c("Group", "Slice", "Day")]) |
                  duplicated(dendra[c("Group", "Slice", "Day")], fromLast = TRUE),
                c("Group", "Slice", "Day")]
  stop("Duplicate included Group + Slice ID + Day rows were found:\n",
       paste(apply(unique(bad), 1, paste, collapse = " | "), collapse = "\n"))
}

prep_check <- dendra %>% group_by(Group, Slice) %>%
  summarise(n_prep = n_distinct(PreparationID[nzchar(PreparationID)]), .groups = "drop") %>%
  filter(n_prep > 1)
if (nrow(prep_check))
  stop("A slice cannot belong to multiple Preparation IDs. Fix: ",
       paste(paste(prep_check$Group, prep_check$Slice, sep = " | "), collapse = ", "))
dendra <- dendra %>% group_by(Group, Slice) %>%
  mutate(PreparationID = { p <- unique(PreparationID[nzchar(PreparationID)]);
                           if (length(p)) p[1] else "" }) %>% ungroup()

.levels <- unique(dendra$Group)
REFERENCE_GROUP_AUTO <- FALSE
if (!is.na(REFERENCE_GROUP) && REFERENCE_GROUP %in% .levels) {
  .levels <- c(REFERENCE_GROUP, setdiff(.levels, REFERENCE_GROUP))
} else {
  ## The configured reference name is absent -- the normal case once group
  ## names get more descriptive ("Control" becomes "WT-tau GFAP Control").
  .ctrl <- .levels[grepl("control|untreated|vehicle|wild[ -]?type", .levels,
                         ignore.case = TRUE)]
  .chosen <- if (length(.ctrl)) .ctrl[1] else .levels[1]
  REFERENCE_GROUP_AUTO <- TRUE
  message("NOTE: reference group '", REFERENCE_GROUP, "' is not in this sheet.\n",
          "      Using '", .chosen, "' as the reference for all contrasts.\n",
          "      Set 'Reference group' on the Analysis Settings sheet to choose ",
          "explicitly.")
  if (!length(.ctrl))
    warning("No group name looks like a control, so '", .chosen, "' was used as ",
            "the reference only because it is listed first. Set 'Reference group' ",
            "on the Analysis Settings sheet.")
  REFERENCE_GROUP <- .chosen
  .levels <- c(.chosen, setdiff(.levels, .chosen))
}

## If the groups form a Background x Control/HSPB1 grid, order the levels so
.bg_of <- function(g) trimws(sub("(?i)\\s*GFAP.*$", "", g, perl = TRUE))
.cons_of <- function(g) ifelse(grepl("control", g, ignore.case = TRUE), "Control",
                        ifelse(grepl("hspb1", g, ignore.case = TRUE), "HSPB1", NA))
.bg <- .bg_of(.levels); .cons <- .cons_of(.levels)
FACTORIAL <- !anyNA(.cons) && length(unique(.cons)) == 2
if (FACTORIAL) {
  .bg_order <- unique(.bg)                       # reference's background first
  .levels <- .levels[order(match(.bg, .bg_order),
                           match(.cons, c("Control", "HSPB1")))]
}

dendra <- dendra %>%
  mutate(Group = factor(Group, levels = .levels),
         SliceID = paste(Group, Slice, sep = " | ")) %>%
  arrange(Group, Slice, Day)

GROUPS <- levels(dendra$Group)
N_GROUP <- length(GROUPS)
TWO_GROUP <- N_GROUP == 2
if (N_GROUP < 2) stop("At least two groups are required.")

## Factorial lookups + the construct palette shared by every figure script:
## Control is ALWAYS blue and HSPB1 ALWAYS vermillion, in every figure.
GROUP_BG <- setNames(.bg_of(GROUPS), GROUPS)
GROUP_CONS <- setNames(.cons_of(GROUPS), GROUPS)
BACKGROUNDS <- unique(unname(GROUP_BG))
CONS_PAL <- c(Control = "#0072B2", HSPB1 = "#D55E00")
BG_SHAPES <- setNames(c(21, 24, 22, 23, 25)[seq_along(BACKGROUNDS)], BACKGROUNDS)

prep_by_slice <- dendra %>% distinct(Group, SliceID, PreparationID)
PREP_COMPLETE <- all(nzchar(prep_by_slice$PreparationID))
N_PREP <- if (PREP_COMPLETE) n_distinct(prep_by_slice$PreparationID) else 0L
prep_groups <- if (PREP_COMPLETE) prep_by_slice %>%
  group_by(PreparationID) %>% summarise(n_groups = n_distinct(Group), .groups = "drop") else NULL
ADJUST_FOR_PREP <- PREP_COMPLETE && N_PREP >= 2 && any(prep_groups$n_groups > 1)
PREP_STATUS <- if (!PREP_COMPLETE) {
  "Not recorded. Slices are treated as independent; fill in PreparationID before publication."
} else if (ADJUST_FOR_PREP) {
  paste0("Recorded. PreparationID is included as a fixed blocking term (", N_PREP, " levels).")
} else {
  "Recorded, but confounded with group, so a preparation-adjusted group effect cannot be estimated."
}
if (PREP_COMPLETE) dendra$PreparationID <- factor(dendra$PreparationID)

LMM_FIXED_TERMS <- c("Day * Group", if (ADJUST_FOR_PREP) "PreparationID")
lmm_formula <- function(response, random_term)
  reformulate(c(LMM_FIXED_TERMS, random_term), response = response)

## ---- 5. ANALYSIS WINDOWS ----------------------------------------------------
rfp_model <- dendra %>% filter(Day > BASELINE_RFP)
gfp_model <- dendra %>% filter(Day > BASELINE_GFP)
decay_fit <- dendra %>% filter(Day > BASELINE_RFP) %>%
  mutate(t = Day - BASELINE_RFP, frac = 1 + RFP_dFF0)
decay_auc <- dendra %>% filter(Day >= BASELINE_RFP) %>%
  mutate(t = Day - BASELINE_RFP, frac = 1 + RFP_dFF0)
decay <- decay_fit  # backward-compatible name; fitting scripts use decay_fit explicitly

## ---- 6. FIGURE THEME AND PALETTE -------------------------------------------
OKABE_ITO <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#000000")
if (N_GROUP > length(OKABE_ITO)) stop("At most ", length(OKABE_ITO), " groups are supported.")
PAL <- setNames(OKABE_ITO[seq_len(N_GROUP)], GROUPS)
SHAPES <- setNames(c(21, 22, 24, 23, 25, 21, 22)[seq_len(N_GROUP)], GROUPS)
LINES <- setNames(c("solid", "solid", "22", "42", "1343", "solid", "22")[seq_len(N_GROUP)], GROUPS)
W1 <- 89; W2 <- 183

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.3, colour = "grey20"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey20"),
          axis.text = element_text(colour = "grey20", size = rel(0.95)),
          axis.title = element_text(colour = "black"), legend.key.size = unit(3.2, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.title = element_blank(),
          legend.text = element_text(size = rel(0.95)), legend.position = "bottom",
          legend.background = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0, size = rel(1)),
          plot.title = element_text(face = "bold", size = rel(1), hjust = 0),
          plot.subtitle = element_text(size = rel(0.9), colour = "grey35", margin = margin(b = 4)),
          plot.caption = element_text(size = rel(0.85), colour = "grey35", hjust = 0),
          plot.tag = element_text(face = "bold", size = rel(1.25)), plot.margin = margin(2, 2, 2, 2))
}
## Day-axis tick labels. Two jobs: drop the trailing ".0" so integer days read
## as "3" not "3.0", and make the photoconversion timepoint self-explanatory --
## the fractional RFP baseline (1.1) exists only to order the post-conversion
## scan after the pre-conversion one within the same imaging day; it is not a
## tenth of a day, so the axis shows the calendar day and the second line
## carries the pre/post distinction.
day_labels <- function(x) {
  is_conv <- abs(x - BASELINE_RFP) < 1e-8
  d <- ifelse(is_conv, floor(x), x)
  lab <- ifelse(d == round(d), format(round(d)), format(d))
  ifelse(is_conv, paste0(lab, "\nconverted"), lab)
}
## the same day written for prose: "day 1", not "day 1.1"
BASELINE_RFP_DISPLAY <- format(floor(BASELINE_RFP))

x_group_scale <- function(narrow_panel = FALSE) {
  longest <- max(nchar(GROUPS))
  if (narrow_panel && (N_GROUP > 2 || longest > 8))
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.8)))
  else if (N_GROUP > 2 || longest > 8)
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.85)))
  else theme()
}
fig_width <- function(min_mm = W1, per_group_mm = 16) min(W2, max(min_mm, per_group_mm * N_GROUP + 34))
wrap_text <- function(x, width_mm) paste(strwrap(x, width = max(40, round(width_mm / 1.7))), collapse = "\n")
save_fig <- function(plot, name, width_mm, height_mm) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(FIGURES, paste0(name, ".pdf")), plot, width = width_mm,
         height = height_mm, units = "mm", device = pdf_device)
  ggsave(file.path(FIGURES, paste0(name, ".png")), plot, width = width_mm,
         height = height_mm, units = "mm", dpi = 300)
  message("  figure: ", name, ".pdf / .png")
}
save_table <- function(df, name) {
  write.csv(df, file.path(TABLES, paste0(name, ".csv")), row.names = FALSE)
  message("  table:  ", name, ".csv")
}

## ---- 7. PROVENANCE ----------------------------------------------------------
describe_run <- function() {
  n <- dendra %>% distinct(Group, Slice) %>% count(Group)
  cat("Workbook      : ", WORKBOOK, "\n", sep = "")
  cat("Workbook MD5  : ", unname(tools::md5sum(WORKBOOK)), "\n", sep = "")
  cat("Analysed      : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("R version     : ", R.version.string, "\n", sep = "")
  cat("Primary test  : ", PRIMARY_ANALYSIS, "\n", sep = "")
  cat("Groups        : ", paste(sprintf("%s n=%d", n$Group, n$n), collapse = ",  "), "\n", sep = "")
  cat("Reference     : ", GROUPS[1],
      if (REFERENCE_GROUP_AUTO)
        "   <- CHOSEN AUTOMATICALLY; set 'Reference group' on Analysis Settings to override"
      else "", "\n", sep = "")
  cat("Preparation  : ", PREP_STATUS, "\n", sep = "")
  cat("Baselines     : GFP day ", BASELINE_GFP, "; RFP day ", BASELINE_RFP, "\n", sep = "")
  cat("Rows          : ", nrow(dendra), " included; ", nrow(excluded_records), " excluded\n", sep = "")
  cat("Results       : ", RESULTS, "\n", sep = "")
}
