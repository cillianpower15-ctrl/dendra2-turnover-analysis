## 01_statistics.R  --  all statistical tests. Run this first.
##   Rscript 01_statistics.R

source("00_setup.R")

.log <- file.path(RESULTS, "analysis_log.txt")
con <- file(.log, open = "wt"); sink(con, split = TRUE)

hdr <- function(x) cat("\n\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
sub <- function(x) cat("\n-- ", x, " ", strrep("-", max(0, 72 - nchar(x))), "\n", sep = "")
show_df <- function(d, digits = 4) print(as.data.frame(d), digits = digits, row.names = FALSE)

hdr("DENDRA2 TURNOVER ANALYSIS")
describe_run()

hdr("0. INPUT QUALITY CONTROL")
expected_days <- sort(unique(dendra$Day))
slice_qc <- dendra %>% group_by(Group, Slice, SliceID) %>%
  summarise(
    n_rows = n(),
    days_present = paste(sort(unique(Day)), collapse = ", "),
    missing_days = paste(setdiff(expected_days, Day), collapse = ", "),
    missing_measurements = sum(is.na(GFP_Raw) | is.na(GFP_dFF0) |
                               is.na(RFP_Raw) | is.na(RFP_dFF0)),
    GFP_at_baseline = {x <- GFP_dFF0[Day == BASELINE_GFP]; if (length(x)) x[1] else NA_real_},
    RFP_at_baseline = {x <- RFP_dFF0[Day == BASELINE_RFP]; if (length(x)) x[1] else NA_real_},
    retained_z_planes = paste(sort(unique(RetainedZPlanes[!is.na(RetainedZPlanes)])), collapse = ", "),
    registration_qc = paste(sort(unique(RegistrationQC)), collapse = ", "),
    minimum_SBR = if (all(is.na(SignalBackgroundRatio))) NA_real_ else min(SignalBackgroundRatio, na.rm = TRUE),
    status = case_when(
      nzchar(missing_days) ~ "WARNING: missing timepoint(s)",
      missing_measurements > 0 ~ "WARNING: missing measurement(s)",
      is.na(GFP_at_baseline) | is.na(RFP_at_baseline) ~ "WARNING: baseline missing",
      abs(GFP_at_baseline) > 1e-6 | abs(RFP_at_baseline) > 1e-6 ~ "WARNING: baseline dF/F0 is not zero",
      any(tolower(RegistrationQC) == "fail") ~ "WARNING: registration QC failed",
      !is.na(minimum_SBR) && minimum_SBR < MIN_SBR ~ "WARNING: signal/background below threshold",
      TRUE ~ "PASS"),
    .groups = "drop")

at_day <- function(x, day, target) { y <- x[day == target]; if (length(y)) y[1] else NA_real_ }
photoconversion_qc <- dendra %>% group_by(Group, Slice, SliceID) %>%
  summarise(
    GFP_pre = at_day(GFP_Raw, Day, BASELINE_GFP),
    GFP_post = at_day(GFP_Raw, Day, BASELINE_RFP),
    RFP_pre = at_day(RFP_Raw, Day, BASELINE_GFP),
    RFP_post = at_day(RFP_Raw, Day, BASELINE_RFP),
    GFP_loss_fraction = ifelse(is.finite(GFP_pre) && GFP_pre > 0, 1 - GFP_post / GFP_pre, NA_real_),
    RFP_fold_increase = ifelse(is.finite(RFP_pre) && RFP_pre > 0, RFP_post / RFP_pre, NA_real_),
    photoconversion_status = case_when(
      is.na(GFP_loss_fraction) | is.na(RFP_fold_increase) ~ "NOT ASSESSABLE: raw pre/post signal missing",
      GFP_loss_fraction < MIN_GFP_LOSS ~ "REVIEW: GFP loss below threshold",
      RFP_fold_increase < MIN_RFP_FOLD ~ "REVIEW: RFP fold increase below threshold",
      TRUE ~ "PASS"),
    .groups = "drop")
slice_qc <- left_join(slice_qc, photoconversion_qc, by = c("Group", "Slice", "SliceID"))
show_df(slice_qc)
save_table(slice_qc, "00_input_QC")
save_table(photoconversion_qc, "00_photoconversion_QC")

## Things a reviewer will ask about that the workbook cannot answer by itself.
## Wording is deliberately plain: these are questions to answer, not audit codes.
readiness <- data.frame(
  item = c("How many animals or preparations?",
           "Are excluded slices explained?",
           "Z-planes kept per slice",
           "Registration checked",
           "Signal-to-background ratio",
           "Photoconversion worked",
           "Primary analysis fixed in advance"),
  status = c(
    if (ADJUST_FOR_PREP) "Recorded" else "Missing",
    if (all(nzchar(excluded_records$ExclusionReason)) || !nrow(excluded_records)) "Yes" else "Reasons missing",
    if (all(!is.na(dendra$RetainedZPlanes))) "Recorded" else "Not recorded",
    if (all(tolower(dendra$RegistrationQC) %in% c("pass"))) "Passed" else "Not recorded",
    if (all(!is.na(dendra$SignalBackgroundRatio)) && all(dendra$SignalBackgroundRatio >= MIN_SBR)) "Passed" else "Not recorded",
    if (all(photoconversion_qc$photoconversion_status == "PASS")) "Passed" else "Check",
    if (identical(PRIMARY_ANALYSIS, "Constrained nonlinear RFP rate model")) "Yes" else "Check"),
  detail = c(PREP_STATUS,
             paste(n_distinct(excluded_records$Group, excluded_records$Slice), "excluded slices of",
                   n_distinct(dendra_all$Group, dendra_all$Slice)),
             "Enter the common cropped Z-plane count for every day; do not infer it from plot axes.",
             "Allowed workbook values: Pass, Review, Fail, Not recorded.",
             paste("Prespecified minimum:", MIN_SBR),
             paste("Thresholds: GFP loss >=", MIN_GFP_LOSS, "; RFP fold increase >=", MIN_RFP_FOLD),
             PRIMARY_ANALYSIS), stringsAsFactors = FALSE)
save_table(readiness, "00_publication_readiness")
cat("\nBEFORE REPORTING THESE RESULTS\n"); show_df(readiness)
if (any(!readiness$status %in% c("Recorded", "Yes", "Passed")))
  cat("\nThe statistics are sound, but the items above are not yet recorded.\n",
      "Fill them in or explain them in the manuscript before submitting.\n", sep = "")

n_by_group <- dendra %>% distinct(Group, Slice) %>% count(Group, name = "n_slices")
if (any(n_by_group$n_slices < 2))
  stop("Every group needs at least two independent slices. Check 00_input_QC.csv.")
if (any(decay_fit$frac <= 0, na.rm = TRUE))
  stop("At least one post-conversion RFP fraction (1 + RFP_dFF0) is zero or negative. ",
       "The decay model cannot fit this value. Check 00_input_QC.csv and the source data.")

statistics_manifest <- data.frame(
  order = 1:9,
  analysis = c(
    "Descriptive statistics",
    "Linear mixed-effects model: RFP",
    "Linear mixed-effects model: GFP",
    "Adjusted pairwise contrasts",
    "Per-slice nonlinear RFP half-life",
    "Per-slice AUC",
    "Per-slice group comparisons",
    "Sensitivity and leave-one-slice-out checks",
    "Nonlinear mixed-effects RFP decay model"),
  question = c(
    "Mean, SD and SEM by group and day",
    "Do RFP level and linear rate differ by group?",
    "Do GFP level and linear rate differ by group?",
    "Which groups differ in level or slope?",
    "What is each slice's decay rate and half-life?",
    "What is the model-free integrated response per slice?",
    "Do per-slice half-life/AUC values differ between groups?",
    "Is the conclusion driven by model choice or one slice?",
    "Does the exponential RFP rate constant K differ by group?"),
  method = c(
    "n, mean, SD, SEM",
    "dF/F0 ~ Day * Group + prespecified estimable random slope; optional preparation block; REML",
    "dF/F0 ~ Day * Group + prespecified estimable random slope; optional preparation block; REML",
    "emmeans/emtrends; Tukey adjustment",
    "Constrained one-phase decay with plateau; structural t=0 excluded from likelihood",
    "Trapezoidal integration from the channel-specific structural baseline",
    "Welch and rank tests; Holm adjustment across the three primary metrics",
    "Alternative random effects, baseline rules and leave-one-slice-out",
    "nlme on log(K) and logit(P); ML likelihood-ratio test; Holm rate-ratio tests"),
  role = c("Descriptive", "Secondary linear trajectory", "Secondary linear trajectory",
           "Adjusted follow-up", "Supporting", "Supporting", "Supporting",
           "Sensitivity", "Primary RFP turnover-rate"),
  stringsAsFactors = FALSE)
save_table(statistics_manifest, "00_statistics_performed")

## 1. DESCRIPTIVE SUMMARY
hdr("1. DESCRIPTIVE SUMMARY  (mean, SD, SEM per group per day)")

describe <- function(df, col) {
  df %>% group_by(Group, Day) %>%
    summarise(n = sum(!is.na(.data[[col]])),
              mean = mean(.data[[col]], na.rm = TRUE),
              sd   = sd(.data[[col]], na.rm = TRUE),
              sem  = sd / sqrt(n), .groups = "drop") %>%
    mutate(channel = col, .before = 1)
}
desc <- bind_rows(describe(dendra, "RFP_dFF0"), describe(dendra, "GFP_dFF0"))
show_df(desc)
save_table(desc, "01_descriptive_summary")

hdr("2. PRIMARY TEST -- LINEAR MIXED-EFFECTS MODELS")
cat("
Model:   dF/F0 ~ Day * Group + <random structure>
Fitted by REML; p-values from lmerTest using Satterthwaite degrees of freedom.
The random effect makes the SLICE the unit of replication, not the observation.

  Group      = do the groups sit at different overall levels?
  Day        = the rate of change common to all groups
  Day:Group  = do the groups change at DIFFERENT rates?
               (the term that tests 'this condition alters turnover')

With more than two groups the Day:Group and Group rows of the ANOVA are
OMNIBUS tests -- 'do the groups differ at all'. Which groups differ from which
is answered by the pairwise contrasts in section 2c.

The random-slope structure is PRESPECIFIED: (Day | SliceID) is used whenever
it converges and is non-singular. It falls back to (1 | SliceID) only when the
data cannot estimate a between-slice slope variance. The likelihood-ratio
comparison is diagnostic and is not used as a p-value-driven selection rule.

Preparation block: ", if (ADJUST_FOR_PREP) "PreparationID is included as a fixed blocking term."
else PREP_STATUS, "

Timepoints: RFP days above ", BASELINE_RFP, ", GFP days above ", BASELINE_GFP,
". Baseline days are excluded because dF/F0 is zero there by construction --
see README note 2.
")

fit_lmm <- function(response, data, label) {
  f_slope <- lmm_formula(response, "(Day | SliceID)")
  f_int   <- lmm_formula(response, "(1 | SliceID)")
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

  m_slope <- suppressMessages(suppressWarnings(
    try(lmer(f_slope, data = data, control = ctrl), silent = TRUE)))
  m_int <- suppressMessages(lmer(f_int, data = data, control = ctrl))

  conv_bad <- function(m) {
    msg <- m@optinfo$conv$lme4$messages
    length(msg) && any(grepl("failed to converge", msg, ignore.case = TRUE))
  }
  ok <- !inherits(m_slope, "try-error") && !isSingular(m_slope) && !conv_bad(m_slope)
  lrt <- NULL; keep_slope <- ok
  if (ok) {
    lrt <- suppressMessages(suppressWarnings(anova(
      lmer(f_int,   data = data, REML = FALSE, control = ctrl),
      lmer(f_slope, data = data, REML = FALSE, control = ctrl))))
  }

  sub(paste0(label, " -- choice of random-effects structure"))
  if (inherits(m_slope, "try-error")) cat("  (Day | SliceID) failed to fit.\n")
  else cat(sprintf("  (Day | SliceID)  singular = %s\n", isSingular(m_slope)))
  if (!is.null(lrt)) {
    cat(sprintf("  likelihood-ratio vs (1 | SliceID): chisq = %.2f, df = %d, p = %.3g\n",
                lrt$Chisq[2], lrt$Df[2], lrt$`Pr(>Chisq)`[2]))
    cat(sprintf("  AIC: intercept-only %.1f, +slope %.1f\n", AIC(m_int), AIC(m_slope)))
  }
  cat(if (keep_slope)
        "  -> USING the prespecified (Day | SliceID) structure.\n"
      else "  -> FALLBACK to (1 | SliceID): the random slope was not estimable.\n")
  m <- if (keep_slope) m_slope else m_int
  attr(m, "structure") <- if (keep_slope) "(Day | SliceID)" else "(1 | SliceID)"
  m
}

report_lmm <- function(model, label) {
  sub(paste0(label, " -- fixed effects"))
  co <- as.data.frame(summary(model)$coefficients)
  ci <- suppressMessages(confint(model, method = "Wald", parm = rownames(co)))
  out <- data.frame(model = label, term = rownames(co), estimate = co[[1]],
                    std_error = co[[2]], ci_low = ci[, 1], ci_high = ci[, 2],
                    df = co[[3]], t = co[[4]], p = co[[5]], row.names = NULL)
  show_df(out)

  sub(paste0(label, " -- Type III ANOVA (Satterthwaite)"))
  av <- anova(model); print(av, digits = 4)

  sub(paste0(label, " -- random effects and fit"))
  cat("  structure = ", attr(model, "structure"), "\n", sep = "")
  print(VarCorr(model), digits = 3)
  cat(sprintf("  observations = %d, slices = %d\n", nobs(model), ngrps(model)))
  cat(if (isSingular(model))
    "  *** SINGULAR FIT -- random effect not supported by the data ***\n"
    else "  Fit is non-singular.\n")

  list(coefs = out, anova = data.frame(model = label, term = rownames(av), av, row.names = NULL))
}

m_rfp <- fit_lmm("RFP_dFF0", rfp_model, "RFP")
r_rfp <- report_lmm(m_rfp, "RFP")
m_gfp <- fit_lmm("GFP_dFF0", gfp_model, "GFP")
r_gfp <- report_lmm(m_gfp, "GFP")

save_table(bind_rows(r_rfp$coefs, r_gfp$coefs), "02_mixed_model_coefficients")
save_table(bind_rows(r_rfp$anova, r_gfp$anova), "03_mixed_model_ANOVA")

## ---- 2c. pairwise contrasts ------------------------------------------------
hdr("2c. PAIRWISE CONTRASTS BETWEEN GROUPS")
cat("
Estimated by emmeans from the fitted mixed models, Tukey-adjusted for the
number of pairs. Two questions, and they are different:

  RATE      -- does the slope over time differ between two groups?
               This is the turnover question.
  LEVEL     -- averaged over the timecourse, does one group sit higher or
               lower than another?

With ", N_GROUP, " groups there are ", choose(N_GROUP, 2), " pairs, so the
adjustment matters. Unadjusted p-values are given alongside for transparency.
")

contrasts_for <- function(model, label) {
  tr <- suppressMessages(emtrends(model, pairwise ~ Group, var = "Day", adjust = "tukey"))
  lv <- suppressMessages(emmeans(model, pairwise ~ Group, adjust = "tukey"))
  tr_u <- suppressMessages(summary(emtrends(model, pairwise ~ Group, var = "Day",
                                            adjust = "none")$contrasts))
  lv_u <- suppressMessages(summary(emmeans(model, pairwise ~ Group,
                                           adjust = "none")$contrasts))
  sub(paste0(label, " -- rate of change per group (slope per day)"))
  show_df(as.data.frame(summary(tr$emtrends)))
  sub(paste0(label, " -- RATE differences between groups"))
  a <- as.data.frame(summary(tr$contrasts))
  a$p_unadjusted <- tr_u$p.value
  show_df(a)
  sub(paste0(label, " -- LEVEL differences between groups"))
  b <- as.data.frame(summary(lv$contrasts))
  b$p_unadjusted <- lv_u$p.value
  show_df(b)
  bind_rows(data.frame(model = label, comparison = "rate (slope/day)", a),
            data.frame(model = label, comparison = "level (mean dF/F0)", b))
}
ctr <- bind_rows(contrasts_for(m_rfp, "RFP"), contrasts_for(m_gfp, "GFP"))
save_table(ctr, "04_pairwise_contrasts")

## ---- 2d. sensitivity to the timepoint-exclusion rule -----------------------
hdr("2d. SENSITIVITY -- effect of including the by-construction baseline zeros")
cat("
NOT the primary results. Here so the effect of the timepoint-exclusion
decision is visible and can be defended to a reviewer.
")
sens <- function(response, data, label, random_term = "(1 | SliceID)") {
  fml <- lmm_formula(response, random_term)
  m <- suppressMessages(suppressWarnings(try(lmer(fml, data = data,
        control = lmerControl(optimizer = "bobyqa")), silent = TRUE)))
  if (inherits(m, "try-error"))
    return(data.frame(variant = label, term = "MODEL_FAILED", F_value = NA_real_,
                      p = NA_real_, singular = NA))
  av <- anova(m)
  data.frame(variant = label, term = rownames(av),
             F_value = av$`F value`, p = av$`Pr(>F)`,
             singular = isSingular(m), row.names = NULL)
}
sens_tbl <- bind_rows(
  sens("RFP_dFF0", rfp_model, "RFP PRIMARY (baseline excluded)"),
  sens("RFP_dFF0", filter(dendra, Day >= BASELINE_RFP), "RFP + structural zero"),
  sens("RFP_dFF0", dendra, "RFP + pre-conversion day"),
  sens("GFP_dFF0", gfp_model, "GFP PRIMARY (baseline excluded)"),
  sens("GFP_dFF0", dendra, "GFP + structural zero"))
show_df(sens_tbl)

sub("Effect of the random-effects structure (timepoints held at PRIMARY)")
cat("
Watch the Day:Group row. Adding the random slope, where the data support it,
can change this term by orders of magnitude -- the simpler model counts
between-slice variation in rate as evidence for a group difference.
")
rand_tbl <- bind_rows(
  sens("RFP_dFF0", rfp_model, "RFP  (1|SliceID)", "(1 | SliceID)"),
  sens("RFP_dFF0", rfp_model, "RFP  (Day|SliceID)", "(Day | SliceID)"),
  sens("GFP_dFF0", gfp_model, "GFP  (1|SliceID)", "(1 | SliceID)"),
  sens("GFP_dFF0", gfp_model, "GFP  (Day|SliceID)", "(Day | SliceID)"))
show_df(rand_tbl)
cat("  (singular = TRUE means that structure is NOT supported by the data)\n")
save_table(bind_rows(sens_tbl, rand_tbl), "05_model_sensitivity")

hdr("3. PER-REPLICATE RFP HALF-LIFE  (one value per slice)")
cat("
One-phase decay fitted by nonlinear least squares to the red fraction
  frac = 1 + RFP_dFF0   against   t = Day - ", BASELINE_RFP, "

  PRIMARY   frac = (1 - P) * exp(-K*t) + P     free plateau, constrained P >= 0
  SECONDARY frac = exp(-K*t)                   no plateau
  SECONDARY log-linear  lm(log(frac) ~ Day)    the Excel workbook's method

Half-life = ln(2)/K. These estimate different quantities -- with a free
plateau, ln(2)/K is the half-time of approach to the plateau, not the
half-life of the whole pool. All three are reported because the point estimate
is sensitive to the choice while the group ordering usually is not.

NOTE: 04_nonlinear_mixed_model.R answers this same question better, by fitting
all slices at once instead of compressing each to a single number first.
")

fit_decay <- function(t, frac, model) {
  f <- switch(model,
    plateau = try(nls(frac ~ (1 - P) * exp(-K * t) + P,
                      start = list(K = 0.05, P = 0.05), algorithm = "port",
                      lower = c(K = 1e-4, P = 0), upper = c(K = 5, P = 0.99),
                      control = nls.control(maxiter = 500)), silent = TRUE),
    pure    = try(nls(frac ~ exp(-K * t), start = list(K = 0.05),
                      control = nls.control(maxiter = 500)), silent = TRUE))
  if (inherits(f, "try-error")) return(c(K = NA, plateau = NA, half_life = NA, R2 = NA))
  K <- coef(f)["K"]; P <- if ("P" %in% names(coef(f))) coef(f)["P"] else 0
  c(K = unname(K), plateau = unname(P), half_life = unname(log(2) / K),
    R2 = 1 - sum(residuals(f)^2) / sum((frac - mean(frac))^2))
}

hl <- decay_fit %>% group_by(Group, Slice, SliceID) %>%
  summarise(a  = list(fit_decay(t, frac, "plateau")),
            b  = list(fit_decay(t, frac, "pure")),
            ll = list({ m <- lm(log(pmax(frac, 1e-3)) ~ Day); k <- -coef(m)[2]
                        c(half_life = unname(log(2)/k), R2 = summary(m)$r.squared) }),
            .groups = "drop") %>%
  mutate(K = map_dbl(a, "K"), plateau = map_dbl(a, "plateau"),
         half_life = map_dbl(a, "half_life"), R2 = map_dbl(a, "R2"),
         half_life_noplateau = map_dbl(b, "half_life"),
         half_life_loglinear = map_dbl(ll, "half_life")) %>%
  select(-a, -b, -ll)

sub("Primary fit (free plateau, P >= 0)")
show_df(hl %>% select(Group, Slice, K, plateau, half_life, R2))

sub("Fit quality check")
bad <- hl %>% filter(is.na(R2) | R2 < MIN_DECAY_R2)
if (nrow(bad) == 0) {
  cat("  All ", nrow(hl), " fits have R2 >= ", MIN_DECAY_R2, ". Minimum R2 = ",
      sprintf("%.3f", min(hl$R2, na.rm = TRUE)), "\n", sep = "")
} else {
  cat("  *** Fits below the R2 = ", MIN_DECAY_R2,
      " threshold, or failed: ***\n", sep = ""); show_df(bad)
}
save_table(hl, "06_per_replicate_halflife")

## 4. PER-REPLICATE AUC
hdr("4. PER-REPLICATE AREA UNDER THE CURVE  (trapezoidal)")
cat("
Model-free: no curve is fitted, so nothing is assumed about the decay shape.
A more negative AUC means signal fell further and faster over the chase.
Each channel starts at its own structural baseline (GFP day ", BASELINE_GFP,
"; RFP day ", BASELINE_RFP, "). If an intermediate day is missing, trapezoidal
integration bridges that gap with a straight line and the missing-day warning
remains visible in 00_input_QC.csv.
")
auc_trapz <- function(x, y) { o <- order(x); x <- x[o]; y <- y[o]
                              sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) }
auc <- dendra %>% group_by(Group, Slice, SliceID) %>%
  summarise(RFP_AUC = auc_trapz(Day[Day >= BASELINE_RFP], RFP_dFF0[Day >= BASELINE_RFP]),
            GFP_AUC = auc_trapz(Day[Day >= BASELINE_GFP], GFP_dFF0[Day >= BASELINE_GFP]),
            .groups = "drop")
show_df(auc)
save_table(auc, "07_per_replicate_AUC")

hdr("5. GROUP COMPARISONS ON THE PER-REPLICATE METRICS")
cat("
Each slice contributes one number, so n = number of slices per group.

With TWO groups: Welch t-test (does not assume equal variances) and
Mann-Whitney (does not assume normality), plus Hedges' g and Cliff's delta.
With MORE than two: Welch one-way ANOVA and Kruskal-Wallis as omnibus tests,
then all pairwise comparisons, Holm-adjusted within each metric.

With small n, normality cannot be established, so BOTH the parametric and the
rank-based result are reported and the more conservative should be quoted.
Effect sizes matter more than p-values at this sample size.
")

hedges_g <- function(a, b) {
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2))
  ((mean(a) - mean(b)) / sp) * (1 - 3 / (4 * (na + nb) - 9))
}
cliffs_delta <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))

group_desc <- function(tbl, metric, label) {
  tbl %>% group_by(Group) %>%
    summarise(metric = label, n = sum(!is.na(.data[[metric]])),
              mean = mean(.data[[metric]], na.rm = TRUE),
              sd = sd(.data[[metric]], na.rm = TRUE),
              median = median(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    select(metric, everything())
}

compare_metric <- function(tbl, metric, label) {
  d <- tbl %>% select(Group, value = all_of(metric)) %>% filter(!is.na(value))
  sub(paste0(label, " -- per-group summary"))
  show_df(group_desc(tbl, metric, label))

  if (TWO_GROUP) {
    a <- d$value[d$Group == GROUPS[1]]; b <- d$value[d$Group == GROUPS[2]]
    tt <- t.test(a, b); mw <- suppressWarnings(wilcox.test(a, b))
    res <- data.frame(metric = label,
      contrast = paste(GROUPS[1], "-", GROUPS[2]),
      difference = mean(a) - mean(b),
      ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
      welch_p = tt$p.value, rank_p = mw$p.value,
      hedges_g = hedges_g(a, b), cliffs_delta = cliffs_delta(a, b), row.names = NULL)
    sub(paste0(label, " -- group comparison"))
  } else {
    om_w <- oneway.test(value ~ Group, data = d, var.equal = FALSE)
    om_k <- kruskal.test(value ~ Group, data = d)
    sub(paste0(label, " -- omnibus tests"))
    cat(sprintf("  Welch one-way ANOVA : F = %.3f, df = %.1f, %.1f, p = %.4g\n",
                om_w$statistic, om_w$parameter[1], om_w$parameter[2], om_w$p.value))
    cat(sprintf("  Kruskal-Wallis      : chisq = %.3f, df = %d, p = %.4g\n",
                om_k$statistic, om_k$parameter, om_k$p.value))

    pt <- pairwise.t.test(d$value, d$Group, p.adjust.method = "holm", pool.sd = FALSE)
    pw <- suppressWarnings(pairwise.wilcox.test(d$value, d$Group, p.adjust.method = "holm"))
    pairs <- combn(GROUPS, 2, simplify = FALSE)
    res <- do.call(rbind, lapply(pairs, function(p) {
      a <- d$value[d$Group == p[1]]; b <- d$value[d$Group == p[2]]
      tt <- t.test(a, b)
      data.frame(metric = label, contrast = paste(p[1], "-", p[2]),
                 difference = mean(a) - mean(b),
                 ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
                 welch_p = tryCatch(pt$p.value[p[2], p[1]], error = function(e) NA),
                 rank_p  = tryCatch(pw$p.value[p[2], p[1]], error = function(e) NA),
                 hedges_g = hedges_g(a, b), cliffs_delta = cliffs_delta(a, b),
                 row.names = NULL)
    }))
    res$omnibus_welch_p <- om_w$p.value
    res$omnibus_kruskal_p <- om_k$p.value
    sub(paste0(label, " -- pairwise (Holm-adjusted within this metric)"))
  }
  show_df(res)
  res
}

primary <- bind_rows(
  compare_metric(hl,  "half_life", "RFP half-life (days)"),
  compare_metric(auc, "RFP_AUC",   "RFP AUC"),
  compare_metric(auc, "GFP_AUC",   "GFP AUC"))

## With two groups there is one comparison per metric. Control the family of
## three primary per-replicate metrics, while retaining the raw values.
if (TWO_GROUP) {
  primary$welch_p_unadjusted <- primary$welch_p
  primary$rank_p_unadjusted  <- primary$rank_p
  primary$welch_p <- p.adjust(primary$welch_p, method = "holm")
  primary$rank_p  <- p.adjust(primary$rank_p,  method = "holm")
}

hdr("5b. SECONDARY -- half-life under the alternative decay models")
secondary <- bind_rows(
  compare_metric(hl, "half_life_noplateau", "RFP half-life, no plateau"),
  compare_metric(hl, "half_life_loglinear", "RFP half-life, log-linear (workbook)"))
save_table(bind_rows(primary, secondary), "08_group_comparisons")

## 6. HEADLINE SUMMARY
hdr("6. HEADLINE SUMMARY")
star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***"); if (p < 0.01) return("**")
  if (p < 0.05)  return("*");   if (p < 0.1)  return(".")
  "ns"
}
omni <- function(m, term) { av <- anova(m); if (term %in% rownames(av)) av[term, "Pr(>F)"] else NA }

headline <- data.frame(
  test = c("RFP: Group x Day  (do groups differ in RATE?)",
           "RFP: Group        (do groups differ in LEVEL?)",
           "GFP: Group x Day  (do groups differ in RATE?)",
           "GFP: Group        (do groups differ in LEVEL?)"),
  p = c(omni(m_rfp, "Day:Group"), omni(m_rfp, "Group"),
        omni(m_gfp, "Day:Group"), omni(m_gfp, "Group")))
headline <- rbind(headline, data.frame(
  test = paste0(primary$metric, ":  ", primary$contrast, "  (rank test)"),
  p = primary$rank_p))
headline$signif <- vapply(headline$p, star, character(1))
show_df(headline, digits = 3)
save_table(headline, "09_headline_summary")

cat("
Reading this table: the mixed models use every timepoint and respect the
repeated-measures structure, but they are secondary trajectory analyses. The
prespecified primary turnover-rate test is the constrained nonlinear M1-vs-M2
comparison in 04_nonlinear_mixed_model.R. Per-slice metrics are supporting
summaries and should not overrule the primary model.
")
if (!TWO_GROUP)
  cat("With ", N_GROUP, " groups, a significant omnibus Group or Day:Group term says
only that SOME groups differ. Read section 2c for which.\n", sep = "")

hdr("DONE")
cat("Tables : ", normalizePath(TABLES), "\n", sep = "")
cat("Log    : ", normalizePath(.log), "\n", sep = "")
cat("Next   : Rscript 02_figures.R   then   Rscript 04_nonlinear_mixed_model.R\n")

sink(); close(con)
saveRDS(list(dendra = dendra, hl = hl, auc = auc, desc = desc,
             m_rfp = m_rfp, m_gfp = m_gfp, primary = primary,
             groups = GROUPS, two_group = TWO_GROUP),
        file.path(RESULTS, "fitted_objects.rds"))
message("\nDone. Log written to ", .log)
