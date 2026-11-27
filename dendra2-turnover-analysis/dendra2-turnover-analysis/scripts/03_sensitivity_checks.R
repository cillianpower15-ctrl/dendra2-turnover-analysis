## 03_sensitivity_checks.R -- reviewer-facing robustness and accessibility QC.
## Works with any group names; nothing is hard-coded to Control or HSPB1.

source("00_setup.R")
fitted_obj <- readRDS(file.path(RESULTS, "fitted_objects.rds"))
rfp_random <- attr(fitted_obj$m_rfp, "structure")
gfp_random <- attr(fitted_obj$m_gfp, "structure")
if (is.null(rfp_random)) rfp_random <- "(1 | SliceID)"
if (is.null(gfp_random)) gfp_random <- "(1 | SliceID)"

.log <- file.path(RESULTS, "sensitivity_log.txt")
con <- file(.log, open = "wt"); sink(con, split = TRUE)
hdr <- function(x) cat("\n\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

hdr("SENSITIVITY AND METHOD CHECKS")
describe_run()

hdr("1. MIXED MODEL vs ORDINARY LINEAR MODEL")
cat("\nThe mixed model is primary because the same slice is measured repeatedly.\n",
    "The ordinary linear model is shown only to demonstrate the effect of\n",
    "incorrectly treating repeated observations as independent.\n", sep = "")

cmp <- function(response, data, label) {
  f_lmer <- lmm_formula(response, "(1 | SliceID)")
  f_lm   <- reformulate(LMM_FIXED_TERMS, response = response)
  m1 <- suppressMessages(lmer(f_lmer, data = data))
  m2 <- lm(f_lm, data = data)
  c1 <- summary(m1)$coefficients; c2 <- summary(m2)$coefficients
  data.frame(channel = label, term = rownames(c1),
             estimate_mixed = c1[, 1], p_mixed = c1[, ncol(c1)],
             estimate_ordinary = c2[, 1], p_ordinary = c2[, 4],
             row.names = NULL)
}
model_comparison <- bind_rows(
  cmp("RFP_dFF0", rfp_model, "RFP"),
  cmp("GFP_dFF0", gfp_model, "GFP"))
print(as.data.frame(model_comparison), digits = 4, row.names = FALSE)
save_table(model_comparison, "10_mixed_vs_ordinary_model")

hdr("2. RANDOM INTERCEPT vs RANDOM INTERCEPT + SLOPE")
random_structure <- list()
for (nm in c("RFP", "GFP")) {
  d <- if (nm == "RFP") rfp_model else gfp_model
  response <- if (nm == "RFP") "RFP_dFF0" else "GFP_dFF0"
  f_i <- lmm_formula(response, "(1 | SliceID)")
  f_s <- lmm_formula(response, "(Day | SliceID)")
  m_i <- suppressMessages(lmer(f_i, data = d, REML = FALSE))
  m_s <- suppressMessages(suppressWarnings(try(lmer(f_s, data = d, REML = FALSE),
                                                silent = TRUE)))
  cat("\n", nm, ":\n", sep = "")
  if (inherits(m_s, "try-error")) {
    cat("  Random-slope model failed; random intercept is retained.\n")
    random_structure[[nm]] <- data.frame(channel = nm, structure = "random slope",
      AIC = NA, singular = NA, LRT_p = NA, conclusion = "failed to fit")
  } else {
    lr <- suppressMessages(anova(m_i, m_s))
    cat(sprintf("  intercept only : AIC %.1f  singular %s\n", AIC(m_i), isSingular(m_i)))
    cat(sprintf("  + random slope : AIC %.1f  singular %s\n", AIC(m_s), isSingular(m_s)))
    print(lr)
    supported <- !isSingular(m_s)
    random_structure[[nm]] <- data.frame(channel = nm, structure = "random slope",
      AIC = AIC(m_s), singular = isSingular(m_s), LRT_p = lr$`Pr(>Chisq)`[2],
      conclusion = if (supported) "estimable: prespecified slope retained" else "singular: intercept fallback")
  }
}
save_table(bind_rows(random_structure), "11_random_structure_check")

hdr("3. LEAVE-ONE-SLICE-OUT ROBUSTNESS")
cat("\nEach row drops one biological slice and refits the model. Omnibus Group\n",
    "and Day:Group terms are used, so this works with any number or names of groups.\n", sep = "")

get_p <- function(model, term) {
  av <- anova(model)
  if (term %in% rownames(av)) av[term, "Pr(>F)"] else NA_real_
}

loso <- function(data, response, label, term, random_term) {
  fml <- lmm_formula(response, random_term)
  full <- suppressMessages(lmer(fml, data = data))
  p_full <- get_p(full, term)
  out <- map_dfr(unique(data$SliceID), function(s) {
    m <- suppressMessages(suppressWarnings(try(
      lmer(fml, data = filter(data, SliceID != s)), silent = TRUE)))
    data.frame(channel = label, term = term, dropped = s,
               p = if (inherits(m, "try-error")) NA_real_ else get_p(m, term))
  })
  valid <- out$p[is.finite(out$p)]
  cat("\n", label, " -- ", term, "  using ", random_term, "\n", sep = "")
  cat(sprintf("  full model p = %.4g\n", p_full))
  if (length(valid))
    cat(sprintf("  leave-one-out range: %.4g to %.4g\n", min(valid), max(valid)))
  print(out[order(-out$p), ], digits = 3, row.names = FALSE)
  mutate(out, full_model_p = p_full,
         robustness = ifelse(is.finite(p) & ((p_full < 0.05) == (p < 0.05)),
                             "same significance classification", "classification changed"))
}

loso_tbl <- bind_rows(
  loso(rfp_model, "RFP_dFF0", "RFP", "Group", rfp_random),
  loso(rfp_model, "RFP_dFF0", "RFP", "Day:Group", rfp_random),
  loso(gfp_model, "GFP_dFF0", "GFP", "Group", gfp_random),
  loso(gfp_model, "GFP_dFF0", "GFP", "Day:Group", gfp_random))
save_table(loso_tbl, "12_leave_one_slice_out")

hdr("4. FIGURE PALETTE -- COLOUR-VISION-DEFICIENCY CHECK")
if (!requireNamespace("colorspace", quietly = TRUE) ||
    !requireNamespace("farver", quietly = TRUE)) {
  cat("  Optional packages 'colorspace' and/or 'farver' are absent; palette\n",
      "  simulation skipped. Figures still use colour plus shape and line type.\n", sep = "")
} else {
  oklab <- function(hex)
    farver::convert_colour(t(col2rgb(hex)), from = "rgb", to = "oklab")
  dE <- function(a, b) sqrt(sum((oklab(a) - oklab(b))^2)) * 100
  rel_lum <- function(hex) {
    v <- col2rgb(hex)[, 1] / 255
    v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
    sum(v * c(0.2126, 0.7152, 0.0722))
  }
  contrast <- function(a, b) {
    l <- sort(c(rel_lum(a), rel_lum(b)), decreasing = TRUE)
    (l[1] + 0.05) / (l[2] + 0.05)
  }
  pair_idx <- combn(seq_len(N_GROUP), 2, simplify = FALSE)
  palette_tbl <- bind_rows(lapply(pair_idx, function(ix) {
    a <- PAL[ix[1]]; b <- PAL[ix[2]]
    data.frame(group1 = names(a), group2 = names(b), colour1 = unname(a), colour2 = unname(b),
      dE_normal = dE(a, b),
      dE_deuteranopia = dE(colorspace::deutan(a), colorspace::deutan(b)),
      dE_protanopia = dE(colorspace::protan(a), colorspace::protan(b)),
      dE_tritanopia = dE(colorspace::tritan(a), colorspace::tritan(b)),
      contrast1_white = contrast(a, "#FFFFFF"),
      contrast2_white = contrast(b, "#FFFFFF"))
  }))
  print(palette_tbl, digits = 3, row.names = FALSE)
  save_table(palette_tbl, "13_palette_accessibility")
}

hdr("DONE")
cat("Log: ", normalizePath(.log), "\n", sep = "")
sink(); close(con)
message("\nSensitivity checks complete. Log written to ", .log)
