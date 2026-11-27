## 07_factorial_figure.R -- the clearest view of a crossed design.

source("00_setup.R")
suppressPackageStartupMessages(library(nlme))

parse_groups <- function(g) {
  bg <- trimws(sub("(?i)\\s*GFAP.*$", "", g, perl = TRUE))
  cons <- ifelse(grepl("control", g, ignore.case = TRUE), "Control",
          ifelse(grepl("hspb1",   g, ignore.case = TRUE), "HSPB1", NA))
  data.frame(Group = g, Background = bg, Construct = cons, stringsAsFactors = FALSE)
}
gmap <- parse_groups(GROUPS)

if (anyNA(gmap$Construct) || length(unique(gmap$Background)) < 2 ||
    length(unique(gmap$Construct)) != 2) {
  message("07_factorial_figure.R: group names are not a Background x Control/HSPB1 ",
          "grid; skipping this figure.")
  quit(save = "no", status = 0)
}
message("Factorial figure: ", length(unique(gmap$Background)), " backgrounds x 2 constructs")

## same colour for a construct in every panel, so panels can be compared
CONS_PAL <- c(Control = "#0072B2", HSPB1 = "#D55E00")
bg_levels <- unique(gmap$Background)

## left_join against a character-keyed lookup drops Group's factor levels, and
## nlme then cannot build contrasts for it. Restore the factor explicitly.
dat <- decay_fit %>% as.data.frame() %>%
  mutate(Group = as.character(Group)) %>%
  left_join(gmap, by = "Group") %>%
  mutate(Group = factor(Group, levels = GROUPS),
         Background = factor(Background, levels = bg_levels),
         Construct = factor(Construct, levels = c("Control", "HSPB1")))

## measured points per day: geometric mean +/- SEM on the log scale, to match
## the log axis and the multiplicative model
obs <- dat %>%
  filter(is.finite(frac), frac > 0) %>%
  group_by(Background, Construct, Day) %>%
  summarise(n = dplyr::n(), m = exp(mean(log(frac))),
            se = sd(log(frac)) / sqrt(dplyr::n()), .groups = "drop") %>%
  mutate(lo = m * exp(-se), hi = m * exp(se))

## ---- refit the primary model so the fitted curves carry uncertainty --------
grp <- groupedData(frac ~ t | SliceID, data = dat)
m2 <- nlme(frac ~ (1 - plogis(logitP)) * exp(-exp(logK) * t) + plogis(logitP),
           data = grp,
           fixed = list(logK = reformulate("Group", response = "logK"),
                        logitP = reformulate("1", response = "logitP")),
           random = logK ~ 1 | SliceID,
           start = c(log(0.15), rep(0, N_GROUP - 1), qlogis(0.25)),
           method = "ML", na.action = na.fail,
           control = nlmeControl(maxIter = 300, msMaxIter = 300,
                                 pnlsMaxIter = 30, pnlsTol = 0.05,
                                 returnObject = TRUE))
beta <- fixef(m2); V <- vcov(m2)
k_names <- grep("^logK", names(beta), value = TRUE)
p_name <- if ("logitP" %in% names(beta)) "logitP" else "logitP.(Intercept)"
Xg <- model.matrix(~ Group, data.frame(Group = factor(GROUPS, levels = GROUPS)))
rownames(Xg) <- GROUPS
tcrit <- qt(0.975, summary(m2)$fixDF$X[[1]])

band <- bind_rows(lapply(GROUPS, function(g) {
  xk <- Xg[g, ]
  K <- exp(sum(xk * beta[k_names])); P <- plogis(unname(beta[p_name]))
  tt <- seq(0, max(dat$t), length.out = 200)
  fit <- (1 - P) * exp(-K * tt) + P
  dK <- -(1 - P) * tt * K * exp(-K * tt); dP <- P * (1 - P) * (1 - exp(-K * tt))
  se <- vapply(seq_along(tt), function(i) {
    gv <- setNames(rep(0, length(beta)), names(beta))
    gv[k_names] <- dK[i] * xk; gv[p_name] <- dP[i]
    sqrt(as.numeric(t(gv) %*% V %*% gv)) }, numeric(1))
  data.frame(Group = g, Day = tt + BASELINE_RFP, fit = fit,
             ## floored above zero so the band survives a log y-axis
             lo = pmax(0.02, fit - 1.96 * se), hi = pmin(1, fit + 1.96 * se))
})) %>% mutate(Group = as.character(Group)) %>%
  left_join(gmap, by = "Group") %>%
  mutate(Background = factor(Background, levels = bg_levels),
         Construct = factor(Construct, levels = c("Control", "HSPB1")))

## area between the two fitted curves = the HSPB1 effect
gap <- band %>%
  select(Background, Day, Construct, fit) %>%
  tidyr::pivot_wider(names_from = Construct, values_from = fit) %>%
  filter(!is.na(Control), !is.na(HSPB1))

## model half-life per group; colour keys each line to its own curve
hl_lab <- bind_rows(lapply(bg_levels, function(b) {
  gs <- gmap$Group[gmap$Background == b]
  cons <- gmap$Construct[match(gs, gmap$Group)]
  hl <- log(2) / exp(vapply(gs, function(g) sum(Xg[g, ] * beta[k_names]), numeric(1)))
  data.frame(Background = b, Construct = cons,
             ## plotmath, not the ½ glyph: the glyph is missing from some
             ## default fonts and silently renders as "t..". Drawn with parse = TRUE.
             lab = sprintf("t[1/2]~%.1f~d", hl),
             row = match(cons, c("Control", "HSPB1")))
}))
hl_lab$Background <- factor(hl_lab$Background, levels = bg_levels)

## Per-panel y-ranges, set via geom_blank + free_y. Each background is zoomed to
## its own data; note this means vertical gaps are NOT comparable between panels.
rng <- dat %>% group_by(Background) %>%
  summarise(dmin = min(frac, na.rm = TRUE), .groups = "drop")
## headroom above 1.0 holds the per-day p-values in two fixed rows
LAB_ROW_LO <- 1.18
LAB_ROW_HI <- 1.52
panel_lims <- bind_rows(
  rng %>% transmute(Background, frac = dmin * 0.93),
  rng %>% transmute(Background, frac = 1.66))
## rows 1-2 are the half-lives; row 3 is left for the overall verdict below
hl_lab <- hl_lab %>%
  left_join(rng, by = "Background") %>%
  mutate(y = dmin * 0.99 * (1.20 ^ (3 - row)),
         Construct = factor(Construct, levels = c("Control", "HSPB1")))

## ---- within-background HSPB1 vs Control contrasts --------------------------
## Not the same as the against-a-single-reference contrasts in 06.
contrasts_bg <- bind_rows(lapply(bg_levels, function(b) {
  gc <- gmap$Group[gmap$Background == b & gmap$Construct == "Control"]
  gh <- gmap$Group[gmap$Background == b & gmap$Construct == "HSPB1"]
  if (!length(gc) || !length(gh)) return(NULL)
  gv <- setNames(rep(0, length(beta)), names(beta))
  gv[k_names] <- Xg[gh[1], ] - Xg[gc[1], ]        # log K(HSPB1) - log K(Control)
  est <- sum(gv * beta); se <- sqrt(as.numeric(t(gv) %*% V %*% gv))
  data.frame(Background = b, ratio = exp(est),
             lo = exp(est - tcrit * se), hi = exp(est + tcrit * se),
             p_raw = 2 * pt(-abs(est / se), summary(m2)$fixDF$X[[1]]))
}))
contrasts_bg$p_adj <- p.adjust(contrasts_bg$p_raw, "holm")
contrasts_bg$Background <- factor(contrasts_bg$Background, levels = bg_levels)
stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else
  if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
contrasts_bg$label <- sprintf("%.2f (%.2f-%.2f) %s", contrasts_bg$ratio,
                              contrasts_bg$lo, contrasts_bg$hi,
                              vapply(contrasts_bg$p_adj, stars, character(1)))
save_table(contrasts_bg, "17_within_background_HSPB1_vs_Control")

cat("\nHSPB1 vs Control within each background (>1 = HSPB1 clears faster):\n")
print(contrasts_bg[, c("Background", "ratio", "lo", "hi", "p_raw", "p_adj")],
      digits = 3, row.names = FALSE)

## ---- does the HSPB1 effect DIFFER between backgrounds? ----------------------
## The two contrasts above each test one background against its own control.
## Neither tests whether the effects differ from each other, which is what a
## claim of background-dependence asserts. Significant in one background and
## not in the other is not evidence of a difference (Gelman & Stern 2006).
## Tested two ways, which are algebraically equivalent here because logK ~ Group
## is saturated for a 2 x 2 design:
##   (a) a Wald contrast on the fitted model -- the difference of differences;
##   (b) an ML likelihood-ratio test of the additive model against this one.
interaction_rows <- list()
if (length(bg_levels) >= 2) {
  pairs <- utils::combn(bg_levels, 2, simplify = FALSE)
  for (pr in pairs) {
    b1 <- pr[1]; b2 <- pr[2]
    g <- function(bg, cons) gmap$Group[gmap$Background == bg & gmap$Construct == cons][1]
    need <- c(g(b1, "Control"), g(b1, "HSPB1"), g(b2, "Control"), g(b2, "HSPB1"))
    if (anyNA(need)) next
    gv <- setNames(rep(0, length(beta)), names(beta))
    ## (effect in b2) - (effect in b1), on the log scale
    gv[k_names] <- (Xg[need[4], ] - Xg[need[3], ]) - (Xg[need[2], ] - Xg[need[1], ])
    est <- sum(gv * beta); se <- sqrt(as.numeric(t(gv) %*% V %*% gv))
    interaction_rows[[length(interaction_rows) + 1]] <- data.frame(
      comparison = paste0("HSPB1 effect in ", b2, " vs in ", b1),
      ratio_of_ratios = exp(est),
      lo = exp(est - tcrit * se), hi = exp(est + tcrit * se),
      p_wald = 2 * pt(-abs(est / se), summary(m2)$fixDF$X[[1]]))
  }
}

## additive model: background and construct both act on logK, but the construct
## effect is forced to be the same size in every background
lrt_p <- NA_real_; lrt_chi <- NA_real_; lrt_df <- NA_integer_
if (length(bg_levels) >= 2) {
  n_add <- ncol(model.matrix(~ Background + Construct, dat))
  m_add <- try(nlme(frac ~ (1 - plogis(logitP)) * exp(-exp(logK) * t) + plogis(logitP),
                    data = grp,
                    fixed = list(logK = logK ~ Background + Construct,
                                 logitP = logitP ~ 1),
                    random = logK ~ 1 | SliceID,
                    start = c(log(0.15), rep(0, n_add - 1), qlogis(0.25)),
                    method = "ML", na.action = na.fail,
                    control = nlmeControl(maxIter = 300, msMaxIter = 300,
                                          pnlsMaxIter = 30, pnlsTol = 0.05,
                                          returnObject = TRUE)),
               silent = TRUE)
  if (!inherits(m_add, "try-error")) {
    cmp <- anova(m_add, m2)
    lrt_chi <- cmp[["L.Ratio"]][2]
    lrt_df <- cmp[["df"]][2] - cmp[["df"]][1]
    lrt_p <- cmp[["p-value"]][2]
  }
}

if (length(interaction_rows)) {
  interaction_tbl <- bind_rows(interaction_rows)
  interaction_tbl$p_LRT_overall <- lrt_p
  interaction_tbl$LRT_chisq <- lrt_chi
  interaction_tbl$LRT_df <- lrt_df
  interaction_tbl$note <-
    "Ratio of ratios: >1 means HSPB1 accelerates turnover more in the first background named. Tests whether the HSPB1 effect DIFFERS between backgrounds; this is the test a background-dependence claim requires."
  save_table(interaction_tbl, "19_background_x_construct_interaction")

  cat("\nDoes the HSPB1 effect differ between backgrounds?\n")
  print(interaction_tbl[, c("comparison", "ratio_of_ratios", "lo", "hi", "p_wald")],
        digits = 3, row.names = FALSE)
  if (is.finite(lrt_p))
    cat(sprintf("Likelihood-ratio test of interaction: chisq = %.2f, df = %d, p = %.3f\n",
                lrt_chi, lrt_df, lrt_p))
  if (is.finite(lrt_p) && lrt_p >= 0.05)
    cat("NOTE: the interaction is not significant. The within-background results\n",
        "  above do not by themselves establish that the HSPB1 effect is\n",
        "  background-dependent; report that separation as suggested, not shown.\n", sep = "")
}

## ---- per-day cross-sectional tests, HSPB1 vs Control within background -----
perday <- dat %>%
  filter(Day >= 3) %>%
  group_by(Background, Day) %>%
  summarise(
    n_control = sum(Construct == "Control"),
    n_hspb1 = sum(Construct == "HSPB1"),
    mean_control = mean(frac[Construct == "Control"]),
    mean_hspb1 = mean(frac[Construct == "HSPB1"]),
    p_welch = tryCatch(t.test(frac ~ Construct)$p.value,
                       error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(note = "Unadjusted Welch t-test at a single day; descriptive only")
save_table(perday, "18_perday_RFP_HSPB1_vs_Control")

## Labels name their day: days 3/5/7 sit close together and the labels have to
## be staggered vertically, which breaks the visual link back to the x-tick.
ann <- perday %>%
  filter(!is.na(p_welch)) %>%
  mutate(lab = paste0("Day ", format(Day, trim = TRUE), "\n",
                      ifelse(p_welch < 0.001, "p<0.001", sprintf("p=%.2f", p_welch))),
         sig = p_welch < 0.05)

## Fixed rows rather than band-relative heights: data never exceeds 1.0, so
## rows above that cannot collide with it.
ann <- ann %>% group_by(Background) %>%
  mutate(y = ifelse(row_number() %% 2 == 0, LAB_ROW_HI, LAB_ROW_LO)) %>%
  ungroup()

## overall within-background verdict, written out rather than starred
ov_lab <- contrasts_bg %>%
  transmute(Background,
            lab = sprintf("overall p = %s %s",
                          ifelse(p_adj < 0.001, "<0.001", sprintf("%.3f", p_adj)),
                          vapply(p_adj, stars, character(1))),
            sig = p_adj < 0.05) %>%
  left_join(rng, by = "Background") %>%
  mutate(y = dmin * 0.99)
ov_lab$Background <- factor(ov_lab$Background, levels = bg_levels)

## ---- PANEL A: decay, one small panel per background ------------------------
brk <- sort(unique(c(BASELINE_RFP, dat$Day)))
pA <- ggplot() +
  ## the effect itself, shaded: the area between the two fitted curves
  geom_ribbon(data = gap, aes(Day, ymin = HSPB1, ymax = Control),
              fill = "grey35", alpha = 0.22, colour = NA) +
  geom_ribbon(data = band, aes(Day, ymin = lo, ymax = hi, fill = Construct),
              alpha = 0.13, colour = NA) +
  ## individual slices: visible enough to trace, subordinate to the fit
  geom_line(data = dat, aes(Day, frac, group = SliceID, colour = Construct),
            linewidth = 0.22, alpha = 0.30) +
  geom_line(data = band, aes(Day, fit, colour = Construct), linewidth = 1.05) +
  ## the measurements themselves: geometric mean +/- SEM at each imaged day
  geom_errorbar(data = obs, aes(Day, ymin = lo, ymax = hi, colour = Construct),
                width = 0.55, linewidth = 0.35, show.legend = FALSE,
                position = position_dodge(width = 0.55), inherit.aes = FALSE) +
  geom_point(data = obs, aes(Day, m, fill = Construct), shape = 21, size = 1.75,
             colour = "white", stroke = 0.4, show.legend = FALSE,
             position = position_dodge(width = 0.55), inherit.aes = FALSE) +
  ## Uniform weight, deliberately: these are unadjusted single-day Welch tests,
  ## so bolding p < 0.05 would imply a claim the adjusted analysis does not make.
  geom_text(data = ann, aes(Day, y, label = lab), colour = "grey40",
            size = 1.85, lineheight = 0.9, show.legend = FALSE,
            inherit.aes = FALSE) +
  ## half-life per construct, colour-keyed to its curve
  geom_text(data = hl_lab, aes(x = BASELINE_RFP, y = y, label = lab,
                               colour = Construct),
            hjust = 0, vjust = 0, size = 2.3, fontface = "bold", parse = TRUE,
            show.legend = FALSE, inherit.aes = FALSE) +
  ## overall verdict, spelled out under the half-lives
  geom_text(data = ov_lab, aes(x = BASELINE_RFP, y = y, label = lab,
                               fontface = I(ifelse(sig, "bold", "plain"))),
            hjust = 0, vjust = 0, size = 2.1, colour = "grey25",
            show.legend = FALSE, inherit.aes = FALSE) +
  geom_blank(data = panel_lims, aes(x = BASELINE_RFP, y = frac)) +
  facet_wrap(~Background, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = CONS_PAL, name = NULL) +
  scale_fill_manual(values = CONS_PAL, name = NULL) +
  scale_x_continuous(breaks = brk, labels = day_labels,
                     expand = expansion(mult = c(0.02, 0.07))) +
  ## Log y-axis: constant rate ratio = constant vertical gap. Limits are set
  ## per panel by geom_blank + free_y, so each background fills its own range.
  scale_y_log10(breaks = c(0.2, 0.3, 0.5, 0.7, 1.0),
                labels = c("0.2", "0.3", "0.5", "0.7", "1.0")) +
  labs(x = "Day", y = expression(Photoconverted~signal~retained~(F/F[0]))) +
  theme_pub() + theme(legend.position = "top")

## ---- PANEL B: the within-background effect ---------------------------------
xr <- range(c(contrasts_bg$lo, contrasts_bg$hi, 1))
pB <- ggplot(contrasts_bg, aes(ratio, Background)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.35, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.10,
                linewidth = 0.45, colour = "grey25") +
  geom_point(size = 2.8, shape = 21, fill = CONS_PAL[["HSPB1"]],
             colour = "white", stroke = 0.4) +
  geom_text(aes(label = label), vjust = -1.3, size = 2.5, colour = "grey20") +
  scale_x_log10(limits = c(xr[1] * 0.6, xr[2] * 1.7)) +
  scale_y_discrete(expand = expansion(add = 0.42)) +
  labs(x = "Turnover rate ratio, HSPB1 / Control (log scale)", y = NULL) +
  theme_pub()

## Tagged C and D: panels A and B of the assembled figure are the schematic and
## the representative images, neither produced here.
n_slice <- dendra %>% distinct(Group, Slice) %>% count(Group)
n_cap <- if (length(unique(n_slice$n)) == 1)
  sprintf("n = %d slices per group", n_slice$n[1]) else
  paste0("n = ", paste(range(n_slice$n), collapse = "-"), " slices per group")
fig <- (pA / pB) + plot_layout(heights = c(1, 0.55)) +
  plot_annotation(tag_levels = list(c("C", "D")),
                  caption = paste0(n_cap,
                    ". Points, geometric mean +/- SEM; per-day P values unadjusted."),
                  theme = theme(plot.caption = element_text(
                    size = 6, colour = "grey35", hjust = 0)))
save_fig(fig, "Figure5_by_background", W2, 122)
save_fig(pA + theme(plot.margin = margin(3, 3, 3, 3)),
         "Figure5A_decay_by_background", W2, 72)
save_fig(pB + theme(plot.margin = margin(3, 3, 3, 3)),
         "Figure5B_HSPB1_effect_by_background", max(W1, W2 * 0.6),
         30 + 13 * nrow(contrasts_bg))

## ---- figure legend, as manuscript text --------------------------------
n_by <- dendra %>% distinct(Group, Slice) %>% count(Group)
n_txt <- paste(sprintf("%s, n = %d slices", n_by$Group, n_by$n), collapse = "; ")
contrast_txt <- paste(sprintf(
  "%s, rate ratio (HSPB1 / Control) %.2f (95%% CI %.2f-%.2f, P = %s)",
  contrasts_bg$Background, contrasts_bg$ratio, contrasts_bg$lo, contrasts_bg$hi,
  format.pval(contrasts_bg$p_adj, digits = 2, eps = 1e-4)), collapse = "; ")

## Title states only what the analysis supports. A background-dependence claim
## needs the interaction, which is not significant (table 19), so the earlier
## title contradicted a table shipped in the same folder.
legend_txt <- c(
  "FIGURE 1. Astrocytic HSPB1 accelerates turnover of the photoconverted tau pool.",
  "",
  paste0("(C) Fraction of photoconverted signal retained over the chase, split by ",
         "genetic background (panels) and construct (colour). Points, geometric mean ",
         "+/- SEM of the measured slices at each imaged day; thin lines, individual ",
         "slices; thick lines, population fit of the constrained nonlinear ",
         "mixed-effects model; shaded bands, 95% confidence intervals of the fit; grey ",
         "shading between the fitted curves, the HSPB1 effect. The day-",
         BASELINE_RFP_DISPLAY,
         " post-conversion point defines fraction = 1 and was excluded from the ",
         "likelihood. Note the ",
         "independent y-axis range per panel. Day-labelled p-values are unadjusted ",
         "Welch t-tests at single timepoints and are descriptive only; the ",
         "\"overall p\" in each panel is the model-based within-background contrast ",
         "reported in (D)."),
  "",
  paste0("(D) Turnover rate ratio (HSPB1 / Control) estimated WITHIN each background, ",
         "with 95% CI (log scale). Values > 1 indicate HSPB1 clears the photoconverted ",
         "pool faster than the matched control; < 1 indicates slower clearance. ",
         contrast_txt, "."),
  "",
  paste0(n_txt, ". Rate constants were compared by a prespecified maximum-likelihood ",
         "likelihood-ratio test on the group effect on log(K); the within-background ",
         "contrasts in (D) are Holm-adjusted across backgrounds. Asterisks: *P < 0.05, ",
         "**P < 0.01, ***P < 0.001; n.s., not significant."),
  "",
  ## the interaction and the plateau-model sensitivity both temper the
  ## background-dependence reading, so they belong in the legend, not only in
  ## the tables where a reader has to go looking for them
  if (exists("interaction_tbl")) paste0(
    "The HSPB1 effect did not differ significantly between backgrounds ",
    "(background x construct interaction on log(K): ratio of ratios ",
    sprintf("%.2f, 95%% CI %.2f-%.2f, P = %.2f", interaction_tbl$ratio_of_ratios[1],
            interaction_tbl$lo[1], interaction_tbl$hi[1], interaction_tbl$p_wald[1]),
    "), so the difference between backgrounds is suggested but not established."),
  "",
  paste0("Sensitivity: the primary model shares one plateau across groups. ",
         "In the secondary model allowing group-specific plateaus, which fits better, ",
         "the WT-tau HSPB1 rate effect also reaches significance, indicating that ",
         "rate and plateau are partially confounded in this dataset."),
  "",
  "NOTE: the slice is the repeated-measures unit. Report the number of animals or",
  "preparations separately; slices sharing a preparation are not independent.")
writeLines(legend_txt, file.path(RESULTS, "Figure1_MAIN_LEGEND.txt"))
message("  legend: Figure1_MAIN_LEGEND.txt")

## ---- short legend -------------------------------------------------------
## Same panels in four lines, for caption limits and for talks.
n_per <- min(n_by$n)
short_txt <- c(
  "FIGURE 1. Astrocytic HSPB1 speeds up clearance of photoconverted tau.",
  "",
  paste0("(C) Photoconverted (red) signal remaining after the pulse, as a fraction of ",
         "the signal at conversion. Blue, control; orange, HSPB1. Points are the ",
         "measured means with SEM, thick lines are the model fit, shaded bands are ",
         "95% confidence intervals, and thin lines are individual slices. Grey ",
         "shading marks the gap between the two fits. t½ is the fitted half-life. ",
         "P values above each day compare the two constructs at that day only. ",
         "Note the y-axis range differs between panels."),
  "",
  paste0("(D) How much faster HSPB1 clears the pool than its matched control, with ",
         "95% confidence intervals. A value of 1 means no difference."),
  "",
  paste0("The effect was significant in Mut-tau and not in WT-tau, but the two were ",
         "not significantly different from each other, so the difference between ",
         "backgrounds is suggested rather than shown."),
  "",
  paste0("n = ", n_per, " slices per group. *P < 0.05; n.s., not significant."))
writeLines(short_txt, file.path(RESULTS, "Figure1_MAIN_LEGEND_SHORT.txt"))
message("  legend: Figure1_MAIN_LEGEND_SHORT.txt (plain-language version)")

message("\nFactorial figures written to ", normalizePath(FIGURES))
