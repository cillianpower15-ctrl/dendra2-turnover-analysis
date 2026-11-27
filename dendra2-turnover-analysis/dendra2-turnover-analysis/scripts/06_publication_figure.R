## 06_publication_figure.R -- the composite main-text figure.

source("00_setup.R")
suppressPackageStartupMessages(library(nlme))

message("Building the composite publication figure...")

## ---- refit the primary model (nlme objects are not carried in the RDS) -----
dn <- decay_fit %>% as.data.frame()
grp <- groupedData(frac ~ t | SliceID, data = dn)
ctrl <- nlmeControl(maxIter = 300, msMaxIter = 300, pnlsMaxIter = 30,
                    pnlsTol = 0.05, returnObject = TRUE)
fixed_m2 <- list(logK = reformulate("Group", response = "logK"),
                 logitP = reformulate("1", response = "logitP"))
m2 <- nlme(frac ~ (1 - plogis(logitP)) * exp(-exp(logK) * t) + plogis(logitP),
           data = grp, fixed = fixed_m2, random = logK ~ 1 | SliceID,
           start = c(log(0.15), rep(0, N_GROUP - 1), qlogis(0.25)),
           method = "ML", na.action = na.fail, control = ctrl)

beta <- fixef(m2); V <- vcov(m2)
k_names <- grep("^logK", names(beta), value = TRUE)
p_name <- if ("logitP" %in% names(beta)) "logitP" else "logitP.(Intercept)"
df_fix <- summary(m2)$fixDF$X[[1]]
tcrit <- qt(0.975, df_fix)

Xg <- model.matrix(~ Group, data = data.frame(Group = factor(GROUPS, levels = GROUPS)))
rownames(Xg) <- GROUPS

## ---- fitted curve + 95% band by the delta method ---------------------------
curve_for <- function(g, tt) {
  xk <- Xg[g, ]
  K <- exp(sum(xk * beta[k_names])); P <- plogis(unname(beta[p_name]))
  fit <- (1 - P) * exp(-K * tt) + P
  d_logK <- -(1 - P) * tt * K * exp(-K * tt)
  d_logitP <- P * (1 - P) * (1 - exp(-K * tt))
  se <- vapply(seq_along(tt), function(i) {
    gv <- setNames(rep(0, length(beta)), names(beta))
    gv[k_names] <- d_logK[i] * xk; gv[p_name] <- d_logitP[i]
    sqrt(as.numeric(t(gv) %*% V %*% gv))
  }, numeric(1))
  data.frame(Group = g, Day = tt + BASELINE_RFP, fit = fit,
             lo = pmax(0, fit - 1.96 * se), hi = pmin(1, fit + 1.96 * se))
}
band <- bind_rows(lapply(GROUPS, curve_for, tt = seq(0, max(dn$t), length.out = 200)))
band$Group <- factor(band$Group, levels = GROUPS)

## ---- per-group rate, half-life, and pairwise ratios ------------------------
rates <- bind_rows(lapply(GROUPS, function(g) {
  gv <- setNames(rep(0, length(beta)), names(beta)); gv[k_names] <- Xg[g, ]
  est <- sum(Xg[g, ] * beta[k_names]); se <- sqrt(as.numeric(t(gv) %*% V %*% gv))
  data.frame(Group = g, half_life = log(2) / exp(est),
             hl_lo = log(2) / exp(est + tcrit * se),
             hl_hi = log(2) / exp(est - tcrit * se))
}))
rates$Group <- factor(rates$Group, levels = GROUPS)

pair_list <- combn(GROUPS, 2, simplify = FALSE)
pairs_df <- bind_rows(lapply(pair_list, function(p) {
  gv <- setNames(rep(0, length(beta)), names(beta))
  gv[k_names] <- Xg[p[2], ] - Xg[p[1], ]
  est <- sum(gv * beta); se <- sqrt(as.numeric(t(gv) %*% V %*% gv))
  data.frame(group1 = p[1], group2 = p[2], ratio = exp(est),
             lo = exp(est - tcrit * se), hi = exp(est + tcrit * se),
             p_raw = 2 * pt(-abs(est / se), df_fix))
}))
pairs_df$p_adj <- p.adjust(pairs_df$p_raw, "holm")

## ---- per-slice half-life ---------------------------------------------------
per_slice <- coef(m2)
slice_map <- dn %>% distinct(SliceID, Group)
slice_map <- slice_map[match(rownames(per_slice), slice_map$SliceID), ]
logK_cols <- grep("^logK", names(per_slice), value = TRUE)
per_slice_K <- data.frame(
  SliceID = rownames(per_slice), Group = slice_map$Group,
  logK = vapply(seq_len(nrow(per_slice)), function(i)
    sum(as.numeric(per_slice[i, logK_cols]) * Xg[as.character(slice_map$Group[i]), ]),
    numeric(1)))
per_slice_K$half_life <- log(2) / exp(per_slice_K$logK)

## ---- significance annotation ------------------------------------------------
## Conventional asterisk coding. Exact p-values also go in the legend text,
## because asterisks alone are not an acceptable substitute in most journals.
stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
}

## Day-axis labels: drop trailing ".0" so integer days read as "3" not "3.0".
## With six timepoints the decimals alone are enough to make the axis collide.

## Push direct end-of-curve labels apart vertically when curves converge.
declash <- function(y, min_gap) {
  o <- order(y); y2 <- y[o]
  for (i in seq_along(y2)[-1])
    if (y2[i] - y2[i - 1] < min_gap) y2[i] <- y2[i - 1] + min_gap
  out <- numeric(length(y)); out[o] <- y2; out
}

## Brackets are drawn only for comparisons that are significant after
sig_pairs <- if (TWO_GROUP) {
  pairs_df
} else if (N_GROUP <= 3) {
  pairs_df[pairs_df$p_adj < 0.05, , drop = FALSE]
} else {
  pairs_df[0, , drop = FALSE]
}

sig_layer <- function(sig, yvals) {
  if (!nrow(sig)) return(list())
  span <- diff(range(yvals)); top <- max(yvals)
  step <- span * 0.13
  out <- list()
  for (i in seq_len(nrow(sig))) {
    x1 <- match(sig$group1[i], GROUPS); x2 <- match(sig$group2[i], GROUPS)
    y <- top + step * i
    out <- c(out, list(
      annotate("segment", x = x1, xend = x2, y = y, yend = y,
               linewidth = 0.3, colour = "grey20"),
      annotate("segment", x = x1, xend = x1, y = y, yend = y - span * 0.025,
               linewidth = 0.3, colour = "grey20"),
      annotate("segment", x = x2, xend = x2, y = y, yend = y - span * 0.025,
               linewidth = 0.3, colour = "grey20"),
      annotate("text", x = (x1 + x2) / 2, y = y + span * 0.012,
               label = stars(sig$p_adj[i]), size = 2.9, vjust = 0, colour = "grey10")))
  }
  out
}

## ---- PANEL A ----------------------------------------------------------------
## Curves are labelled DIRECTLY at their right-hand end, not by colour alone.
anchor <- dn %>% distinct(Group) %>% mutate(Day = BASELINE_RFP, frac = 1)
many <- N_GROUP > 3
ribbon_alpha <- if (many) 0.10 else 0.16
slice_alpha <- if (many) 0.12 else 0.22

end_lab <- band %>% group_by(Group) %>% slice_max(Day, n = 1) %>% ungroup()
end_lab$y_lab <- declash(end_lab$fit, min_gap = 0.052)
x_pad <- diff(range(band$Day)) * (if (many) 0.46 else 0.26)

pA <- ggplot() +
  geom_ribbon(data = band, aes(Day, ymin = lo, ymax = hi, fill = Group),
              alpha = ribbon_alpha, colour = NA) +
  geom_line(data = dn, aes(Day, frac, group = SliceID, colour = Group),
            linewidth = 0.2, alpha = slice_alpha) +
  geom_point(data = dn, aes(Day, frac, colour = Group), size = 0.55,
             alpha = if (many) 0.25 else 0.4) +
  geom_line(data = band, aes(Day, fit, colour = Group, linetype = Group), linewidth = 0.8) +
  geom_point(data = anchor, aes(Day, frac, colour = Group), shape = 1, size = 1.6) +
  geom_hline(yintercept = plogis(unname(beta[p_name])), linetype = "22",
             linewidth = 0.28, colour = "grey45") +
  geom_segment(data = end_lab,
               aes(x = Day, xend = Day + diff(range(band$Day)) * 0.012,
                   y = fit, yend = y_lab, colour = Group),
               linewidth = 0.25, alpha = 0.8) +
  geom_text(data = end_lab, aes(Day, y_lab, label = Group, colour = Group),
            hjust = 0, nudge_x = diff(range(band$Day)) * 0.018,
            size = 2.2, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
  scale_linetype_manual(values = LINES) +
  scale_x_continuous(breaks = sort(unique(c(BASELINE_RFP, dn$Day))),
                     labels = day_labels,
                     limits = c(min(band$Day), max(band$Day) + x_pad)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  ## "Photoconverted SIGNAL", not "protein": without a photobleaching control
  labs(x = "Day", y = expression(Photoconverted~signal~retained~(F/F[0]))) +
  ## The direct labels ARE the key, so the colour legend is redundant here.
  guides(linetype = "none", fill = "none", colour = "none") +
  theme_pub()

## ---- PANEL B ----------------------------------------------------------------
yv <- c(per_slice_K$half_life, rates$hl_lo, rates$hl_hi)
pB <- ggplot(per_slice_K, aes(Group, half_life)) +
  geom_errorbar(data = rates, aes(Group, ymin = hl_lo, ymax = hl_hi), width = 0.14,
                linewidth = 0.3, colour = "grey45", inherit.aes = FALSE) +
  geom_errorbar(data = rates, aes(Group, ymin = half_life, ymax = half_life),
                width = 0.26, linewidth = 0.45, colour = "grey20", inherit.aes = FALSE) +
  geom_point(aes(fill = Group, shape = Group), colour = "white", stroke = 0.3, size = 2,
             position = position_jitter(width = 0.06, height = 0, seed = 1)) +
  sig_layer(sig_pairs, yv) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_shape_manual(values = SHAPES, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  labs(x = NULL, y = "Half-life (days)") +
  theme_pub() + x_group_scale(narrow_panel = TRUE)

## ---- PANEL C ----------------------------------------------------------------
fitted_obj <- readRDS(file.path(RESULTS, "fitted_objects.rds"))
gfp_av <- anova(fitted_obj$m_gfp)
gfp_p <- if ("Day:Group" %in% rownames(gfp_av)) gfp_av["Day:Group", "Pr(>F)"] else NA

## The green pool is NOT purely new synthesis: after photoconversion ~30% of
.gfp_ylab <- expression(Unconverted~signal~(Delta*F/F[0]~vs~"pre-conversion"))

if (FACTORIAL) {
  ## Same layout as the main decay figure: one small panel per background,
  gfp_f <- gfp_model %>%
    mutate(Background = factor(GROUP_BG[as.character(Group)], levels = BACKGROUNDS),
           Construct = factor(GROUP_CONS[as.character(Group)], levels = c("Control", "HSPB1")))
  gfp_s <- gfp_f %>% group_by(Background, Construct, Day) %>%
    summarise(n = sum(!is.na(GFP_dFF0)), mean = mean(GFP_dFF0, na.rm = TRUE),
              sem = sd(GFP_dFF0, na.rm = TRUE) / sqrt(n), .groups = "drop")
  pC <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
    geom_line(data = gfp_f, aes(Day, GFP_dFF0, group = SliceID, colour = Construct),
              linewidth = 0.2, alpha = 0.25) +
    geom_errorbar(data = gfp_s, aes(Day, ymin = mean - sem, ymax = mean + sem,
                                    colour = Construct),
                  width = 0, linewidth = 0.4) +
    geom_line(data = gfp_s, aes(Day, mean, colour = Construct), linewidth = 0.7) +
    geom_point(data = gfp_s, aes(Day, mean, fill = Construct), shape = 21,
               colour = "white", size = 1.8, stroke = 0.3) +
    facet_wrap(~Background, nrow = 1) +
    scale_colour_manual(values = CONS_PAL, name = NULL) +
    scale_fill_manual(values = CONS_PAL, name = NULL) +
    scale_x_continuous(breaks = sort(unique(gfp_model$Day)), labels = day_labels) +
    labs(x = "Day", y = .gfp_ylab,
         subtitle = paste0("Group x day interaction: P = ",
                           format.pval(gfp_p, digits = 2, eps = 1e-4))) +
    guides(colour = "none", fill = "none") +
    theme_pub()
} else {
  gfp_s <- gfp_model %>% group_by(Group, Day) %>%
    summarise(n = sum(!is.na(GFP_dFF0)), mean = mean(GFP_dFF0, na.rm = TRUE),
              sem = sd(GFP_dFF0, na.rm = TRUE) / sqrt(n), .groups = "drop")
  pC <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
    geom_line(data = gfp_model, aes(Day, GFP_dFF0, group = SliceID, colour = Group),
              linewidth = 0.2, alpha = 0.22) +
    geom_errorbar(data = gfp_s, aes(Day, ymin = mean - sem, ymax = mean + sem, colour = Group),
                  width = 0, linewidth = 0.4) +
    geom_line(data = gfp_s, aes(Day, mean, colour = Group, linetype = Group), linewidth = 0.7) +
    geom_point(data = gfp_s, aes(Day, mean, fill = Group, shape = Group),
               colour = "white", size = 1.8, stroke = 0.3) +
    annotate("text", x = min(gfp_s$Day), y = max(gfp_s$mean + gfp_s$sem),
             label = paste0("group x day ", stars(gfp_p)), hjust = 0, vjust = 1,
             size = 2.6, colour = "grey25") +
    scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
    scale_shape_manual(values = SHAPES) + scale_linetype_manual(values = LINES) +
    scale_x_continuous(breaks = sort(unique(gfp_model$Day)), labels = day_labels) +
    labs(x = "Day", y = .gfp_ylab) +
    guides(linetype = "none", fill = "none", shape = "none", colour = "none") +
    theme_pub()
}

## ---- PANEL D ----------------------------------------------------------------
ref <- GROUPS[1]
fp <- pairs_df[pairs_df$group1 == ref, , drop = FALSE]
if (!nrow(fp)) fp <- pairs_df
fp$contrast <- paste(fp$group2, "vs", fp$group1)
fp$label <- sprintf("%.2f (%.2f-%.2f) %s", fp$ratio, fp$lo, fp$hi,
                    vapply(fp$p_adj, stars, character(1)))
## Labels sit above their point, so the x-range must leave room on BOTH sides
## or a label on an extreme estimate gets clipped at the panel edge.
lab_pad <- diff(range(c(fp$lo, fp$hi, 1))) * 0.55
xmin <- min(0.85, min(fp$lo) - lab_pad * 0.5)
xmax <- max(1.2, max(fp$hi) + lab_pad * 0.5)

pD <- ggplot(fp, aes(ratio, contrast)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.35, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.09,
                linewidth = 0.45, colour = "grey25") +
  geom_point(size = 2.8, shape = 21, fill = PAL[2], colour = "white", stroke = 0.4) +
  geom_text(aes(label = label), vjust = -1.3, size = 2.5, colour = "grey20") +
  scale_x_continuous(limits = c(xmin, xmax)) +
  ## Tighter vertical padding when there is only one contrast, otherwise a
  ## single row floats in a sea of white space when saved standalone.
  scale_y_discrete(expand = expansion(add = if (nrow(fp) == 1) 0.32 else 0.6)) +
  labs(x = "Turnover rate ratio", y = NULL) +
  theme_pub() + theme(axis.text.y = element_text(size = rel(0.95)))

## ---- assemble ---------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none")
save_fig(fig, "Figure4_MAIN_publication", W2, 150)

## ---- each panel on its own, for assembling a custom figure -----------------
panel_h <- 68 + 3 * max(0, N_GROUP - 2)          # taller as groups are added
cat_w <- fig_width(W1, 17)                        # categorical panels widen
line_w <- if (many) W2 * 0.62 else W1             # curve panels need label room

## Panel C carries no group labels inside the composite, because panel A
.leg_rows <- if (FACTORIAL) 1 else ceiling(N_GROUP / 2)
pC_standalone <- pC +
  guides(colour = guide_legend(nrow = .leg_rows, byrow = TRUE),
         fill = guide_legend(nrow = .leg_rows, byrow = TRUE)) +
  theme(legend.position = "bottom")
.pC_w <- if (FACTORIAL) W2 else line_w

single_panels <- list(
  list(p = pA, name = "Figure4A_turnover_curves",  w = line_w, h = panel_h),
  list(p = pB, name = "Figure4B_halflife",         w = cat_w,  h = panel_h),
  list(p = pC_standalone, name = "Figure4C_new_synthesis",
       w = .pC_w, h = panel_h + 8 * .leg_rows),
  list(p = pD, name = "Figure4D_effect_size",
       w = max(W1, W2 * 0.55), h = 26 + 11 * nrow(fp)))

for (sp in single_panels)
  save_fig(sp$p + theme(plot.margin = margin(3, 3, 3, 3)), sp$name, sp$w, sp$h)

## ---- the figure legend, as manuscript text ---------------------------------
n_by <- dendra %>% distinct(Group, Slice) %>% count(Group)
n_txt <- paste(sprintf("%s, n = %d slices", n_by$Group, n_by$n), collapse = "; ")
plateau_val <- plogis(unname(beta[p_name]))
hl_txt <- paste(sprintf("%s %.1f days (95%% CI %.1f-%.1f)",
                        rates$Group, rates$half_life, rates$hl_lo, rates$hl_hi),
                collapse = "; ")
ratio_txt <- paste(sprintf("%s versus %s, rate ratio %.2f (95%% CI %.2f-%.2f, P = %s)",
                           fp$group2, fp$group1, fp$ratio, fp$lo, fp$hi,
                           format.pval(fp$p_adj, digits = 2, eps = 1e-4)),
                   collapse = "; ")

legend_txt <- c(
  "FIGURE 4. Turnover of the photoconverted Dendra2 pool.",
  "",
  paste0("(A) Fraction of photoconverted (red) signal retained over the chase. ",
         "Thin lines and small points, individual slices; open symbols at day ",
         BASELINE_RFP_DISPLAY, " post-conversion denote the structural normalisation ",
         "anchor, which defines ",
         "fraction = 1 and was excluded from the likelihood; thick lines, population ",
         "fit of the constrained nonlinear mixed-effects model; shaded bands, 95% ",
         "confidence intervals of the fit; dashed horizontal line, shared plateau (",
         sprintf("%.2f", plateau_val), ")."),
  "",
  paste0("(B) Half-life of the decaying component per slice. Points, individual ",
         "slices; horizontal bar, model-based group estimate; whiskers, 95% CI. ",
         "Half-lives: ", hl_txt, "."),
  "",
  paste0("(C) Recovery of the unconverted (green) signal, expressed as dF/F0. ",
         "Thin lines, individual slices; symbols, group mean +/- s.e.m. ",
         "Annotation gives the group x day interaction from the linear mixed model ",
         "(P = ", format.pval(gfp_p, digits = 2, eps = 1e-4), ")."),
  "",
  paste0("(D) Turnover rate ratio with 95% CI. ", ratio_txt, "."),
  "",
  paste0(n_txt, ". Rate constants were compared by a prespecified maximum-likelihood ",
         "likelihood-ratio test on the group effect on log(K); pairwise contrasts are ",
         "Holm-adjusted. Asterisks: *P < 0.05, **P < 0.01, ***P < 0.001; n.s., not ",
         "significant. Exact P values are given above and in the accompanying tables."),
  "",
  "NOTE: the slice is the repeated-measures unit. Report the number of animals or",
  "preparations separately; slices sharing a preparation are not independent.")

writeLines(legend_txt, file.path(RESULTS, "Figure4_LEGEND.txt"))
message("  legend: Figure4_LEGEND.txt")

save_table(rates, "16_main_figure_group_rates")
save_table(pairs_df, "16_main_figure_rate_ratios")
message("\nComposite publication figure written to ", normalizePath(FIGURES))
