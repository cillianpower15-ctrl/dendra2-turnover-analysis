## 02_figures.R  --  publication figures. Run AFTER 01_statistics.R.
##   Rscript 02_figures.R

source("00_setup.R")

fitted_obj <- readRDS(file.path(RESULTS, "fitted_objects.rds"))
hl <- fitted_obj$hl; auc <- fitted_obj$auc; primary <- fitted_obj$primary

message("Building figures for ", N_GROUP, " groups...")

grp_scale <- list(
  scale_colour_manual(values = PAL, name = NULL),
  scale_fill_manual(values = PAL, name = NULL),
  scale_shape_manual(values = SHAPES, name = NULL),
  scale_linetype_manual(values = LINES, name = NULL))

fmt_p <- function(p) {
  if (is.na(p)) return("italic(P)==NA")
  if (p < 0.001) "italic(P)<0.001" else sprintf("italic(P)==%.3f", p)
}

## FIGURE 1 -- decay trajectories
## Individual slices thin in the background so the reader sees the real spread
## and the n; group mean +/- SEM on top.
summ <- function(df, col) {
  df %>% group_by(Group, Day) %>%
    summarise(n = sum(!is.na(.data[[col]])),
              mean = mean(.data[[col]], na.rm = TRUE),
              sem  = sd(.data[[col]], na.rm = TRUE) / sqrt(n), .groups = "drop")
}

panel_decay <- function(df, col, ylab, title, subtitle, baseline, include_baseline) {
  d <- if (include_baseline) df %>% filter(Day >= baseline) else df %>% filter(Day > baseline)
  s <- summ(d, col)
  brk <- sort(unique(d$Day))
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
    geom_line(data = d, aes(Day, .data[[col]], group = SliceID, colour = Group),
              linewidth = 0.22, alpha = 0.25) +
    geom_errorbar(data = s, aes(Day, ymin = mean - sem, ymax = mean + sem, colour = Group),
                  width = 0, linewidth = 0.4) +
    geom_line(data = s, aes(Day, mean, colour = Group, linetype = Group), linewidth = 0.65) +
    geom_point(data = s, aes(Day, mean, fill = Group, shape = Group),
               colour = "white", size = 1.9, stroke = 0.35) +
    grp_scale +
    scale_x_continuous(breaks = brk, labels = day_labels,
                       expand = expansion(mult = c(0.03, 0.04))) +
    labs(x = "Day", y = ylab, title = title, subtitle = subtitle) +
    theme_pub()
}

f1a <- panel_decay(dendra, "RFP_dFF0", expression(RFP~Delta*F/F[0]),
                   "Red - photoconverted pool (chase)",
                   "Decay of protein present at photoconversion",
                   BASELINE_RFP, TRUE)
## Not "new synthesis": ~30% of green survives photoconversion, so the green
## signal is residual unconverted protein PLUS anything newly made. Only the
## rise above the post-conversion level is unambiguously new.
f1b <- panel_decay(dendra, "GFP_dFF0", expression(Unconverted~signal~(Delta*F/F[0])),
                   "Green - unconverted pool (residual + new synthesis)",
                   "Recovery above the post-conversion level indicates new protein",
                   BASELINE_GFP, FALSE)

n_txt <- dendra %>% distinct(Group, Slice) %>% count(Group) %>%
  mutate(lab = paste0(Group, " n=", n)) %>% pull(lab) %>% paste(collapse = ", ")

fig1_w <- W2
fig1 <- (f1a | f1b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
    caption = wrap_text(paste0("Thin lines, individual slices; points, group mean +/- SEM. ", n_txt, "."),
                        fig1_w)) &
  theme_pub() & theme(legend.position = "bottom")
save_fig(fig1, "Figure1_decay_trajectories", fig1_w,
         78 + 4 * max(0, N_GROUP - 2) + 5 * ceiling(nchar(n_txt) / (fig1_w / 1.7)))

## FIGURE 2 -- per-replicate metrics, one point per slice
panel_points <- function(tbl, metric, ylab, title) {
  d <- tbl %>% select(Group, value = all_of(metric)) %>% filter(!is.na(value))
  s <- d %>% group_by(Group) %>%
    summarise(mean = mean(value), sd = sd(value), n = n(),
              lo = mean - qt(0.975, n - 1) * sd / sqrt(n),
              hi = mean + qt(0.975, n - 1) * sd / sqrt(n), .groups = "drop")
  rows <- primary[primary$metric == title, , drop = FALSE]
  rng <- range(c(d$value, s$lo, s$hi)); pad <- diff(rng) * 0.18

  ## In a factorial design, colour follows the CONSTRUCT (Control always blue,
  if (FACTORIAL) {
    d$Construct <- GROUP_CONS[as.character(d$Group)]
    d$Background <- GROUP_BG[as.character(d$Group)]
    pts <- list(
      geom_point(aes(fill = Construct, shape = Background), colour = "white",
                 stroke = 0.3, size = 2.0,
                 position = position_jitter(width = 0.07, height = 0, seed = 1)),
      scale_fill_manual(values = CONS_PAL, guide = "none"),
      scale_shape_manual(values = BG_SHAPES, guide = "none"))
  } else {
    pts <- list(
      geom_point(aes(fill = Group, shape = Group), colour = "white", stroke = 0.3,
                 size = 2.0,
                 position = position_jitter(width = 0.07, height = 0, seed = 1)),
      scale_fill_manual(values = PAL, guide = "none"),
      scale_shape_manual(values = SHAPES, guide = "none"))
  }
  p <- ggplot(d, aes(Group, value)) +
    geom_errorbar(data = s, aes(Group, ymin = lo, ymax = hi), width = 0.16,
                  linewidth = 0.3, colour = "grey45", inherit.aes = FALSE) +
    geom_errorbar(data = s, aes(Group, ymin = mean, ymax = mean), width = 0.30,
                  linewidth = 0.45, colour = "grey20", inherit.aes = FALSE) +
    pts +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub() + x_group_scale(narrow_panel = TRUE)

  if (TWO_GROUP && nrow(rows)) {
    ytop <- rng[2] + pad * 0.85
    p <- p +
      annotate("segment", x = 1, xend = 2, y = ytop, yend = ytop,
               linewidth = 0.3, colour = "grey30") +
      annotate("text", x = 1.5, y = ytop + pad * 0.10, label = fmt_p(rows$rank_p[1]),
               parse = TRUE, size = 2.5, colour = "grey20", vjust = 0)
  } else if (nrow(rows) && "omnibus_kruskal_p" %in% names(rows)) {
    p <- p + labs(subtitle = sprintf("Kruskal-Wallis P = %.3g", rows$omnibus_kruskal_p[1]))
  }
  p
}

f2a <- panel_points(hl,  "half_life", "Half-life (days)", "RFP half-life (days)")
f2b <- panel_points(auc, "RFP_AUC",   "AUC",              "RFP AUC")
f2c <- panel_points(auc, "GFP_AUC",   "AUC",              "GFP AUC")

cap2 <- if (TWO_GROUP) {
  "Supporting summaries: points, individual slices; bar, mean; whiskers, 95% CI. P, rank-based test with Holm adjustment across the three per-slice metrics."
} else if (FACTORIAL) {
  paste("Points, individual slices (blue = Control, orange = HSPB1); bar, mean; whiskers, 95% CI of the mean.",
        "Omnibus Kruskal-Wallis P shown; pairwise comparisons are in tables/08_group_comparisons.csv.")
} else {
  paste("Points, individual slices; bar, mean; whiskers, 95% CI of the mean.",
        "Omnibus P shown; all pairwise comparisons are in tables/08_group_comparisons.csv.")
}

## NOTE: each panel already carries a complete theme_pub()+x_group_scale().
fig2_w <- W2
fig2 <- (f2a | f2b | f2c) +
  plot_annotation(tag_levels = "A", caption = wrap_text(cap2, fig2_w)) &
  theme(plot.caption = element_text(size = rel(0.85), colour = "grey35", hjust = 0))
save_fig(fig2, "Figure2_per_replicate_metrics", fig2_w,
         72 + 4 * max(0, N_GROUP - 2) + ifelse(N_GROUP > 2, 14, 0))

## FIGURE S1 -- model diagnostics (supplement)
diag_df <- function(m, label)
  data.frame(model = label, fitted = fitted(m), resid = resid(m),
             std = resid(m) / sd(resid(m)))
dd <- bind_rows(diag_df(fitted_obj$m_rfp, "RFP model"), diag_df(fitted_obj$m_gfp, "GFP model"))

s1a <- ggplot(dd, aes(fitted, resid)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_point(size = 1.1, alpha = 0.65, colour = "grey25") +
  facet_wrap(~model, scales = "free") +
  labs(x = "Fitted value", y = "Residual", title = "Residuals vs fitted") + theme_pub()
s1b <- ggplot(dd, aes(sample = std)) +
  stat_qq(size = 1.1, alpha = 0.65, colour = "grey25") +
  stat_qq_line(linewidth = 0.35, colour = PAL[1]) +
  facet_wrap(~model, scales = "free") +
  labs(x = "Theoretical quantile", y = "Standardised residual", title = "Normal Q-Q") + theme_pub()
save_fig((s1a / s1b) + plot_annotation(tag_levels = "A") & theme_pub(),
         "FigureS1_model_diagnostics", W2, 110)

hl_long <- hl %>%
  select(Group, Slice, `Free plateau (primary)` = half_life,
         `No plateau` = half_life_noplateau,
         `Log-linear (workbook)` = half_life_loglinear) %>%
  pivot_longer(-c(Group, Slice), names_to = "method", values_to = "half_life") %>%
  mutate(method = factor(method, levels = c("Free plateau (primary)", "No plateau",
                                            "Log-linear (workbook)")))
hl_summ <- hl_long %>% group_by(method, Group) %>%
  summarise(mean = mean(half_life, na.rm = TRUE), .groups = "drop")

figS2 <- ggplot(hl_long, aes(Group, half_life)) +
  geom_errorbar(data = hl_summ, aes(Group, ymin = mean, ymax = mean), width = 0.32,
                linewidth = 0.45, colour = "grey20", inherit.aes = FALSE) +
  geom_point(aes(fill = Group, shape = Group), colour = "white", stroke = 0.3, size = 1.9,
             position = position_jitter(width = 0.07, height = 0, seed = 1)) +
  facet_wrap(~method) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_shape_manual(values = SHAPES, guide = "none") +
  labs(x = NULL, y = "RFP half-life (days)",
       title = "Half-life depends on the decay model; the group ordering usually does not") +
  theme_pub() + x_group_scale()
save_fig(figS2, "FigureS2_halflife_model_sensitivity", W2, 72 + 6 * max(0, N_GROUP - 2))

message("\nFigures written to ", normalizePath(FIGURES))
