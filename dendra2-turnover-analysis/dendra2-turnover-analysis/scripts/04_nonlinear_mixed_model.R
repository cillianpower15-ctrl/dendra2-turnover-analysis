## 04_nonlinear_mixed_model.R -- primary test of photoconverted RFP turnover.

source("00_setup.R")
suppressPackageStartupMessages(library(nlme))

.log <- file.path(RESULTS, "nonlinear_model_log.txt")
con <- file(.log, open = "wt"); sink(con, split = TRUE)
hdr <- function(x) cat("\n\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

hdr("PRIMARY ANALYSIS: CONSTRAINED NONLINEAR MIXED-EFFECTS RFP MODEL")
describe_run()
cat("
Model: fraction retained = (1 - plogis(logitP)) * exp(-exp(logK) * t)
                           + plogis(logitP)

  t       = Day - ", BASELINE_RFP, "
  K       = exp(logK), constrained above zero
  P       = plogis(logitP), constrained between zero and one
  half-life of the decaying component = ln(2)/K
  slice random effect = logK ~ 1 | SliceID

The post-conversion baseline defines fraction(0)=1 mathematically. It is a
structural normalisation point and is EXCLUDED from the likelihood; only days
strictly after ", BASELINE_RFP, " are fitted.

M1: group-shared logK and P
M2: group-specific logK, shared P       [PRESPECIFIED PRIMARY COMPARISON]
M3: group-specific logK and P           [secondary]
All fixed-effect comparisons use maximum likelihood (ML).
Preparation handling: ", PREP_STATUS, "
")

dn <- decay_fit %>% as.data.frame()
if (nrow(dn) < 10 || any(table(dn$SliceID) < 3))
  stop("The nonlinear model requires at least three post-baseline measurements per slice.")
if (any(!is.finite(dn$frac)) || any(dn$frac <= 0))
  stop("Post-baseline RFP fractions must be finite and greater than zero.")
grp <- groupedData(frac ~ t | SliceID, data = dn)

parameter_formula <- function(parameter, group_effect = FALSE) {
  rhs <- c(if (group_effect) "Group", if (ADJUST_FOR_PREP) "PreparationID")
  ## c() of all-NULL returns NULL, and reformulate(NULL) errors out. An empty
  ## right-hand side must be "1", which gives "<parameter> ~ 1". character(0)
  ## is only accepted by reformulate() from R 4.4.0 onward.
  if (!length(rhs)) rhs <- "1"
  reformulate(rhs, response = parameter)
}
f_logK_shared <- parameter_formula("logK", FALSE)
f_logK_group <- parameter_formula("logK", TRUE)
f_logitP_shared <- reformulate("1", response = "logitP")
f_logitP_group <- reformulate("Group", response = "logitP")

start_for <- function(fixed, data) {
  bases <- c(logK = log(0.15), logitP = qlogis(0.25))
  unlist(lapply(names(fixed), function(p) {
    mm <- model.matrix(delete.response(terms(fixed[[p]])), data = data)
    c(unname(bases[p]), rep(0, ncol(mm) - 1))
  }), use.names = FALSE)
}

ctrl <- nlmeControl(maxIter = 300, msMaxIter = 300, pnlsMaxIter = 30,
                    pnlsTol = 0.05, returnObject = TRUE)
fit_model <- function(fixed, data = grp) {
  try(nlme(frac ~ (1 - plogis(logitP)) * exp(-exp(logK) * t) + plogis(logitP),
           data = data, fixed = fixed, random = logK ~ 1 | SliceID,
           start = start_for(fixed, as.data.frame(data)), method = "ML",
           na.action = na.fail, control = ctrl), silent = TRUE)
}

fixed_m1 <- list(logK = f_logK_shared, logitP = f_logitP_shared)
fixed_m2 <- list(logK = f_logK_group,  logitP = f_logitP_shared)
fixed_m3 <- list(logK = f_logK_group,  logitP = f_logitP_group)
m1 <- fit_model(fixed_m1)
m2 <- fit_model(fixed_m2)
m3 <- fit_model(fixed_m3)

show_model <- function(m, label) {
  cat("\n-- ", label, " --\n", sep = "")
  if (inherits(m, "try-error")) {
    cat("FAILED TO CONVERGE\n"); return(invisible(NULL))
  }
  print(summary(m)$tTable, digits = 4)
  cat(sprintf("AIC %.2f; BIC %.2f; logLik %.2f\n", AIC(m), BIC(m), logLik(m)))
}
hdr("1. MODEL FITS")
show_model(m1, "M1: shared rate and plateau")
show_model(m2, "M2: group-specific rate, shared plateau")
show_model(m3, "M3: group-specific rate and plateau")

if (inherits(m1, "try-error") || inherits(m2, "try-error")) {
  stop("Primary nonlinear M1 or M2 did not converge. Review nonlinear_model_log.txt; ",
       "do not report a primary turnover-rate p-value from this run.")
}

hdr("2. PRESPECIFIED PRIMARY MODEL COMPARISON")
cmp12 <- anova(m1, m2)
print(cmp12)
comparison_table <- data.frame(
  role = "PRIMARY",
  comparison = "Does turnover rate differ between groups? (M1 vs M2, log K)",
  likelihood_ratio = cmp12[["L.Ratio"]][2],
  df = cmp12[["df"]][2] - cmp12[["df"]][1],
  p = cmp12[["p-value"]][2])
if (!inherits(m3, "try-error")) {
  cmp23 <- anova(m2, m3)
  comparison_table <- rbind(comparison_table, data.frame(
    role = "SECONDARY",
    comparison = "Does the plateau differ between groups? (M2 vs M3, logit P)",
    likelihood_ratio = cmp23[["L.Ratio"]][2],
    df = cmp23[["df"]][2] - cmp23[["df"]][1],
    p = cmp23[["p-value"]][2]))
}
save_table(comparison_table, "14_nonlinear_model_comparisons")

## Adjusted group estimates and contrasts from M2.
beta <- fixef(m2)
V <- vcov(m2)
df_candidates <- summary(m2)$fixDF$X
df_resid <- suppressWarnings(min(df_candidates[is.finite(df_candidates)]))
if (!is.finite(df_resid)) df_resid <- max(1, nrow(dn) - length(beta))
tcrit <- qt(0.975, df_resid)

logK_vector <- function(group) {
  grid <- if (ADJUST_FOR_PREP) {
    expand.grid(Group = factor(group, levels = GROUPS),
                PreparationID = levels(dendra$PreparationID))
  } else data.frame(Group = factor(group, levels = GROUPS))
  mm <- model.matrix(delete.response(terms(f_logK_group)), data = grid)
  x <- colMeans(mm)
  out <- setNames(rep(0, length(beta)), names(beta))
  expected <- paste0("logK.", names(x))
  hit <- match(expected, names(out))
  if (anyNA(hit)) stop("Could not map nonlinear logK design columns to coefficients.")
  out[hit] <- x
  out
}

rate_rows <- lapply(GROUPS, function(g) {
  x <- logK_vector(g)
  est <- as.numeric(sum(x * beta))
  se <- as.numeric(sqrt(t(x) %*% V %*% x))
  lo <- est - tcrit * se; hi <- est + tcrit * se
  data.frame(Group = g, logK = est, SE_logK = se,
             K_per_day = exp(est), K_CI_low = exp(lo), K_CI_high = exp(hi),
             half_life_days = log(2) / exp(est),
             half_life_CI_low = log(2) / exp(hi),
             half_life_CI_high = log(2) / exp(lo))
})
group_rates <- bind_rows(rate_rows)

hdr("3. ADJUSTED GROUP RATE CONSTANTS")
print(group_rates, digits = 4, row.names = FALSE)
save_table(group_rates, "12_nonlinear_group_rates")

## nlme names a parameter with an intercept-only formula plain ("logitP"), but
## names it "logitP.(Intercept)" as soon as that parameter has predictors.
## Accept either so the plateau is found in both M2 and M3 parameterisations.
p_index <- match("logitP", names(beta))
if (is.na(p_index)) p_index <- match("logitP.(Intercept)", names(beta))
if (is.na(p_index)) stop("Could not locate the shared plateau coefficient.")
p_logit <- unname(beta[p_index])
p_se <- sqrt(V[p_index, p_index])
plateau_table <- data.frame(
  logitP = p_logit, SE_logitP = p_se, plateau = plogis(p_logit),
  plateau_CI_low = plogis(p_logit - tcrit * p_se),
  plateau_CI_high = plogis(p_logit + tcrit * p_se))
save_table(plateau_table, "12_nonlinear_shared_plateau")

hdr("4. HOLM-ADJUSTED PAIRWISE RATE RATIOS")
pairs <- combn(GROUPS, 2, simplify = FALSE)
pairwise_rates <- bind_rows(lapply(pairs, function(p) {
  x <- logK_vector(p[1]) - logK_vector(p[2])
  est <- as.numeric(sum(x * beta))
  se <- as.numeric(sqrt(t(x) %*% V %*% x))
  stat <- est / se
  p_raw <- 2 * pt(-abs(stat), df_resid)
  data.frame(group1 = p[1], group2 = p[2], log_rate_ratio = est, SE = se,
             df = df_resid, t = stat, rate_ratio = exp(est),
             ratio_CI_low = exp(est - tcrit * se), ratio_CI_high = exp(est + tcrit * se),
             p_raw = p_raw)
}))
pairwise_rates$p_holm <- p.adjust(pairwise_rates$p_raw, "holm")
print(pairwise_rates, digits = 4, row.names = FALSE)
cat("\nrate_ratio > 1 means faster fitted turnover in group1 than group2.\n")
save_table(pairwise_rates, "11_nonlinear_pairwise_rate_ratios")

hdr("5. ROBUSTNESS: LEAVE ONE SLICE OUT")
loso <- bind_rows(lapply(unique(dn$SliceID), function(s) {
  d2 <- droplevels(filter(dn, SliceID != s))
  if (n_distinct(d2$Group) < N_GROUP)
    return(data.frame(dropped = s, p = NA_real_, status = "group lost"))
  g2 <- groupedData(frac ~ t | SliceID, data = d2)
  a <- fit_model(fixed_m1, g2); b <- fit_model(fixed_m2, g2)
  if (inherits(a, "try-error") || inherits(b, "try-error"))
    return(data.frame(dropped = s, p = NA_real_, status = "fit failed"))
  data.frame(dropped = s, p = anova(a, b)[["p-value"]][2], status = "fit")
}))
print(loso[order(-loso$p), ], digits = 4, row.names = FALSE)
save_table(loso, "13_nonlinear_leave_one_out")

res <- resid(m2, type = "normalized")
shapiro_p <- if (length(res) >= 3 && length(res) <= 5000) shapiro.test(res)$p.value else NA_real_
diagnostics <- data.frame(
  model = "M2", n_observations = nobs(m2), n_slices = n_distinct(dn$SliceID),
  AIC = AIC(m2), BIC = BIC(m2), logLik = as.numeric(logLik(m2)),
  residual_SD = sd(res), residual_Shapiro_p = shapiro_p,
  preparation_adjusted = ADJUST_FOR_PREP,
  structural_baseline_in_likelihood = FALSE)
save_table(diagnostics, "15_nonlinear_model_diagnostics")

## ---- publication figure ----------------------------------------------------
pred <- bind_rows(lapply(seq_len(nrow(group_rates)), function(i) {
  tt <- seq(0, max(dn$t), length.out = 200)
  data.frame(Group = group_rates$Group[i], t = tt,
             Day = tt + BASELINE_RFP,
             frac = (1 - plateau_table$plateau) * exp(-group_rates$K_per_day[i] * tt) +
                    plateau_table$plateau)
}))
pred$Group <- factor(pred$Group, levels = GROUPS)
anchor <- data.frame(Group = factor(GROUPS, levels = GROUPS),
                     Day = BASELINE_RFP, frac = 1)
primary_p <- comparison_table$p[comparison_table$role == "PRIMARY"][1]
subtitle <- if (TWO_GROUP) {
  sprintf("Adjusted K %.3f vs %.3f /day; rate ratio %.2f (Holm P = %.3g); omnibus ML LRT P = %.3g",
          group_rates$K_per_day[1], group_rates$K_per_day[2],
          pairwise_rates$rate_ratio[1], pairwise_rates$p_holm[1], primary_p)
} else sprintf("Prespecified omnibus ML likelihood-ratio test across %d groups: P = %.3g",
               N_GROUP, primary_p)

figN <- ggplot() +
  geom_line(data = dn, aes(Day, frac, group = SliceID, colour = Group),
            linewidth = 0.22, alpha = 0.25) +
  geom_point(data = dn, aes(Day, frac, colour = Group), size = 0.65, alpha = 0.45) +
  geom_point(data = anchor, aes(Day, frac, colour = Group), shape = 1, size = 1.8) +
  geom_line(data = pred, aes(Day, frac, colour = Group, linetype = Group), linewidth = 0.85) +
  scale_colour_manual(values = PAL) + scale_linetype_manual(values = LINES) +
  scale_x_continuous(breaks = sort(unique(c(BASELINE_RFP, dn$Day)))) +
  labs(x = "Day", y = "Fraction of photoconverted signal retained",
       title = "Constrained nonlinear mixed-effects model of RFP turnover",
       subtitle = wrap_text(subtitle, fig_width(W1 + 20, 20)),
       caption = wrap_text(paste("Thin lines/points: included slices.",
                                 "Thick lines: adjusted group fits.",
                                 "Open t=0 points are structural anchors and were not fitted."),
                           fig_width(W1 + 20, 20))) +
  theme_pub() + guides(linetype = "none")
if (is.finite(plateau_table$plateau))
  figN <- figN + geom_hline(yintercept = plateau_table$plateau, linetype = "22",
                            linewidth = 0.3, colour = "grey45")
save_fig(figN, "Figure3_nonlinear_mixed_model", fig_width(W1 + 20, 20), 92)

hdr("DONE")
cat("Primary comparison p = ", format(primary_p, digits = 4), "\n", sep = "")
cat("Log: ", normalizePath(.log), "\n", sep = "")
sink(); close(con)
message("\nDone. Log written to ", .log)
