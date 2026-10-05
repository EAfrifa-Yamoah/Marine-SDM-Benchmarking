# =====================================================================
# 10_theory_analysis.R
# Analyses motivated by the theoretical foundation, computed from the
# existing pilot fits. Each block names the result it tests.
#
#   A  Theorem 1 / Proposition 1 identities verified against the
#      benchmark's own co observed metrics (checks the scoring code)
#   B  Harmonisation restructured: exact vs empirical kind, interdecile
#      normalisation (Corollary 2), prevalence conditioning (Theorem 1b),
#      within context residual and the pooling bound (Theorem 3, Cor. 3)
#   C  Cluster bootstrap over species within survey: class stability
#   D  Proposition 3: predicted vs observed directional asymmetry
#   E  Theorem 2(iii): calibration slope as the residual of AUC -> Brier
#   F  Proposition 1: MAE - 2 Brier as a calibration diagnostic
#   G  Corollary 5: Brier optimism against log sample size by class
#   H  Lemma 5: marginal spatial advantage against covariate count
#   I  fit count arithmetic
#
# Analyses that need refitted models (same fit contrast, near oracle
# intrinsic term, isotonic recalibration, fill distance and separation of
# each draw) are NOT here; they require a driver change and a rerun.
# =====================================================================
suppressMessages({ library(data.table); library(mgcv) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT")); source(file.path("R", "config.R"))   # design constants for the fit count
RES <- "results"; dir.create(RES, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")

d0 <- fread(file.path(RES, "fit_metrics.csv"))
d  <- d0[converged == TRUE & is.finite(AUC) & is.finite(R2) & is.finite(Brier)]
say("=== theory analysis over %d converged fits ===", nrow(d))

# ---------------------------------------------------------------------
# A. identities against the benchmark's own scores
# ---------------------------------------------------------------------
pi <- d$test_prev; n <- d$n_test; n1 <- round(n * pi); n0 <- n - n1
ds <- d[is.finite(Spearman)]
pis <- ds$test_prev; ns <- ds$n_test; n1s <- round(ns * pis); n0s <- ns - n1s
res_1c <- abs(ds$Spearman - (ds$AUC - 0.5) * sqrt(12 * n0s * n1s / (ns^2 - 1)))
ident <- data.table(
  result = c("Theorem 1a  Brier = RMSE^2",
             "Theorem 1b  R2 = 1 - Brier / (pi (1 - pi)), pi = evaluation prevalence",
             "Theorem 1c  Spearman = (AUC - 1/2) sqrt(12 n0 n1 / (n^2 - 1)), no tie form",
             "Proposition 1  MAE^2 <= Brier <= MAE"),
  n_fits = c(nrow(d), nrow(d), nrow(ds), nrow(d)),
  max_abs_dev = c(max(abs(d$Brier - d$RMSE^2)),
                  max(abs(d$R2 - (1 - d$Brier / (pi * (1 - pi))))),
                  max(res_1c),
                  max(pmax(d$MAE^2 - d$Brier, d$Brier - d$MAE, 0))),
  median_abs_dev = c(median(abs(d$Brier - d$RMSE^2)),
                     median(abs(d$R2 - (1 - d$Brier / (pi * (1 - pi))))),
                     median(res_1c), 0),
  n_dev_over_1e3 = c(sum(abs(d$Brier - d$RMSE^2) > 1e-3),
                     sum(abs(d$R2 - (1 - d$Brier / (pi * (1 - pi)))) > 1e-3),
                     sum(res_1c > 1e-3),
                     sum(d$MAE^2 > d$Brier + 1e-12 | d$Brier > d$MAE + 1e-12)))
fwrite(ident, file.path(RES, "table_identity_verification.csv"))
print(ident[, .(result, n_fits, max_abs_dev = signif(max_abs_dev, 3), n_dev_over_1e3)])
ties_by_method <- ds[, .(n = .N, n_tie_affected = sum(res_1c > 1e-3),
                         max_dev = signif(max(res_1c), 3)), by = method][order(-n_tie_affected)]
fwrite(ties_by_method, file.path(RES, "table_spearman_tie_effect_by_method.csv"))

# ---------------------------------------------------------------------
# B. harmonisation restructured by the theory
# ---------------------------------------------------------------------
METS <- c("AUC", "TSS", "Spearman", "R2", "RMSE", "Brier", "MAE")
pairs <- list(c("RMSE","Brier"), c("Brier","RMSE"),
              c("R2","RMSE"), c("RMSE","R2"), c("R2","Brier"), c("Brier","R2"),
              c("AUC","Spearman"), c("Spearman","AUC"),
              c("AUC","TSS"), c("TSS","AUC"), c("TSS","Spearman"),
              c("AUC","RMSE"), c("AUC","Brier"), c("AUC","R2"),
              c("MAE","Brier"), c("Brier","MAE"))
kind_of <- function(a, b) {
  err <- c("R2", "RMSE", "Brier")
  if (a %in% err && b %in% err) return("exact given prevalence")
  if (setequal(c(a, b), c("AUC", "Spearman"))) return("exact given score distribution")
  if (setequal(c(a, b), c("MAE", "Brier"))) return("empirical, bounded (Prop. 1)")
  "empirical"
}
grid_n <- 50L
fit_gam <- function(dat, from, to, xg) {
  m <- tryCatch(gam(as.formula(sprintf("%s ~ s(%s, k=5)", to, from)),
                    data = dat, method = "REML"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(pred = as.numeric(predict(m, setNames(data.frame(xg), from))),
       resid = residuals(m), fitted = fitted(m))
}
idr <- function(x) diff(quantile(x, c(.1, .9), names = FALSE))

harm_one <- function(dat, from, to, boot = FALSE) {
  sub <- dat[is.finite(get(from)) & is.finite(get(to))]
  if (nrow(sub) < 200) return(NULL)
  rng <- quantile(sub[[from]], c(.02, .98), names = FALSE)
  xg  <- seq(rng[1], rng[2], length.out = grid_n)
  ov  <- fit_gam(sub, from, to, xg); if (is.null(ov)) return(NULL)
  # shifts across inference target
  co <- fit_gam(sub[target == "conditional"], from, to, xg)
  ma <- fit_gam(sub[target == "marginal"],    from, to, xg)
  tshift <- if (!is.null(co) && !is.null(ma)) mean(abs(co$pred - ma$pred)) else NA_real_
  # shifts across survey (each survey on its own support)
  dss <- sort(unique(sub$dataset))
  P <- vapply(dss, function(s) {
    sd_ <- sub[dataset == s]; if (nrow(sd_) < 50) return(rep(NA_real_, grid_n))
    q <- quantile(sd_[[from]], c(.02, .98), names = FALSE)
    f <- fit_gam(sd_, from, to, xg); if (is.null(f)) return(rep(NA_real_, grid_n))
    p <- f$pred; p[xg < q[1] | xg > q[2]] <- NA_real_; p
  }, numeric(grid_n))
  ok <- stats::complete.cases(P)
  sshift <- if (sum(ok) >= 5) mean(apply(P[ok, , drop = FALSE], 1, function(r) mean(dist(r)))) else NA_real_
  # source weighted survey shift (Corollary 3 version): weight grid by pooled source density
  wgt <- if (sum(ok) >= 5) { h <- hist(sub[[from]], breaks = c(-Inf, xg[-1] - diff(xg)/2, Inf), plot = FALSE)$counts; h[ok] / sum(h[ok]) } else NA
  sshift_w <- if (sum(ok) >= 5) sum(wgt * apply(P[ok, , drop = FALSE], 1, function(r) mean(dist(r)))) else NA_real_
  # within context residual (target x survey), pooled variance
  wr <- sub[, {
    f <- if (.N >= 50) fit_gam(.SD, from, to, xg) else NULL
    if (is.null(f)) list(ss = NA_real_, nn = NA_integer_) else list(ss = sum(f$resid^2), nn = .N)
  }, by = .(target, dataset)]
  within_sd <- sqrt(sum(wr$ss, na.rm = TRUE) / sum(wr$nn, na.rm = TRUE))
  resid_sd <- sd(ov$resid)
  list(from = from, to = to, kind = kind_of(from, to),
       pearson = cor(sub[[from]], sub[[to]]),
       resid_sd = resid_sd,
       rel_resid_range = resid_sd / diff(range(sub[[to]])),
       rel_resid_idr   = resid_sd / idr(sub[[to]]),
       within_context_sd = within_sd,
       target_shift = tshift, survey_shift = sshift, survey_shift_weighted = sshift_w,
       pool_bound = sqrt(3/8) * sshift_w,          # Corollary 3, K = 4
       eta2 = var(ov$fitted) / var(sub[[to]]))
}
classify <- function(rel, tshift, sshift) {
  s <- max(c(tshift, sshift), na.rm = TRUE); if (!is.finite(s)) s <- Inf
  if (rel < 0.10 && s < 0.05) "reliable" else if (rel < 0.18 && s < 0.10) "context dependent" else "unreliable"
}
# The restructured table and the prevalence conditioning are the costly
# parts of this script at full scale (about 250 REML smooths on up to 1.2
# million fits). They are computed once and cached so that the resumable
# bootstrap below can be continued across several invocations without
# recomputing them. Delete the cache after any change to fit_metrics.rds.
cache_file <- file.path(RES, "theory_prebootstrap_cache.rds")
if (file.exists(cache_file)) {
  cc <- readRDS(cache_file); H <- cc$H; prev_cond <- cc$prev_cond
  say("loaded cached restructured harmonisation table (%d pairs) and prevalence conditioning", nrow(H))
} else {
H <- rbindlist(lapply(pairs, function(p) { r <- harm_one(d, p[1], p[2]); if (is.null(r)) NULL else as.data.table(r) }))
H[, binding_shift := ifelse(survey_shift > target_shift, "survey", "inference target")]
H[, class_range := mapply(classify, rel_resid_range, target_shift, survey_shift)]
H[, class_idr   := mapply(classify, rel_resid_idr,   target_shift, survey_shift)]
H[, pooled_exceeds_bound := resid_sd >= pool_bound]
fwrite(H, file.path(RES, "table_harmonisation_theory.csv"))
say("restructured harmonisation table written (bootstrap stability columns follow)")

# prevalence conditioning (Theorem 1b): shifts within prevalence deciles
d[, prev_bin := cut(test_prev, quantile(test_prev, seq(0, 1, .1)), include.lowest = TRUE, labels = FALSE)]
prev_cond <- rbindlist(lapply(list(c("R2","RMSE"), c("R2","Brier"), c("AUC","R2"), c("AUC","Brier")), function(p) {
  bb <- rbindlist(lapply(sort(unique(d$prev_bin)), function(b) {
    r <- harm_one(d[prev_bin == b], p[1], p[2]); if (is.null(r)) NULL else
      data.table(bin = b, target_shift = r$target_shift, survey_shift = r$survey_shift, rel_resid_idr = r$rel_resid_idr)
  }))
  data.table(from = p[1], to = p[2],
             target_shift_pooled = H[from == p[1] & to == p[2], target_shift],
             target_shift_within_prev = mean(bb$target_shift, na.rm = TRUE),
             survey_shift_pooled = H[from == p[1] & to == p[2], survey_shift],
             survey_shift_within_prev = mean(bb$survey_shift, na.rm = TRUE),
             n_bins = nrow(bb))
}))
saveRDS(list(H = H, prev_cond = prev_cond), cache_file)
}
fwrite(prev_cond, file.path(RES, "table_prevalence_conditioning.csv"))
say("\n--- prevalence conditioning (Theorem 1b): shifts pooled vs within prevalence deciles ---")
print(prev_cond[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

# ---------------------------------------------------------------------
# C. cluster bootstrap over species within survey (class stability)
# A lean classifier: only the seven fits that determine a class (overall,
# two targets, four surveys), with the faster GCV smoother.
# ---------------------------------------------------------------------
fit_fast <- function(dat, from, to, xg) {
  m <- tryCatch(gam(as.formula(sprintf("%s ~ s(%s, k=5)", to, from)), data = dat,
                    method = "GCV.Cp"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(pred = as.numeric(predict(m, setNames(data.frame(xg), from))), resid = residuals(m))
}
class_fast <- function(dat, from, to) {
  sub <- dat[is.finite(get(from)) & is.finite(get(to))]
  if (nrow(sub) < 200) return(c(NA, NA))
  rng <- quantile(sub[[from]], c(.02, .98), names = FALSE); xg <- seq(rng[1], rng[2], length.out = grid_n)
  ov <- fit_fast(sub, from, to, xg); if (is.null(ov)) return(c(NA, NA))
  co <- fit_fast(sub[target == "conditional"], from, to, xg); ma <- fit_fast(sub[target == "marginal"], from, to, xg)
  ts <- if (!is.null(co) && !is.null(ma)) mean(abs(co$pred - ma$pred)) else NA_real_
  P <- vapply(sort(unique(sub$dataset)), function(sv) {
    sd_ <- sub[dataset == sv]; if (nrow(sd_) < 50) return(rep(NA_real_, grid_n))
    q <- quantile(sd_[[from]], c(.02, .98), names = FALSE)
    f <- fit_fast(sd_, from, to, xg); if (is.null(f)) return(rep(NA_real_, grid_n))
    pp <- f$pred; pp[xg < q[1] | xg > q[2]] <- NA_real_; pp }, numeric(grid_n))
  ok <- stats::complete.cases(P)
  ss <- if (sum(ok) >= 5) mean(apply(P[ok, , drop = FALSE], 1, function(r) mean(dist(r)))) else NA_real_
  rs <- sd(ov$resid)
  c(classify(rs / diff(range(sub[[to]])), ts, ss), classify(rs / idr(sub[[to]]), ts, ss))
}
NB <- as.integer(Sys.getenv("SMD_BOOT_REPS", unset = "100"))
# wall clock budget for this invocation; replicates are checkpointed one at a
# time in boot_classes.rds, so the script can be re run until NB is reached
BOOT_SECS <- as.numeric(Sys.getenv("SMD_BOOT_SECS", unset = "200"))
boot_file <- file.path(RES, "boot_classes.rds")
boot_part <- file.path(RES, "boot_partial.rds")     # pairs finished in the replicate in progress
boot_cls <- if (file.exists(boot_file)) readRDS(boot_file) else list()
sp_by_ds <- d[, .(sp = unique(species)), by = dataset]
# Each replicate r draws its species resample from its own seed, so a
# replicate interrupted by the time budget is resumed with the same resample
# and only its pending pairs are refitted (at full scale one replicate is
# longer than a single invocation on a small machine).
t0 <- Sys.time(); elapsed <- function() as.numeric(Sys.time() - t0, units = "secs")
while (length(boot_cls) < NB && elapsed() < BOOT_SECS) {
  r <- length(boot_cls) + 1L
  set.seed(20260916 + r)
  samp <- sp_by_ds[, .(species = sample(sp, length(sp), replace = TRUE)), by = dataset]
  samp[, cid := seq_len(.N)]
  db <- merge(d, samp, by = c("dataset", "species"), allow.cartesian = TRUE)
  part <- if (file.exists(boot_part)) readRDS(boot_part) else list(rep = r, cls = list())
  if (!identical(part$rep, r)) part <- list(rep = r, cls = list())
  done_pairs <- vapply(part$cls, function(x) paste(x$from, x$to), "")
  for (p in pairs) {
    if (paste(p[1], p[2]) %in% done_pairs) next
    if (elapsed() >= BOOT_SECS) break
    k <- class_fast(db, p[1], p[2])
    part$cls[[length(part$cls) + 1]] <- data.table(from = p[1], to = p[2], class_range = k[1], class_idr = k[2])
    saveRDS(part, boot_part)
  }
  if (length(part$cls) < length(pairs)) {
    say("bootstrap replicate %d paused with %d of %d pairs done", r, length(part$cls), length(pairs)); break
  }
  boot_cls[[r]] <- rbindlist(part$cls)
  saveRDS(boot_cls, boot_file); unlink(boot_part)
  say("bootstrap replicate %d complete (%.1f min elapsed)", r, elapsed() / 60)
}
say("bootstrap replicates completed: %d of %d", length(boot_cls), NB)
if (length(boot_cls)) {
  bc <- rbindlist(boot_cls, idcol = "rep")
  stab <- merge(H[, .(from, to, class_range, class_idr)],
                bc[, .(stab_range = mean(class_range == class_range[1]),
                       stab_idr   = mean(class_idr == class_idr[1]), n_rep = .N), by = .(from, to)],
                by = c("from", "to"))
  # stability = share of replicates agreeing with the FULL DATA class
  stab2 <- merge(H[, .(from, to, class_range, class_idr)], bc, by = c("from", "to"), suffixes = c("", "_b"))
  stab <- stab2[, .(n_rep = .N,
                    stab_range = mean(class_range_b == class_range),
                    stab_idr   = mean(class_idr_b   == class_idr)), by = .(from, to)]
  H <- merge(H, stab, by = c("from", "to"), all.x = TRUE)
}
setorder(H, kind, -pearson)
fwrite(H, file.path(RES, "table_harmonisation_theory.csv"))
say("\n--- harmonisation restructured ---")
print(H[, .(pair = paste(from, "->", to), kind = substr(kind, 1, 14), r = round(pearson, 2),
            rel_range = round(rel_resid_range, 3), rel_idr = round(rel_resid_idr, 3),
            within = round(within_context_sd, 4), pooled = round(resid_sd, 4),
            tshift = round(target_shift, 3), sshift = round(survey_shift, 3),
            bound = round(pool_bound, 4), ok = pooled_exceeds_bound,
            cls_range = class_range, cls_idr = class_idr,
            stab = if ("stab_idr" %in% names(H)) round(stab_idr, 2) else NA)])

# ---------------------------------------------------------------------
# D. Proposition 3: directional asymmetry
# ---------------------------------------------------------------------
occ <- function(x, nominal) c(sd_over_range = sd(x) / diff(range(x)),
                              idr_over_nominal = idr(x) / nominal)
asym <- rbindlist(lapply(list(c("AUC","TSS"), c("AUC","Spearman"), c("TSS","Spearman"), c("AUC","Brier")), function(p) {
  a <- p[1]; b <- p[2]
  fwd <- H[from == a & to == b]; rev <- H[from == b & to == a]
  if (!nrow(fwd) || !nrow(rev)) return(NULL)
  sub <- d[is.finite(get(a)) & is.finite(get(b))]
  va <- var(sub[[a]]) / diff(range(sub[[a]]))^2; vb <- var(sub[[b]]) / diff(range(sub[[b]]))^2
  pred_ratio <- (va / vb) * (1 - rev$eta2) / (1 - fwd$eta2)   # r^2(b->a) / r^2(a->b)
  obs_ratio  <- (rev$rel_resid_range / fwd$rel_resid_range)^2
  data.table(forward = paste(a, "->", b), reverse = paste(b, "->", a),
             var_over_range2_from = va, var_over_range2_to = vb,
             eta2_fwd = fwd$eta2, eta2_rev = rev$eta2,
             predicted_ratio_sq = pred_ratio, observed_ratio_sq = obs_ratio)
}))
fwrite(asym, file.path(RES, "table_directional_asymmetry.csv"))
say("\n--- Proposition 3: directional asymmetry ---")
print(asym[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
say("range occupancy: AUC sd/range %.3f | TSS %.3f | Spearman %.3f",
    sd(d$AUC)/diff(range(d$AUC)), sd(d$TSS, na.rm=TRUE)/diff(range(d$TSS, na.rm=TRUE)),
    sd(d$Spearman, na.rm=TRUE)/diff(range(d$Spearman, na.rm=TRUE)))

# ---------------------------------------------------------------------
# E. Theorem 2(iii): calibration explains the AUC -> Brier residual
# ---------------------------------------------------------------------
sub <- d[is.finite(calib_slope) & is.finite(calib_intercept)]
rng <- quantile(sub$AUC, c(.02, .98), names = FALSE); xg <- seq(rng[1], rng[2], length.out = grid_n)
f <- fit_gam(sub, "AUC", "Brier", xg)
sub[, res_ab := f$resid]
sub[, log_slope := log(pmax(calib_slope, 1e-3))]
m_cal  <- lm(res_ab ~ log_slope + I(log_slope^2) + test_prev, data = sub)
m_cal2 <- lm(res_ab ~ log_slope + I(log_slope^2) + calib_intercept + I(calib_intercept^2) + test_prev, data = sub)
brt_share <- sub[, .(share_fits = mean(method == "BRT"),
                     share_rss  = sum(res_ab[method == "BRT"]^2) / sum(res_ab^2))]
cal_tab <- data.table(quantity = c("R2 of AUC->Brier residual on log calibration slope (quadratic) and prevalence",
                                   "R2 adding calibration intercept (quadratic)",
                                   "BRT share of fits", "BRT share of residual sum of squares"),
                      value = c(summary(m_cal)$r.squared, summary(m_cal2)$r.squared,
                                brt_share$share_fits, brt_share$share_rss))
fwrite(cal_tab, file.path(RES, "table_calibration_explains_residual.csv"))
say("\n--- Theorem 2(iii) ---"); print(cal_tab[, .(quantity, value = round(value, 3))])

# ---------------------------------------------------------------------
# F. Proposition 1: MAE - 2 Brier as calibration diagnostic
# ---------------------------------------------------------------------
d[, mae_2brier := MAE - 2 * Brier]
d[, abs_log_slope := abs(log(pmax(calib_slope, 1e-3)))]
f_tab <- d[is.finite(calib_slope), .(median_mae_minus_2brier = median(mae_2brier),
                                      median_calib_slope = median(calib_slope),
                                      spearman_with_abs_log_slope = cor(mae_2brier, abs_log_slope, method = "spearman"),
                                      n = .N), by = method][order(median_mae_minus_2brier)]
fwrite(f_tab, file.path(RES, "table_mae_minus_2brier_by_method.csv"))
say("\n--- Proposition 1 diagnostic: MAE - 2 Brier by method ---")
print(f_tab[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
say("overall Spearman(MAE - 2 Brier, |log calibration slope|) = %.3f",
    d[is.finite(calib_slope), cor(mae_2brier, abs_log_slope, method = "spearman")])

# ---------------------------------------------------------------------
# G. Corollary 5: Brier optimism against log n (matched draws)
# ---------------------------------------------------------------------
key <- c("dataset", "species", "n", "coverage", "rep", "method", "spatial")
w <- dcast(d, dataset + species + n + coverage + rep + method + spatial ~ target, value.var = "Brier")
w <- w[is.finite(conditional) & is.finite(marginal)]
w[, opt_brier := marginal - conditional]       # sign reversed: larger is worse
w[, log10n := log10(n)]
m_g <- lm(opt_brier ~ log10n * spatial, data = w)
co <- summary(m_g)$coefficients
g_tab <- data.table(term = rownames(co), estimate = co[, 1], se = co[, 2])
slope_ns <- co["log10n", 1]; slope_sp <- co["log10n", 1] + co["log10n:spatialTRUE", 1]
say("\n--- Corollary 5: slope of Brier optimism on log10 n ---")
say("non spatial %.4f (se %.4f) | spatial %.4f | interaction %.4f (se %.4f, p = %.3g)",
    slope_ns, co["log10n", 2], slope_sp, co["log10n:spatialTRUE", 1], co["log10n:spatialTRUE", 2], co["log10n:spatialTRUE", 4])
by_n <- w[, .(mean_opt_brier = mean(opt_brier), n_pairs = .N), by = .(spatial, n)][order(spatial, n)]
fwrite(g_tab, file.path(RES, "table_brier_optimism_slope.csv"))
fwrite(by_n, file.path(RES, "table_brier_optimism_by_n.csv"))
print(dcast(by_n, n ~ spatial, value.var = "mean_opt_brier")[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

# ---------------------------------------------------------------------
# H. Lemma 5: marginal spatial advantage against covariate count
# ---------------------------------------------------------------------
n_env <- c(EBS = 3L, `NS-IBTS` = 1L, GMEX = 3L, NEUS = 3L)
adv <- d[target == "marginal", .(AUC = mean(AUC)), by = .(dataset, species, spatial)]
adv <- dcast(adv, dataset + species ~ spatial, value.var = "AUC")
setnames(adv, c("FALSE", "TRUE"), c("non_spatial", "spatial"))
adv[, advantage := spatial - non_spatial]; adv[, n_env := n_env[dataset]]
h_tab <- data.table(quantity = c("Spearman(marginal spatial advantage, n covariates) over survey x species",
                                 "n survey x species combinations",
                                 "mean advantage, depth only survey (NS-IBTS)",
                                 "mean advantage, three covariate surveys"),
                    value = c(cor(adv$advantage, adv$n_env, method = "spearman"), nrow(adv),
                              adv[dataset == "NS-IBTS", mean(advantage)], adv[dataset != "NS-IBTS", mean(advantage)]))
fwrite(adv, file.path(RES, "table_spatial_advantage_by_species.csv"))
fwrite(h_tab, file.path(RES, "table_lemma5_covariate_count.csv"))
say("\n--- Lemma 5 ---"); print(h_tab[, .(quantity, value = round(value, 3))])

# ---------------------------------------------------------------------
# I. fit count arithmetic
# ---------------------------------------------------------------------
sf <- fread(file.path(RES, "split_failures.csv"))
n_planned <- length(DATASETS) * N_SPECIES_PER_DATASET * 30 * N_REPLICATES * length(METHODS)
say("\n--- fit count: %d x %d x 30 x %d x %d = %d; degenerate splits = %d; minus %d x %d = %d; rows in results = %d",
    length(DATASETS), N_SPECIES_PER_DATASET, N_REPLICATES, length(METHODS), n_planned,
    nrow(sf), nrow(sf), length(METHODS), n_planned - nrow(sf) * length(METHODS), nrow(d0))
say("theory analysis complete")
