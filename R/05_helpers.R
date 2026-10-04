# =====================================================================
# 05_helpers.R
# Shared cell logic for the serial and parallel drivers. Assumes config.R,
# 02_subsample.R, 03_methods.R and 04_metrics.R are sourced.
#
# Beyond the pilot's metrics, every fit now records:
#   REF, CAL      refinement and calibration components of the Brier score,
#                 by isotonic recalibration of the test predictions
#                 (Theorem 1(d), Theorem 2(ii))
#   interp_*      under the marginal target, the SAME fitted predictor
#                 evaluated on an interpolation test set inside the training
#                 window, giving the same fit optimism contrast (Theorem 4)
#   fill_km       fill distance of the training design (Proposition 5)
#   separation_km separation of the held out block from the training
#                 locations (Theorem 6)
#   cov_shift     energy distance between training and test covariates
#                 (Corollary 4)
# =====================================================================

design <- expand.grid(n = SAMPLE_SIZES, coverage = COVERAGE_LEVELS,
                      target = INFERENCE_TARGETS, stringsAsFactors = FALSE)
design <- design[order(design$target, design$coverage, design$n), ]
n_cells <- nrow(design)
reg <- method_registry; rownames(reg) <- reg$method

cell_seed <- function(di, si, ci, rep)
  (GLOBAL_SEED + di * 1000003L + si * 10007L + ci * 101L + rep) %% .Machine$integer.max

# refinement and calibration by isotonic regression of y on p (pool
# adjacent violators through stats::isoreg); REF = mean kappa(1 - kappa),
# CAL = BS - REF, so that BS = CAL + REF holds exactly on the test set
ref_cal <- function(y, p) {
  ok <- is.finite(p) & is.finite(y)
  if (sum(ok) < 10 || length(unique(y[ok])) < 2) return(c(REF = NA_real_, CAL = NA_real_))
  o <- order(p[ok]); ir <- stats::isoreg(p[ok][o], y[ok][o])
  kap <- ir$yf
  ref <- mean(kap * (1 - kap)); bs <- mean((y[ok] - p[ok])^2)
  c(REF = ref, CAL = bs - ref)
}

# energy distance between two covariate samples (Szekely and Rizzo), on
# standardised covariates; a base R implementation
energy_distance <- function(A, B, max_n = 300L) {
  A <- as.matrix(A); B <- as.matrix(B)
  if (nrow(A) > max_n) A <- A[sample.int(nrow(A), max_n), , drop = FALSE]
  if (nrow(B) > max_n) B <- B[sample.int(nrow(B), max_n), , drop = FALSE]
  Z <- rbind(A, B); mu <- colMeans(Z); sdv <- apply(Z, 2, sd); sdv[sdv == 0] <- 1
  A <- sweep(sweep(A, 2, mu), 2, sdv, "/"); B <- sweep(sweep(B, 2, mu), 2, sdv, "/")
  md <- function(X, Y) mean(sqrt(pmax(outer(rowSums(X^2), rowSums(Y^2), "+") - 2 * X %*% t(Y), 0)))
  2 * md(A, B) - md(A, A) - md(B, B)
}

train_metrics <- function(y, p) {
  m <- compute_metrics(y, p)
  c(train_AUC = unname(m["AUC"]), train_R2 = unname(m["R2"]), train_Brier = unname(m["Brier"]))
}
na_te <- setNames(rep(NA_real_, length(METRIC_NAMES)), METRIC_NAMES)
na_tr <- c(train_AUC = NA_real_, train_R2 = NA_real_, train_Brier = NA_real_)
na_in <- c(interp_AUC = NA_real_, interp_Brier = NA_real_, interp_TSS = NA_real_,
           interp_R2 = NA_real_, interp_REF = NA_real_, interp_CAL = NA_real_,
           interp_prev = NA_real_, n_interp = NA_integer_)

make_fit_row <- function(code, region, focal, srow, nn, cov, tgt, rep, edge, mth,
                         ntr, nte, tr_prev, te_prev, te_m, tr_m, rc, inm, geo, conv, dt) {
  data.frame(
    dataset = code, region = region, species = focal, band = srow$band,
    prevalence = srow$prevalence, niche = srow$niche, trophic = srow$trophic,
    n = nn, coverage = cov, coverage_frac = COVERAGE_FRAC[[cov]],
    target = tgt, rep = rep, edge = edge,
    method = mth, family = reg[mth, "family"],
    spatial = reg[mth, "spatial"], multisp = reg[mth, "multisp"],
    n_train = ntr, n_test = nte, train_prev = tr_prev, test_prev = te_prev,
    AUC = unname(te_m["AUC"]), TSS = unname(te_m["TSS"]),
    RMSE = unname(te_m["RMSE"]), MAE = unname(te_m["MAE"]),
    R2 = unname(te_m["R2"]), Spearman = unname(te_m["Spearman"]),
    Brier = unname(te_m["Brier"]),
    calib_intercept = unname(te_m["calib_intercept"]),
    calib_slope = unname(te_m["calib_slope"]),
    REF = unname(rc["REF"]), CAL = unname(rc["CAL"]),
    train_AUC = unname(tr_m["train_AUC"]), train_R2 = unname(tr_m["train_R2"]),
    train_Brier = unname(tr_m["train_Brier"]),
    gap_AUC = unname(tr_m["train_AUC"]) - unname(te_m["AUC"]),
    interp_AUC = unname(inm["interp_AUC"]), interp_Brier = unname(inm["interp_Brier"]),
    interp_TSS = unname(inm["interp_TSS"]), interp_R2 = unname(inm["interp_R2"]),
    interp_REF = unname(inm["interp_REF"]), interp_CAL = unname(inm["interp_CAL"]),
    interp_prev = unname(inm["interp_prev"]), n_interp = unname(inm["n_interp"]),
    fill_km = geo$fill_km, separation_km = geo$separation_km,
    cov_shift = geo$cov_shift,
    converged = conv, fit_seconds = dt, stringsAsFactors = FALSE)
}

# process one design cell across all replicates and methods
run_cell <- function(code, di, si, ci, r, reps = seq_len(N_REPLICATES)) {
  hauls <- r$hauls; occ <- r$occ; env_used <- r$env_used
  sp_all <- r$species$accepted_name; spinfo <- r$species
  region <- unname(DATASETS[code]); focal <- sp_all[si]
  yf <- occ[[focal]][match(hauls$haul_id, occ$haul_id)]
  srow <- spinfo[spinfo$accepted_name == focal, ]
  occ_lookup <- function(df) as.data.frame(occ[match(df$haul_id, occ$haul_id), ..sp_all])
  nn <- design$n[ci]; cov <- design$coverage[ci]; tgt <- design$target[ci]
  fit_rows <- list(); fail_rows <- list(); fi <- 0L; qi <- 0L

  for (rep in reps) {
    set.seed(cell_seed(di, si, ci, rep))
    sp <- make_split(hauls, yf, env_used, n = nn, coverage = cov,
                     target = tgt, rep_id = rep, test_size = TEST_SIZE)
    if (!isTRUE(sp$ok)) {
      qi <- qi + 1L
      fail_rows[[qi]] <- data.frame(dataset = code, region = region, species = focal,
        band = srow$band, n = nn, coverage = cov, target = tgt, rep = rep,
        reason = sp$reason, stringsAsFactors = FALSE)
      next
    }
    co_tr <- occ_lookup(sp$train); co_te <- occ_lookup(sp$test)
    ntr <- nrow(sp$train); nte <- nrow(sp$test)
    tr_prev <- mean(sp$train$y); te_prev <- mean(sp$test$y)
    geo <- list(fill_km = sp$fill_km, separation_km = sp$separation_km,
                cov_shift = tryCatch(energy_distance(sp$train[, env_used, drop = FALSE],
                                                     sp$test[, env_used, drop = FALSE]),
                                     error = function(e) NA_real_))
    has_interp <- !is.null(sp$interp) && nrow(sp$interp) >= 10 && length(unique(sp$interp$y)) == 2

    for (mth in METHODS) {
      t0 <- Sys.time()
      # fit once; predict on test, on train, and (marginal) on the interpolation set
      res <- tryCatch({
        if (mth == "JointSDM")
          METHOD_FUNS[[mth]](sp$train, sp$test, env_used, co_occ_train = co_tr,
                             co_occ_test = co_te, focal = focal)
        else METHOD_FUNS[[mth]](sp$train, sp$test, env_used)
      }, error = function(e) list(prob = NULL, err = conditionMessage(e)))
      inm <- na_in
      if (has_interp && is.null(res$err) && !is.null(res$prob) && !is.null(res$model)) {
        pin <- tryCatch(res$predict(res$model, sp$interp), error = function(e) NULL)
        if (!is.null(pin)) {
          m_in <- compute_metrics(sp$interp$y, pin); rc_in <- ref_cal(sp$interp$y, pin)
          inm <- c(interp_AUC = unname(m_in["AUC"]), interp_Brier = unname(m_in["Brier"]),
                   interp_TSS = unname(m_in["TSS"]), interp_R2 = unname(m_in["R2"]),
                   interp_REF = unname(rc_in["REF"]), interp_CAL = unname(rc_in["CAL"]),
                   interp_prev = mean(sp$interp$y), n_interp = nrow(sp$interp))
        }
      }
      dt <- as.numeric(Sys.time() - t0, units = "secs")
      if (!is.null(res$err) || is.null(res$prob)) {
        te_m <- na_te; tr_m <- na_tr; rc <- c(REF = NA_real_, CAL = NA_real_); conv <- FALSE
      } else {
        te_m <- compute_metrics(sp$test$y, res$prob)
        tr_m <- if (!is.null(res$prob_train)) train_metrics(sp$train$y, res$prob_train) else na_tr
        rc <- ref_cal(sp$test$y, res$prob); conv <- isTRUE(res$converged)
      }
      fi <- fi + 1L
      fit_rows[[fi]] <- make_fit_row(code, region, focal, srow, nn, cov, tgt, rep, sp$edge,
        mth, ntr, nte, tr_prev, te_prev, te_m, tr_m, rc, inm, geo, conv, dt)
    }
  }
  list(fits  = if (length(fit_rows))  data.table::rbindlist(fit_rows)  else data.table::data.table(),
       fails = if (length(fail_rows)) data.table::rbindlist(fail_rows) else data.table::data.table(),
       n_fits = fi, n_fail = qi)
}
