# =====================================================================
# 05_run_sim.R
# Main simulation driver.
#
# For every dataset, focal species, design cell (sample size x coverage x
# inference target) and replicate, one train/test split is drawn and all
# seven methods are fitted on it. Nine metrics are computed on the honest
# test set, and a compact set on the training set, so the analysis can
# quantify overfitting (train minus test gap) and validation optimism
# (conditional minus marginal target). Split failures and method non
# convergence are recorded so the failure mode analysis sees genuine data
# limited breakdowns rather than silent gaps.
#
# The same code runs the pilot and the full study; only config.R changes.
# Progress is written to logs/run_sim.log for monitoring.
# =====================================================================

t_start <- Sys.time()
suppressMessages({ library(data.table) })
src <- function(f) source(file.path("R", f))
src("config.R"); src("02_subsample.R")
src("03_methods.R"); src("04_metrics.R")

log_file <- file.path(PATHS$logs, "run_sim.log")
logline <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE)
}
# wall clock budget per invocation; the driver stops launching new species
# once exceeded, having checkpointed everything finished so far. Repeated
# invocations resume and converge. Default is effectively unlimited.
TIME_BUDGET <- as.numeric(Sys.getenv("SMD_TIME_BUDGET", unset = "1e9"))
logline("---- invocation start (budget %.0fs) ----", TIME_BUDGET)

# deterministic per cell seed so any fit is independently reproducible
cell_seed <- function(di, si, ci, rep) {
  (GLOBAL_SEED + di * 1000003L + si * 10007L + ci * 101L + rep) %% .Machine$integer.max
}

# compact training metrics for the overfitting gap
train_metrics <- function(y, p) {
  m <- compute_metrics(y, p)
  c(train_AUC = unname(m["AUC"]), train_R2 = unname(m["R2"]),
    train_Brier = unname(m["Brier"]))
}

design <- expand.grid(n = SAMPLE_SIZES, coverage = COVERAGE_LEVELS,
                      target = INFERENCE_TARGETS,
                      stringsAsFactors = FALSE)
design <- design[order(design$target, design$coverage, design$n), ]
n_cells <- nrow(design)

reg <- method_registry
rownames(reg) <- reg$method

# optional dataset subset (comma separated codes) for chunked execution;
# per dataset outputs make the whole run resumable and cheap to restart
sel <- Sys.getenv("SMD_DATASETS", unset = "")
ds_codes <- if (nzchar(sel)) trimws(strsplit(sel, ",")[[1]]) else names(DATASETS)

parts_dir <- file.path(PATHS$results, "parts")
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)

logline("run start | scale=%s | datasets=%s | cells=%d | reps=%d | methods=%d",
        SCALE, paste(ds_codes, collapse = ","), n_cells,
        N_REPLICATES, length(METHODS))

# build one fit result row (kept as a helper so the cell loop stays legible)
make_fit_row <- function(code, region, focal, srow, nn, cov, tgt, rep, edge,
                         mth, ntr, nte, tr_prev, te_prev, te_m, tr_m, conv, dt) {
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
    train_AUC = unname(tr_m["train_AUC"]), train_R2 = unname(tr_m["train_R2"]),
    train_Brier = unname(tr_m["train_Brier"]),
    gap_AUC = unname(tr_m["train_AUC"]) - unname(te_m["AUC"]),
    converged = conv, fit_seconds = dt, stringsAsFactors = FALSE)
}

na_te <- setNames(rep(NA_real_, length(METRIC_NAMES)), METRIC_NAMES)
na_tr <- c(train_AUC = NA_real_, train_R2 = NA_real_, train_Brier = NA_real_)

for (code in ds_codes) {
  di <- match(code, names(DATASETS))
  pf <- file.path(PATHS$processed, paste0(code, "_processed.rds"))
  if (!file.exists(pf)) { logline("MISSING processed %s, skipping", code); next }
  r <- readRDS(pf)
  hauls <- r$hauls; occ <- r$occ; env_used <- r$env_used
  sp_all <- r$species$accepted_name
  spinfo <- r$species
  region <- unname(DATASETS[code])

  occ_lookup <- function(df) {
    idx <- match(df$haul_id, occ$haul_id)
    as.data.frame(occ[idx, ..sp_all])
  }

  for (si in seq_along(sp_all)) {
    focal <- sp_all[si]
    yf <- occ[[focal]][match(hauls$haul_id, occ$haul_id)]
    srow <- spinfo[spinfo$accepted_name == focal, ]

    for (ci in seq_len(n_cells)) {
      cell_file <- file.path(parts_dir,
                             sprintf("cell_%s_%d_%d.rds", code, si, ci))
      if (file.exists(cell_file)) next                  # cell already done
      if (as.numeric(Sys.time() - t_start, units = "secs") > TIME_BUDGET) {
        logline("time budget reached at %s sp%d cell%d, stopping invocation",
                code, si, ci)
        quit(save = "no", status = 0)
      }
      nn  <- design$n[ci]; cov <- design$coverage[ci]; tgt <- design$target[ci]
      fit_rows <- vector("list", 0L); fail_rows <- vector("list", 0L)
      fi <- 0L; qi <- 0L

      for (rep in seq_len(N_REPLICATES)) {
        set.seed(cell_seed(di, si, ci, rep))
        sp <- make_split(hauls, yf, env_used, n = nn, coverage = cov,
                         target = tgt, rep_id = rep, test_size = TEST_SIZE)
        if (!isTRUE(sp$ok)) {
          qi <- qi + 1L
          fail_rows[[qi]] <- data.frame(
            dataset = code, region = region, species = focal,
            band = srow$band, n = nn, coverage = cov, target = tgt,
            rep = rep, reason = sp$reason, stringsAsFactors = FALSE)
          next
        }
        co_tr <- occ_lookup(sp$train); co_te <- occ_lookup(sp$test)
        ntr <- nrow(sp$train); nte <- nrow(sp$test)
        tr_prev <- mean(sp$train$y); te_prev <- mean(sp$test$y)

        for (mth in METHODS) {
          t0 <- Sys.time()
          res <- tryCatch({
            if (mth == "JointSDM") {
              METHOD_FUNS[[mth]](sp$train, sp$test, env_used,
                                 co_occ_train = co_tr, co_occ_test = co_te,
                                 focal = focal)
            } else METHOD_FUNS[[mth]](sp$train, sp$test, env_used)
          }, error = function(e) list(prob = NULL, err = conditionMessage(e)))
          dt <- as.numeric(Sys.time() - t0, units = "secs")

          if (!is.null(res$err) || is.null(res$prob)) {
            te_m <- na_te; tr_m <- na_tr; conv <- FALSE
          } else {
            te_m <- compute_metrics(sp$test$y, res$prob)
            tr_m <- if (!is.null(res$prob_train))
              train_metrics(sp$train$y, res$prob_train) else na_tr
            conv <- isTRUE(res$converged)
          }
          fi <- fi + 1L
          fit_rows[[fi]] <- make_fit_row(code, region, focal, srow, nn, cov,
            tgt, rep, sp$edge, mth, ntr, nte, tr_prev, te_prev, te_m, tr_m,
            conv, dt)
        }
      }

      saveRDS(list(
        fits  = if (length(fit_rows))  rbindlist(fit_rows)  else data.table(),
        fails = if (length(fail_rows)) rbindlist(fail_rows) else data.table()),
        cell_file)
    }
    logline("%s sp%d (%s) cells checkpointed through %d",
            code, si, srow$band, n_cells)
  }
}

# ---- merge all completed cell parts ----------------------------------
part_files <- list.files(parts_dir, pattern = "^cell_.*\\.rds$",
                         full.names = TRUE)
n_expected <- length(ds_codes) * 3L * n_cells
if (length(part_files)) {
  parts <- lapply(part_files, readRDS)
  fit_metrics <- rbindlist(lapply(parts, `[[`, "fits"), fill = TRUE)
  split_failures <- rbindlist(lapply(parts, `[[`, "fails"), fill = TRUE)
  fwrite(fit_metrics, file.path(PATHS$results, "fit_metrics.csv"))
  saveRDS(fit_metrics, file.path(PATHS$results, "fit_metrics.rds"))
  if (nrow(split_failures))
    fwrite(split_failures, file.path(PATHS$results, "split_failures.csv"))

  el <- as.numeric(Sys.time() - t_start, units = "mins")
  logline("MERGE | cells done=%d/%d | fits=%d | split failures=%d | elapsed=%.1f min",
          length(part_files), n_expected, nrow(fit_metrics),
          nrow(split_failures), el)
  if (length(part_files) == n_expected)
    logline("ALL DONE | converged rate=%.3f | mean fit=%.3fs",
            mean(fit_metrics$converged, na.rm = TRUE),
            mean(fit_metrics$fit_seconds, na.rm = TRUE))
}

