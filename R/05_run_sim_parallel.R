# =====================================================================
# 05_run_sim_parallel.R
# Parallel, resumable, local multi core simulation driver. Each design
# cell (dataset x species x design cell) is an independent task,
# checkpointed atomically to results/parts_<scale>_r<reps>/. A PSOCK
# cluster is used, so the same script runs on Windows, macOS and Linux.
# Re running resumes from wherever it stopped.
#   SMD_SCALE=full SMD_CORES=8 Rscript R/05_run_sim_parallel.R
# Overrides: SMD_CORES, SMD_REPLICATES, SMD_DATASETS (e.g. "EBS,GMEX"),
#            SMD_MAX_CELLS (cap for a trial run)
# =====================================================================
t_start <- Sys.time()
suppressMessages({ library(parallel); library(data.table) })
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
PROJECT_ROOT <- Sys.getenv("SMD_ROOT"); setwd(PROJECT_ROOT)
for (f in c("config.R", "02_subsample.R", "03_methods.R", "04_metrics.R", "05_helpers.R"))
  source(file.path("R", f))
log_file <- file.path(PATHS$logs, "run_sim_parallel.log")
logline <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))
  cat(m, "\n"); cat(m, "\n", file = log_file, append = TRUE) }
parts_dir <- file.path(PATHS$results, sprintf("parts_%s_r%d", SCALE, N_REPLICATES))
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)
sel <- Sys.getenv("SMD_DATASETS", unset = "")
ds_codes <- if (nzchar(sel)) trimws(strsplit(sel, ",")[[1]]) else names(DATASETS)

species_per <- integer(0); tasks <- list(); k <- 0L
for (code in ds_codes) {
  di <- match(code, names(DATASETS))
  pf <- file.path(PATHS$processed, paste0(code, "_processed.rds"))
  if (!file.exists(pf)) { logline("MISSING processed %s; run 01_preprocess.R", code); next }
  nsp <- nrow(readRDS(pf)$species); species_per[code] <- nsp
  for (si in seq_len(nsp)) for (ci in seq_len(n_cells)) {
    cf <- file.path(parts_dir, sprintf("cell_%s_%d_%d.rds", code, si, ci))
    if (!file.exists(cf)) { k <- k + 1L; tasks[[k]] <- list(code = code, di = di, si = si, ci = ci) }
  }
}
n_cells_total <- sum(species_per) * n_cells
fits_per_cell <- N_REPLICATES * length(METHODS)
logline("scale=%s | datasets=%s | cores=%d | reps=%d | methods=%d", SCALE,
        paste(ds_codes, collapse = ","), N_CORES, N_REPLICATES, length(METHODS))
logline("cells total=%d | done=%d | pending=%d | fits at completion=%s",
        n_cells_total, n_cells_total - length(tasks), length(tasks),
        format(n_cells_total * fits_per_cell, big.mark = ","))
maxc <- Sys.getenv("SMD_MAX_CELLS", unset = "")
if (nzchar(maxc) && length(tasks) > as.integer(maxc)) {
  tasks <- tasks[seq_len(as.integer(maxc))]; logline("SMD_MAX_CELLS: processing %d cells", length(tasks)) }

process_task <- function(task) {
  cf <- file.path(parts_dir, sprintf("cell_%s_%d_%d.rds", task$code, task$si, task$ci))
  if (file.exists(cf)) return(list(n_fits = 0L, n_fail = 0L, error = NA_character_))
  # an R error inside a cell is reported back rather than left to kill the
  # worker; the cell stays pending and is retried on the next run
  out <- tryCatch({
    r <- .smd_get_data(task$code)
    o <- run_cell(task$code, task$di, task$si, task$ci, r)
    tmp <- paste0(cf, ".tmp-", Sys.getpid())
    saveRDS(list(fits = o$fits, fails = o$fails), tmp, compress = TRUE)
    file.rename(tmp, cf)
    list(n_fits = o$n_fits, n_fail = o$n_fail, error = NA_character_)
  }, error = function(e) list(n_fits = 0L, n_fail = 0L,
                              error = sprintf("cell %s/%d/%d: %s", task$code, task$si, task$ci, conditionMessage(e))))
  gc(verbose = FALSE)
  out
}

# ---- cluster set up, reusable so a dead cluster can be rebuilt ---------
worker_env <- c(SMD_ROOT = PROJECT_ROOT, SMD_SCALE = SCALE, SMD_REPLICATES = as.character(N_REPLICATES))
setup_cluster <- function(ncores) {
  cl <- makeCluster(ncores, outfile = file.path(PATHS$logs, "workers.log"))
  clusterExport(cl, c("worker_env", "parts_dir"), envir = environment(setup_cluster))
  clusterEvalQ(cl, {
    suppressMessages(library(data.table)); do.call(Sys.setenv, as.list(worker_env))
    setwd(Sys.getenv("SMD_ROOT"))
    for (f in c("config.R", "02_subsample.R", "03_methods.R", "04_metrics.R", "05_helpers.R"))
      source(file.path("R", f))
    .smd_get_data <- local({ cache <- list(); function(code) {
      if (is.null(cache[[code]])) cache[[code]] <<- readRDS(file.path(PATHS$processed, paste0(code, "_processed.rds")))
      cache[[code]] } })
    invisible(TRUE) })
  clusterExport(cl, "process_task", envir = environment(setup_cluster))
  cl
}
environment(setup_cluster) <- environment()

if (length(tasks)) {
  ncores <- N_CORES
  cl <- setup_cluster(ncores); on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  n_total <- length(tasks)
  # smaller chunks so that a lost cluster costs at most a few minutes of work
  chunk_size <- max(ncores, min(24L, ceiling(n_total / 50)))
  chunks <- split(tasks, ceiling(seq_along(tasks) / chunk_size))
  done <- 0L; fits_run <- 0L; fails_run <- 0L; t_run <- Sys.time(); crashes <- 0L
  kk <- 1L
  while (kk <= length(chunks)) {
    res <- tryCatch(parLapplyLB(cl, chunks[[kk]], process_task), error = function(e) e)
    if (inherits(res, "error")) {
      # a worker process died (sleep, out of memory, or a crash in a library):
      # rebuild the cluster and retry the chunk; completed cells are on disk
      crashes <- crashes + 1L
      logline("WORKER LOSS %d at chunk %d: %s", crashes, kk, conditionMessage(res))
      try(stopCluster(cl), silent = TRUE); Sys.sleep(10)
      if (crashes %% 3L == 0L && ncores > 1L) {
        ncores <- max(1L, ncores %/% 2L)
        logline("repeated losses; continuing with %d workers (memory is the usual cause)", ncores)
      }
      if (crashes > 12L) stop("too many worker losses; see logs/workers.log and the notes in README_full_scale.md")
      cl <- setup_cluster(ncores)
      next
    }
    errs <- Filter(Negate(is.na), vapply(res, function(x) x$error, ""))
    for (e in errs) logline("CELL ERROR (left pending): %s", e)
    done <- done + length(chunks[[kk]])
    fits_run <- fits_run + sum(vapply(res, function(x) x$n_fits, 0L))
    fails_run <- fails_run + sum(vapply(res, function(x) x$n_fail, 0L))
    el <- as.numeric(Sys.time() - t_run, units = "mins"); rate <- done / max(el, 1e-6)
    logline("chunk %d/%d | cells %d/%d (%.1f%%) | elapsed %.1f min | ETA %.1f min | fits %d",
            kk, length(chunks), done, n_total, 100 * done / n_total, el, (n_total - done) / rate, fits_run)
    kk <- kk + 1L
  }
  stopCluster(cl)
  logline("worker phase done | %d cells | %d fits | %d split failures | %d worker losses", done, fits_run, fails_run, crashes)
} else logline("nothing pending; all cells complete")

part_files <- list.files(parts_dir, pattern = "^cell_.*\\.rds$", full.names = TRUE)
if (length(part_files)) {
  parts <- lapply(part_files, readRDS)
  fit_metrics <- rbindlist(lapply(parts, `[[`, "fits"), fill = TRUE)
  split_failures <- rbindlist(lapply(parts, `[[`, "fails"), fill = TRUE)
  fwrite(fit_metrics, file.path(PATHS$results, "fit_metrics.csv"))
  saveRDS(fit_metrics, file.path(PATHS$results, "fit_metrics.rds"))
  if (nrow(split_failures)) fwrite(split_failures, file.path(PATHS$results, "split_failures.csv"))
  logline("MERGE | cells %d/%d | fits %d | failures %d | %.1f min total", length(part_files),
          n_cells_total, nrow(fit_metrics), nrow(split_failures), as.numeric(Sys.time() - t_start, units = "mins"))
  if (length(part_files) == n_cells_total)
    logline("ALL DONE | converged %.4f | mean fit %.3fs", mean(fit_metrics$converged, na.rm = TRUE), mean(fit_metrics$fit_seconds, na.rm = TRUE))
  else logline("PARTIAL | %d of %d cells; re run to continue", length(part_files), n_cells_total)
}
