# =====================================================================
# 05_run_sim_shard.R
# Serial, socket free driver for one shard of the full scale run.
# Processes pending replicate blocks whose index satisfies
#   index %% SMD_SHARDS == SMD_SHARD
# so N shards launched as separate processes use N cores without a
# cluster. Each block is checkpointed as the parallel driver does
# (results/parts_<scale>_r<reps>/cell_<code>_<si>_<ci>_b<k>.rds), so the
# two drivers are interchangeable and resumable. Stops taking new blocks
# once SMD_TIME_BUDGET seconds have elapsed (default: unlimited).
#   SMD_SCALE=full SMD_SHARDS=2 SMD_SHARD=0 Rscript R/05_run_sim_shard.R
# =====================================================================
t_start <- Sys.time()
suppressMessages(library(data.table))
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT"))
for (f in c("config.R", "02_subsample.R", "03_methods.R", "04_metrics.R", "05_helpers.R"))
  source(file.path("R", f))
SHARDS <- as.integer(Sys.getenv("SMD_SHARDS", "1")); SHARD <- as.integer(Sys.getenv("SMD_SHARD", "0"))
BUDGET <- as.numeric(Sys.getenv("SMD_TIME_BUDGET", "Inf"))
REP_BLOCK <- as.integer(Sys.getenv("SMD_REP_BLOCK", "25"))
rep_blocks <- split(seq_len(N_REPLICATES), ceiling(seq_len(N_REPLICATES) / REP_BLOCK))
parts_dir <- file.path(PATHS$results, sprintf("parts_%s_r%d", SCALE, N_REPLICATES))
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(PATHS$logs, sprintf("shard_%d.log", SHARD))
logline <- function(...) { m <- sprintf("[%s] shard %d: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), SHARD, sprintf(...))
  cat(m, "\n"); cat(m, "\n", file = log_file, append = TRUE) }
sel <- Sys.getenv("SMD_DATASETS", unset = "")
ds_codes <- if (nzchar(sel)) trimws(strsplit(sel, ",")[[1]]) else names(DATASETS)

done <- list.files(parts_dir)
tasks <- list(); k <- 0L
for (code in ds_codes) {
  di <- match(code, names(DATASETS))
  nsp <- nrow(readRDS(file.path(PATHS$processed, paste0(code, "_processed.rds")))$species)
  for (si in seq_len(nsp)) for (ci in seq_len(n_cells)) {
    if (sprintf("cell_%s_%d_%d.rds", code, si, ci) %in% done) next
    for (b in seq_along(rep_blocks)) {
      k <- k + 1L                                   # global index, so shards partition consistently
      if (k %% SHARDS != SHARD) next
      if (sprintf("cell_%s_%d_%d_b%d.rds", code, si, ci, b) %in% done) next
      tasks[[length(tasks) + 1L]] <- list(code = code, di = di, si = si, ci = ci, b = b)
    }
  }
}
logline("pending blocks in this shard: %d (of %d total units)", length(tasks), k)
cache <- list(); n_done <- 0L; n_fits <- 0L
for (task in tasks) {
  if (as.numeric(Sys.time() - t_start, units = "secs") > BUDGET) { logline("time budget reached"); break }
  cf <- file.path(parts_dir, sprintf("cell_%s_%d_%d_b%d.rds", task$code, task$si, task$ci, task$b))
  if (file.exists(cf)) next
  if (is.null(cache[[task$code]]))
    cache[[task$code]] <- readRDS(file.path(PATHS$processed, paste0(task$code, "_processed.rds")))
  t0 <- Sys.time()
  out <- tryCatch(run_cell(task$code, task$di, task$si, task$ci, cache[[task$code]], reps = rep_blocks[[task$b]]),
                  error = function(e) { logline("CELL ERROR %s/%d/%d b%d: %s", task$code, task$si, task$ci, task$b, conditionMessage(e)); NULL })
  if (is.null(out)) next
  tmp <- paste0(cf, ".tmp-", Sys.getpid())
  saveRDS(list(fits = out$fits, fails = out$fails), tmp, compress = TRUE); file.rename(tmp, cf)
  n_done <- n_done + 1L; n_fits <- n_fits + out$n_fits
  logline("%s sp%d cell%d b%d | %d fits %d fails | %.1f s", task$code, task$si, task$ci, task$b,
          out$n_fits, out$n_fail, as.numeric(Sys.time() - t0, units = "secs"))
}
logline("finished: %d blocks, %d fits, %.1f min", n_done, n_fits, as.numeric(Sys.time() - t_start, units = "mins"))
