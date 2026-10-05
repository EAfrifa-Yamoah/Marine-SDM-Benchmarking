# 05c_merge_parts.R  Merge block checkpoints into results/fit_metrics.{rds,csv}
suppressMessages(library(data.table))
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT")); source(file.path("R", "config.R"))
parts_dir <- file.path(PATHS$results, sprintf("parts_%s_r%d", SCALE, N_REPLICATES))
fs <- list.files(parts_dir, pattern = "^cell_.*\\.rds$", full.names = TRUE)
t0 <- Sys.time(); cat(sprintf("merging %d checkpoint files from %s\n", length(fs), parts_dir))
parts <- lapply(fs, readRDS)
fits  <- rbindlist(lapply(parts, `[[`, "fits"), fill = TRUE)
fails <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, lapply(parts, `[[`, "fails"))
fails <- if (length(fails)) rbindlist(fails, fill = TRUE) else data.table()
setorder(fits, dataset, species, n, coverage, target, rep, method)
saveRDS(fits, file.path(PATHS$results, "fit_metrics.rds"))
fwrite(fits,  file.path(PATHS$results, "fit_metrics.csv"))
if (nrow(fails)) fwrite(fails, file.path(PATHS$results, "split_failures.csv")) else
  unlink(file.path(PATHS$results, "split_failures.csv"))
cat(sprintf("fits %s | split failures %d | converged %.4f | mean fit %.3fs | %.1f min\n",
  format(nrow(fits), big.mark = ","), nrow(fails), mean(fits$converged, na.rm = TRUE),
  mean(fits$fit_seconds, na.rm = TRUE), as.numeric(Sys.time() - t0, units = "mins")))
print(fits[, .N, by = .(dataset, method)][order(dataset, method)])
