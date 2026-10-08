# run_full_local.R   One command, full study, local multi core execution.
#   Rscript run_full_local.R
# Options: SMD_CORES, SMD_REPLICATES (default 100), SMD_DATASETS, SMD_FORCE_PREPROCESS=1
Sys.setenv(SMD_SCALE = Sys.getenv("SMD_SCALE", unset = "full"))
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
ROOT <- Sys.getenv("SMD_ROOT"); setwd(ROOT); dir.create("logs", showWarnings = FALSE)
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = "logs/run_full_local.log", append = TRUE) }
stage <- function(script, label) { say("==== %s ====", label); t0 <- Sys.time()
  st <- system2(RSCRIPT, shQuote(file.path("R", script)))
  if (st != 0L) stop(sprintf("stage '%s' failed (status %d); see logs/", label, st))
  say("%s complete in %.1f min", label, as.numeric(Sys.time() - t0, units = "mins")) }
codes <- c("EBS", "NS-IBTS", "GMEX", "NEUS"); t_all <- Sys.time()
say("FULL LOCAL RUN | cores=%s | replicates=%s", Sys.getenv("SMD_CORES", "auto"), Sys.getenv("SMD_REPLICATES", "100"))
if (!all(file.exists(file.path("data_raw", paste0(codes, "_clean.RData"))))) stage("00_download_data.R", "download raw data")
if (!all(file.exists(file.path("data_processed", paste0(codes, "_processed.rds")))) || nzchar(Sys.getenv("SMD_FORCE_PREPROCESS")))
  stage("01_preprocess.R", "preprocess surveys") else say("processed frames present, skipping preprocessing")
stage("05_run_sim_parallel.R", "parallel simulation (resumable)")
stage("05b_oracle.R",          "near oracle intrinsic term (resumable)")
stage("05d_practical_range.R", "residual variogram range of every species")
stage("06_analysis.R",         "hierarchical analysis")
stage("06b_theory_full.R",     "theory driven analyses")
stage("07_harmonise.R",        "metric harmonisation and package")
stage("08_figures.R",          "result figures")
stage("09_explore_spatial.R",  "spatial exploration figures")
say("FULL LOCAL RUN complete in %.1f min | outputs in results/ figures/ smdMetricHarmonise/", as.numeric(Sys.time() - t_all, units = "mins"))
