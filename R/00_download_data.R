# 00_download_data.R  Download the four cleaned FISHGLOB surveys to data_raw/.
# Skipped for files already present unless SMD_FORCE_DOWNLOAD is set.
# Source: Maureaud et al. (2024), https://github.com/fishglob/FishGlob_data
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
raw_dir <- file.path(Sys.getenv("SMD_ROOT"), "data_raw"); dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
BASE <- "https://raw.githubusercontent.com/fishglob/FishGlob_data/main/outputs/Cleaned_data"
CODES <- c("EBS", "NS-IBTS", "GMEX", "NEUS"); force <- nzchar(Sys.getenv("SMD_FORCE_DOWNLOAD"))
options(timeout = max(600, getOption("timeout")))
for (code in CODES) {
  dest <- file.path(raw_dir, paste0(code, "_clean.RData"))
  if (file.exists(dest) && !force) { cat("[download]", code, "present, skipping\n"); next }
  cat("[download]", code, "\n"); utils::download.file(sprintf("%s/%s_clean.RData", BASE, code), dest, mode = "wb", quiet = TRUE)
}
stopifnot(all(file.exists(file.path(raw_dir, paste0(CODES, "_clean.RData")))))
cat("[download] all four surveys available\n")
