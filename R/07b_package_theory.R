# =====================================================================
# 07b_package_theory.R
# Refresh the theory table inside the smdMetricHarmonise package data:
# the kind of every translation (identity or empirical), its class under
# both normalisation rules and the stability of that class under the
# species bootstrap. 07_harmonise.R reads results/table_harmonisation_theory.csv
# when it exists, but at full scale 07 runs before R/10_theory_analysis.R
# has finished its bootstrap. This step updates the package without
# refitting the smooths of 07. Run after 10.
# =====================================================================
suppressMessages(library(data.table))
if (!nzchar(Sys.getenv("SMD_ROOT"))) Sys.setenv(SMD_ROOT = getwd())
setwd(Sys.getenv("SMD_ROOT")); source(file.path("R", "config.R"))
RES <- PATHS$results
PKG <- file.path(dirname(RES), "smdMetricHarmonise")

th_file <- file.path(RES, "table_harmonisation_theory.csv")
if (!file.exists(th_file)) stop("run R/10_theory_analysis.R first: ", th_file, " not found")
th <- as.data.frame(fread(th_file))
if (!"stab_range" %in% names(th)) stop("the theory table has no bootstrap stability; finish the bootstrap in R/10_theory_analysis.R")

ht <- readRDS(file.path(RES, "harmonise_tables.rds"))
ht$theory <- th
saveRDS(ht, file.path(RES, "harmonise_tables.rds"))

harmonise_tables <- ht
save(harmonise_tables, file = file.path(PKG, "data", "harmonise_tables.rda"), compress = "xz")
cat(sprintf("package data refreshed: %d learned maps, %d translations in the theory table, built %s from %s fits\n",
            length(ht$maps), nrow(th), ht$built, format(ht$n_fits, big.mark = ",")))
