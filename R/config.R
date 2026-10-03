# =====================================================================
# config.R
# Benchmark harmonisation protocol for data limited marine SDM
# Follow up simulation operationalising Section 7.1 of the MEE review
#
# This file encodes the experimental design exactly as specified in the
# study design document, plus a reduced PILOT scale that runs on a single
# workstation. The full study runs by setting SCALE <- "full" and
# executing on an HPC allocation. Everything downstream reads from here,
# so the same code produces the pilot and the full run.
# =====================================================================

# ---- run scale -------------------------------------------------------
# "pilot"  : single core demonstration on real FISHGLOB data
# "full"   : the 1.26 million fit design in the study document (HPC)
SCALE <- Sys.getenv("SMD_SCALE", unset = "pilot")

# ---- design factors (identical across scales) ------------------------
# Factor 1: sample size, spanning the data limited regime
SAMPLE_SIZES   <- c(30, 50, 100, 200, 500)

# Factor 2: spatial coverage of the training observations
#   clustered    = observations concentrated in ~25% of the spatial extent
#   intermediate = ~50% of the extent
#   distributed  = full 100% extent, spread approximately uniformly
COVERAGE_LEVELS <- c("clustered", "intermediate", "distributed")
COVERAGE_FRAC   <- c(clustered = 0.25, intermediate = 0.50, distributed = 1.00)

# Factor 3: inference target
#   conditional = test points drawn from within the training spatial and
#                 temporal extent (interpolation; optimistic analogue)
#   marginal    = test points drawn from a contiguous spatial block held
#                 out before subsampling (extrapolation; honest analogue)
INFERENCE_TARGETS <- c("conditional", "marginal")

# 5 sample sizes x 3 coverage x 2 targets = 30 design cells per dataset
stopifnot(length(SAMPLE_SIZES) * length(COVERAGE_LEVELS) *
            length(INFERENCE_TARGETS) == 30)

# ---- source datasets -------------------------------------------------
# Four harmonised FISHGLOB survey regions spanning distinct large marine
# ecosystems. Three of the four datasets named in the design are FISHGLOB
# constituents: NS-IBTS is an ICES DATRAS survey, and EBS sits inside the
# surveyjoin Northeast Pacific compilation. GMEX supplies a subtropical
# contrast standing in for the very different ecological context that
# SARDARA provides in the full design.
DATASETS <- c(
  EBS       = "Eastern Bering Sea (Northeast Pacific, NOAA)",
  `NS-IBTS` = "North Sea IBTS (Northeast Atlantic, ICES DATRAS)",
  GMEX      = "Gulf of Mexico (subtropical, SEAMAP)",
  NEUS      = "Northeast US shelf (Northwest Atlantic, NEFSC)"
)

# ---- environmental predictors ----------------------------------------
# depth and coordinates are complete across all four surveys; bottom and
# surface temperature are retained only where they are largely observed
# (they are >75% missing in NS-IBTS). Temperature is kept for a survey
# only when its completeness exceeds TEMP_KEEP_THRESHOLD, and residual
# gaps are mean imputed within the training set to avoid leakage.
ENV_CORE          <- c("depth")          # universal environmental driver
ENV_OPTIONAL      <- c("sbt", "sst")     # kept where well observed
TEMP_KEEP_THRESHOLD <- 0.60              # keep covariate if >=60% complete
SPATIAL_PREDICTORS  <- c("lon", "lat")   # used by spatial methods only
TEMPORAL_PREDICTOR  <- "year"

# ---- species selection -----------------------------------------------
# Prevalence bands follow the design: rare <5%, intermediate 5-25%,
# common >25% of hauls. Species must have enough total presences that a
# small subsample can still contain occurrences.
PREV_BANDS <- list(
  rare         = c(0.02, 0.05),
  intermediate = c(0.05, 0.25),
  common       = c(0.25, 0.75)
)
MIN_TOTAL_PRESENCES <- 300   # over the whole survey haul table

# ---- scale dependent knobs -------------------------------------------
if (SCALE == "full") {
  N_SPECIES_PER_DATASET <- 15    # ~12-20 in the design
  N_REPLICATES          <- 100   # subsample replicates per design cell
  TEST_SIZE             <- 400   # honest test set size per fit
} else {                          # pilot
  N_SPECIES_PER_DATASET <- 3     # one rare, one intermediate, one common
  N_REPLICATES          <- 4     # replicates per design cell (pilot)
  TEST_SIZE             <- 300
}
# Optional local override of the replicate count (for example 25 or 50 for a
# faster first pass). Set SMD_REPLICATES.
.rep_override <- Sys.getenv("SMD_REPLICATES", unset = "")
if (nzchar(.rep_override)) N_REPLICATES <- as.integer(.rep_override)

# ---- local parallelism -----------------------------------------------
# Worker processes for the parallel driver: all but one core by default,
# override with SMD_CORES. A PSOCK cluster runs on Windows, macOS and Linux.
N_CORES <- as.integer(Sys.getenv("SMD_CORES",
  unset = as.character(max(1L, parallel::detectCores() - 1L))))

# ---- methods benchmarked ---------------------------------------------
# Seven methods spanning the design's 2x2x2 structure: spatial structure
# (yes/no) x algorithm family x single/multi species. Where a package
# named in the design is not installable in this environment, an
# available package implementing the same modelling principle is used and
# the substitution is recorded in method_registry below.
METHODS <- c("RF", "SpatialRF", "BRT", "MaxEntPO",
             "SpatialGAM", "GeostatGP", "JointSDM")

method_registry <- data.frame(
  method   = METHODS,
  family   = c("tree", "tree", "boosting", "presence_only",
               "gam", "geostatistical", "jsdm"),
  spatial  = c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
  multisp  = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  design_ref = c("Random Forest (ranger)",
                 "Spatial Random Forest (ranger + spatial buffers)",
                 "BRT (gbm / xgboost)",
                 "MaxEnt (presence only baseline)",
                 "Spatial GAM (mgcv spatial smooths)",
                 "sdmTMB (GMRF / SPDE geostatistical)",
                 "tinyVAST / HMSC (multi species geostatistical / JSDM)"),
  implemented_with = c("ranger", "ranger + coordinate and buffer features",
                       "gbm", "penalised presence background GLM (glmnet)",
                       "mgcv thin plate spatial smooth",
                       "mgcv Gaussian process smooth + fields kriging check",
                       "mgcv joint model, shared spatial field + species effects"),
  stringsAsFactors = FALSE
)

# ---- metrics computed on every fit -----------------------------------
METRIC_NAMES <- c("AUC", "TSS", "RMSE", "MAE", "R2",
                  "Spearman", "Brier", "calib_intercept", "calib_slope")

# ---- reproducibility -------------------------------------------------
GLOBAL_SEED <- 20260701

# ---- paths -----------------------------------------------------------
# Resolved relative to the project root so the pipeline runs on any machine:
# the working directory by default, or SMD_ROOT if set.
PROJECT_ROOT <- Sys.getenv("SMD_ROOT", unset = getwd())
PATHS <- list(
  raw        = file.path(PROJECT_ROOT, "data_raw"),
  processed  = file.path(PROJECT_ROOT, "data_processed"),
  results    = file.path(PROJECT_ROOT, "results"),
  figures    = file.path(PROJECT_ROOT, "figures"),
  logs       = file.path(PROJECT_ROOT, "logs"),
  pkg_data   = file.path(PROJECT_ROOT, "smdMetricHarmonise", "data")
)
for (.p in PATHS) dir.create(.p, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[config] scale = %s | %d species/dataset | %d replicates | %d methods\n",
            SCALE, N_SPECIES_PER_DATASET, N_REPLICATES, length(METHODS)))
cat(sprintf("[config] design cells per dataset = %d | datasets = %d\n",
            30L, length(DATASETS)))
