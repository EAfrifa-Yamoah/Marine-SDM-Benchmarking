# =====================================================================
# 01_preprocess.R
# Turn raw FISHGLOB survey records into compact modelling frames.
#
# For each survey this produces:
#   hauls    : one row per unique haul with coordinates, year, and
#              environmental predictors (cleaned and, for temperature,
#              lightly imputed within survey)
#   occ      : haul_id x selected species presence/absence matrix
#   species  : selection table with prevalence, niche breadth, trophic
#              proxy and prevalence band
# and writes one RDS per survey to data_processed/.
# =====================================================================

suppressMessages({
  library(data.table)
})
source("/home/claude/marineSDMbench/R/config.R")
set.seed(GLOBAL_SEED)

# ---- coarse trophic guild from taxonomic order -----------------------
# A documented approximation used only for species stratification; the
# full study would join quantitative trophic levels from FishBase.
trophic_from_order <- function(order) {
  ord <- tolower(ifelse(is.na(order), "", order))
  high <- c("gadiformes", "lophiiformes", "anguilliformes", "squaliformes",
            "carcharhiniformes", "rajiformes", "torpediniformes",
            "myliobatiformes", "scorpaeniformes", "merlucciiformes",
            "beloniformes", "aulopiformes")
  low  <- c("clupeiformes", "osmeriformes", "argentiniformes",
            "atheriniformes", "mugiliformes", "stomiiformes")
  out <- rep(2L, length(ord))               # default: mid trophic
  out[ord %in% high] <- 3L                  # predators / higher trophic
  out[ord %in% low]  <- 1L                  # planktivores / lower trophic
  out
}

# ---- standardised niche breadth (Levins' B on depth) -----------------
# High values indicate a broad depth niche (generalist); low values a
# narrow niche (specialist). Computed from the depth distribution of a
# species' occurrences relative to available depth bins.
niche_breadth_depth <- function(depth_at_presence, depth_breaks) {
  if (length(depth_at_presence) < 5) return(NA_real_)
  b   <- cut(depth_at_presence, breaks = depth_breaks, include.lowest = TRUE)
  p   <- as.numeric(table(b)); p <- p / sum(p)
  p   <- p[p > 0]
  B   <- 1 / sum(p^2)
  nb  <- length(depth_breaks) - 1L
  (B - 1) / (nb - 1)                        # standardised to [0,1]
}

clip_coords <- function(dt) {
  qlon <- quantile(dt$lon, c(0.005, 0.995), na.rm = TRUE)
  qlat <- quantile(dt$lat, c(0.005, 0.995), na.rm = TRUE)
  dt[lon >= qlon[1] & lon <= qlon[2] & lat >= qlat[1] & lat <= qlat[2]]
}

pick_species <- function(cand, n_target) {
  # cand: data.table with prevalence, band, n_pres, niche, trophic
  per_band <- max(1L, round(n_target / length(PREV_BANDS)))
  chosen <- list()
  for (bd in names(PREV_BANDS)) {
    pool <- cand[band == bd][order(-n_pres)]
    if (nrow(pool) == 0) next
    # spread niche breadth: take species at low, mid, high niche where possible
    k <- min(per_band, nrow(pool))
    if (k == 1) {
      chosen[[bd]] <- pool[1]
    } else {
      idx <- unique(round(seq(1, nrow(pool), length.out = k)))
      chosen[[bd]] <- pool[idx]
    }
  }
  out <- rbindlist(chosen)
  if (nrow(out) > n_target) out <- out[order(band, -n_pres)][1:n_target]
  out
}

process_one <- function(code) {
  cat(sprintf("\n[preprocess] %s\n", code))
  e <- new.env()
  load(file.path(PATHS$raw, paste0(code, "_clean.RData")), envir = e)
  d <- as.data.table(get("data", envir = e))
  setnames(d, c("longitude", "latitude"), c("lon", "lat"))

  # some surveys store depth or temperature as character; coerce the
  # numeric predictors up front so completeness is judged correctly
  for (v in c("depth", "sbt", "sst", "lon", "lat", "year")) {
    if (v %in% names(d) && !is.numeric(d[[v]])) {
      suppressWarnings(d[, (v) := as.numeric(as.character(get(v)))])
    }
  }

  # --- haul table ---
  hauls <- unique(d[, .(haul_id, lon, lat, year, depth, sbt, sst)])
  hauls <- clip_coords(hauls)
  hauls <- hauls[!is.na(depth) & !is.na(lon) & !is.na(lat) & !is.na(year)]

  # temperature: keep a covariate only where mostly observed, then mean
  # impute residual gaps within survey
  env_used <- ENV_CORE
  for (v in ENV_OPTIONAL) {
    complete_frac <- mean(!is.na(hauls[[v]]))
    if (complete_frac >= TEMP_KEEP_THRESHOLD) {
      mu <- mean(hauls[[v]], na.rm = TRUE)
      hauls[is.na(get(v)), (v) := mu]
      env_used <- c(env_used, v)
    } else {
      hauls[, (v) := NULL]
    }
  }
  cat(sprintf("  hauls kept = %d | env predictors = %s\n",
              nrow(hauls), paste(env_used, collapse = ", ")))

  # --- species prevalence and niche over the cleaned haul set ---
  d2 <- d[haul_id %in% hauls$haul_id]
  n_hauls <- nrow(hauls)
  depth_breaks <- quantile(hauls$depth, probs = seq(0, 1, length.out = 11),
                           na.rm = TRUE)
  depth_breaks[1] <- depth_breaks[1] - 1e-6

  # presence records: a species is present in a haul if it has any positive
  # record there (num or wgt > 0, or simply listed)
  d2[, pres := 1L]
  sp_hauls <- unique(d2[, .(accepted_name, haul_id, order)])
  # attach depth for niche computation
  sp_hauls <- merge(sp_hauls, hauls[, .(haul_id, depth)], by = "haul_id")

  sp_stats <- sp_hauls[, .(
    n_pres  = uniqueN(haul_id),
    order   = order[1]
  ), by = accepted_name]
  sp_stats[, prevalence := n_pres / n_hauls]
  sp_stats[, trophic := trophic_from_order(order)]

  # niche breadth per species
  nb <- sp_hauls[, .(niche = niche_breadth_depth(depth, depth_breaks)),
                 by = accepted_name]
  sp_stats <- merge(sp_stats, nb, by = "accepted_name")

  # assign prevalence band
  sp_stats[, band := NA_character_]
  for (bd in names(PREV_BANDS)) {
    rng <- PREV_BANDS[[bd]]
    sp_stats[prevalence >= rng[1] & prevalence < rng[2], band := bd]
  }
  cand <- sp_stats[!is.na(band) & n_pres >= MIN_TOTAL_PRESENCES &
                     !is.na(niche)]
  cat(sprintf("  candidate species (banded, >=%d presences) = %d\n",
              MIN_TOTAL_PRESENCES, nrow(cand)))

  sel <- pick_species(cand, N_SPECIES_PER_DATASET)
  cat(sprintf("  selected %d species:\n", nrow(sel)))
  for (i in seq_len(nrow(sel))) {
    cat(sprintf("    %-28s prev=%.3f band=%-12s niche=%.2f trophic=%d n_pres=%d\n",
                sel$accepted_name[i], sel$prevalence[i], sel$band[i],
                sel$niche[i], sel$trophic[i], sel$n_pres[i]))
  }

  # --- occurrence matrix for selected species ---
  sel_names <- sel$accepted_name
  pres_tab <- unique(d2[accepted_name %in% sel_names, .(haul_id, accepted_name)])
  pres_tab[, val := 1L]
  occ <- dcast(pres_tab, haul_id ~ accepted_name, value.var = "val",
               fill = 0L)
  # ensure all selected species present as columns and all hauls represented
  miss_cols <- setdiff(sel_names, names(occ))
  for (mc in miss_cols) occ[[mc]] <- 0L
  occ <- merge(hauls[, .(haul_id)], occ, by = "haul_id", all.x = TRUE)
  for (sn in sel_names) occ[is.na(get(sn)), (sn) := 0L]

  list(code = code, hauls = hauls, occ = occ, species = sel,
       env_used = env_used, n_hauls = n_hauls)
}

# ---- run all datasets ------------------------------------------------
all_species <- list()
for (code in names(DATASETS)) {
  res <- process_one(code)
  saveRDS(res, file.path(PATHS$processed, paste0(code, "_processed.rds")))
  si <- as.data.table(res$species)
  si[, dataset := code]
  all_species[[code]] <- si
  rm(res); gc(FALSE)
}
species_master <- rbindlist(all_species, fill = TRUE)
fwrite(species_master, file.path(PATHS$results, "selected_species.csv"))

cat("\n[preprocess] done. selected species across datasets:\n")
print(species_master[, .(dataset, accepted_name, band, prevalence = round(prevalence,3),
                         niche = round(niche,2), trophic, n_pres)])
