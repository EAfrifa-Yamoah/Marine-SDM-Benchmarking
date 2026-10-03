# =====================================================================
# 02_subsample.R
# Stratified factorial subsampling.
#
# Draws a training subsample and an honest test set from a processed
# survey for a given focal species and design cell. The three design
# factors are varied independently:
#
#   sample size  -> exactly n training hauls
#   coverage     -> the training hauls occupy ~25% / ~50% / 100% of the
#                   available spatial extent (clustered / intermediate /
#                   distributed). Distributed hauls are spread across the
#                   extent by grid stratification.
#   target       -> conditional: test hauls lie within the training
#                     region and year range (interpolation)
#                   marginal:    test hauls lie in a contiguous spatial
#                     block held out before subsampling (extrapolation)
#
# A subsample is rejected (and the failure recorded) when the training or
# test occurrence vector is degenerate, so downstream failure mode
# analysis sees genuine data limited breakdowns.
# =====================================================================

suppressMessages(library(data.table))

# grid stratified indices to spread a sample across a bounding box
.spread_grid <- function(lon, lat, n) {
  g <- max(2L, ceiling(sqrt(n)))
  bx <- cut(lon, breaks = g, labels = FALSE, include.lowest = TRUE)
  by <- cut(lat, breaks = g, labels = FALSE, include.lowest = TRUE)
  cell <- paste(bx, by)
  idx_by_cell <- split(seq_along(cell), cell)
  idx_by_cell <- idx_by_cell[sample(length(idx_by_cell))]  # shuffle cells
  chosen <- integer(0); k <- 1L
  # round robin across occupied cells until n reached
  while (length(chosen) < n) {
    took_any <- FALSE
    for (cc in seq_along(idx_by_cell)) {
      pool <- setdiff(idx_by_cell[[cc]], chosen)
      if (length(pool) > 0) {
        chosen <- c(chosen, pool[sample.int(length(pool), 1L)])
        took_any <- TRUE
        if (length(chosen) >= n) break
      }
    }
    if (!took_any) break
    k <- k + 1L; if (k > n + 5L) break
  }
  head(chosen, n)
}

# a random axis aligned window covering ~frac of a bounding box that
# contains at least min_pts hauls; returns row indices inside the window
.coverage_window <- function(lon, lat, frac, min_pts, max_try = 30L) {
  if (frac >= 0.999) return(seq_along(lon))
  lon_rng <- range(lon); lat_rng <- range(lat)
  s <- sqrt(frac)                                   # side scale
  wlon <- diff(lon_rng) * s; wlat <- diff(lat_rng) * s
  for (t in seq_len(max_try)) {
    ax <- runif(1, lon_rng[1], lon_rng[2] - wlon)
    ay <- runif(1, lat_rng[1], lat_rng[2] - wlat)
    inside <- which(lon >= ax & lon <= ax + wlon &
                      lat >= ay & lat <= ay + wlat)
    if (length(inside) >= min_pts) return(inside)
  }
  # fall back to nearest neighbours around a random anchor
  anchor <- sample.int(length(lon), 1L)
  d <- (lon - lon[anchor])^2 + (lat - lat[anchor])^2
  order(d)[seq_len(max(min_pts, 1L))]
}

# choose the held out block for the marginal target; rotate the edge by
# replicate so different geography is extrapolated to across replicates
.holdout_block <- function(lon, lat, rep_id, block_frac = 0.20) {
  edge <- c("east", "west", "north", "south")[(rep_id %% 4L) + 1L]
  test <- switch(edge,
    east  = lon >= quantile(lon, 1 - block_frac, names = FALSE),
    west  = lon <= quantile(lon, block_frac,     names = FALSE),
    north = lat >= quantile(lat, 1 - block_frac, names = FALSE),
    south = lat <= quantile(lat, block_frac,     names = FALSE))
  list(test = which(test), train_pool = which(!test), edge = edge)
}

# main entry: build one train/test split for a design cell
make_split <- function(hauls, y, env_used, n, coverage, target,
                       rep_id, test_size = 300L,
                       min_pres = 3L, min_abs = 3L, max_retry = 10L) {
  frac <- COVERAGE_FRAC[[coverage]]
  lon <- hauls$lon; lat <- hauls$lat
  N <- length(lon)
  keepcols <- c("haul_id", "lon", "lat", "year", env_used)

  build_df <- function(idx) {
    df <- as.data.frame(hauls[idx, ..keepcols])
    df$y <- y[idx]
    df
  }

  for (attempt in seq_len(max_retry)) {
    if (target == "marginal") {
      hb <- .holdout_block(lon, lat, rep_id + attempt)
      pool <- hb$train_pool
      test_pool <- hb$test
      # coverage window within the training pool
      win_local <- .coverage_window(lon[pool], lat[pool], frac,
                                    min_pts = n + 5L)
      win <- pool[win_local]
      if (length(win) < n) next
      if (coverage == "distributed") {
        tr_local <- .spread_grid(lon[win], lat[win], n)
        train_idx <- win[tr_local]
      } else {
        train_idx <- win[sample.int(length(win), n)]
      }
      te_n <- min(test_size, length(test_pool))
      test_idx <- test_pool[sample.int(length(test_pool), te_n)]
      # interpolation test set for the same fit contrast: drawn from the
      # training window, disjoint from the training set, so that one fitted
      # predictor can be evaluated under both inference targets
      rest_in <- setdiff(win, train_idx)
      interp_idx <- if (length(rest_in) >= 5L)
        rest_in[sample.int(length(rest_in), min(test_size, length(rest_in)))] else integer(0)
    } else {                                   # conditional
      interp_idx <- integer(0)
      win <- .coverage_window(lon, lat, frac, min_pts = n + test_size + 5L)
      if (length(win) < n + 5L) next
      if (coverage == "distributed") {
        tr_local <- .spread_grid(lon[win], lat[win], n)
        train_idx <- win[tr_local]
      } else {
        train_idx <- win[sample.int(length(win), n)]
      }
      rest <- setdiff(win, train_idx)
      if (length(rest) < 1L) next
      te_n <- min(test_size, length(rest))
      test_idx <- rest[sample.int(length(rest), te_n)]
    }

    ytr <- y[train_idx]; yte <- y[test_idx]
    if (sum(ytr) >= min_pres && sum(1 - ytr) >= min_abs &&
        sum(yte) >= 1 && sum(1 - yte) >= 1) {
      # design geometry of the draw (Proposition 5): fill distance of the
      # training design over the sampled extent outside the block, and the
      # separation of the block test points from the training locations
      geo <- .draw_geometry(lon, lat, train_idx, test_idx,
                            if (target == "marginal") pool else seq_len(N))
      return(list(ok = TRUE,
                  train = build_df(train_idx),
                  test  = build_df(test_idx),
                  interp = if (length(interp_idx)) build_df(interp_idx) else NULL,
                  fill_km = geo$fill_km, separation_km = geo$sep_km,
                  reason = NA_character_,
                  edge = if (target == "marginal") hb$edge else NA_character_))
    }
  }
  list(ok = FALSE, train = NULL, test = NULL,
       reason = "degenerate_occurrence", edge = NA_character_)
}


# fill distance and block separation of a draw, in kilometres. The sampled
# extent is represented by a subsample of the haul locations outside the
# block; distances use a local equirectangular approximation.
.draw_geometry <- function(lon, lat, train_idx, test_idx, domain_idx,
                           max_dom = 2000L) {
  latm <- mean(lat[train_idx])
  kx <- 111.32 * cos(latm * pi / 180); ky <- 110.57
  tr <- cbind(lon[train_idx] * kx, lat[train_idx] * ky)
  dom <- if (length(domain_idx) > max_dom)
    domain_idx[sample.int(length(domain_idx), max_dom)] else domain_idx
  dm <- cbind(lon[dom] * kx, lat[dom] * ky)
  te <- cbind(lon[test_idx] * kx, lat[test_idx] * ky)
  nn <- function(a, b) {   # min distance from each row of a to rows of b
    d2 <- outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2
    sqrt(apply(d2, 1, min))
  }
  list(fill_km = max(nn(dm, tr)), sep_km = min(nn(te, tr)))
}
