# =====================================================================
# variogram.R
# Practical range of an exponential variogram of the residuals of an
# environmental additive model, used to express the separation of a held
# out block from its training locations in correlation range units
# (Corollary 4). Sourced by 05b_oracle.R and 05d_practical_range.R.
#
# Model: gamma(h) = c0 + c1 {1 - exp(-3 h / A)}, nugget c0 >= 0, partial sill
# c1 >= 0, practical range A (the lag at which 95% of the sill is reached),
# fitted by least squares to the binned empirical variogram on lags up to
# the median distance between hauls.
#
# The original oracle stage minimised this criterion with one local
# optimiser start. For many species the criterion is flat in A (no spatial
# structure in the residuals, so A is not identified and the optimiser
# returned its starting value) or keeps falling as A grows (no sill within
# the survey, and A ran away to many times the survey extent). Given A the
# model is linear in (c0, c1), so here the criterion is profiled exactly
# over a grid of A from the smallest lag examined to the largest distance
# between the hauls used, with non negative sills, and each species is
# classed by where the profile minimum lies:
#   identified     minimum inside the grid with a partial sill of at least
#                  5% of the total sill; A is reported
#   no structure   partial sill below 5% of the total: residual spatial
#                  autocorrelation not detectable at the lags examined
#   below lags     the fitted variogram reaches 95% of its sill before the
#                  second lag bin, so its rise is not resolved: residual
#                  autocorrelation shorter than the lags examined
#   no sill        minimum at the largest distance: no sill within the survey
# Only identified ranges are used to scale separation. The optimiser result
# of the original stage is also returned (range_optim_km) for comparison.
# =====================================================================
.exp_profile <- function(hh, gh, A) {
  # least squares fit of gamma = c0 + c1 (1 - exp(-3 h / A)) with c0, c1 >= 0
  x1 <- 1 - exp(-3 * hh / A)
  fits <- list(c(mean(gh), 0))                                         # nugget only
  c1 <- max(0, sum(x1 * gh) / sum(x1^2)); fits[[2]] <- c(0, c1)        # sill only
  X <- cbind(1, x1); b <- tryCatch(qr.solve(X, gh), error = function(e) c(-1, -1))
  if (all(b >= 0)) fits[[3]] <- b
  sse <- vapply(fits, function(b) sum((gh - b[1] - b[2] * x1)^2), 0)
  k <- which.min(sse); c(sse = sse[k], c0 = fits[[k]][1], c1 = fits[[k]][2])
}

variogram_fit <- function(df, env_used, max_n = 1500L, n_grid = 400L, min_sill_share = 0.05) {
  na <- list(range_km = NA_real_, regime = NA_character_, sill_share = NA_real_, lag_min_km = NA_real_,
             lag_max_km = NA_real_, dist_max_km = NA_real_, range_optim_km = NA_real_)
  if (nrow(df) > max_n) df <- df[sample.int(nrow(df), max_n), ]
  fm <- as.formula(paste("y ~", paste(sprintf("s(%s, k=6)", env_used), collapse = " + ")))
  m <- tryCatch(mgcv::gam(fm, data = df, family = binomial()), error = function(e) NULL)
  if (is.null(m)) return(na)
  res <- residuals(m, type = "pearson")
  latm <- mean(df$lat); xy <- cbind(df$lon * 111.32 * cos(latm * pi / 180), df$lat * 110.57)
  D <- as.matrix(dist(xy)); G <- outer(res, res, "-")^2 / 2
  ut <- upper.tri(D); h <- D[ut]; g <- G[ut]; dist_max <- max(h)
  keep <- h > 0 & h < quantile(h, 0.5); h <- h[keep]; g <- g[keep]
  bins <- cut(h, breaks = 15); gh <- tapply(g, bins, mean); hh <- tapply(h, bins, mean)
  ok <- is.finite(gh) & is.finite(hh); gh <- as.numeric(gh[ok]); hh <- as.numeric(hh[ok])
  if (length(gh) < 5) return(na)
  # the original single start optimiser, kept for comparison with oracle_intrinsic.csv
  f <- function(par) sum((gh - (par[1] + par[2] * (1 - exp(-hh / par[3]))))^2)
  o <- tryCatch(optim(c(min(gh), diff(range(gh)), mean(hh)), f, method = "L-BFGS-B",
                      lower = c(0, 1e-6, 1e-3)), error = function(e) NULL)
  r_optim <- if (is.null(o)) NA_real_ else 3 * o$par[3]
  # exact profile over the practical range
  A <- exp(seq(log(min(hh)), log(dist_max), length.out = n_grid))
  P <- t(vapply(A, function(a) .exp_profile(hh, gh, a), c(sse = 0, c0 = 0, c1 = 0)))
  k <- which.min(P[, "sse"]); share <- P[k, "c1"] / max(P[k, "c0"] + P[k, "c1"], .Machine$double.eps)
  regime <- if (share < min_sill_share) "no structure" else if (A[k] <= hh[2]) "below lags" else
            if (k == n_grid) "no sill" else "identified"
  list(range_km = if (regime == "identified") A[k] else NA_real_, regime = regime, sill_share = share,
       lag_min_km = min(hh), lag_max_km = max(h), dist_max_km = dist_max, range_optim_km = r_optim)
}

# the identified practical range alone (km; NA when not identified), as used by the oracle stage
practical_range <- function(df, env_used, max_n = 1500L) variogram_fit(df, env_used, max_n)$range_km
