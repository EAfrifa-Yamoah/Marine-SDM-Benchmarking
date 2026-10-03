# =====================================================================
# 04_metrics.R
# All performance metrics computed on every fit, from the observed binary
# outcome y and the predicted occurrence probability p.
#
# Discrimination : AUC, TSS (at the prevalence matched threshold)
# Regression     : RMSE, MAE, R2 (predicted probability vs binary outcome)
# Rank           : Spearman
# Calibration    : Brier, calibration intercept and slope (logistic
#                  recalibration on the evaluation set)
# =====================================================================

# fast Mann Whitney AUC
.auc <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# TSS at prevalence matched threshold
.tss <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  prev <- mean(y)
  thr <- as.numeric(quantile(p, probs = 1 - prev, names = FALSE))
  pred <- as.integer(p >= thr)
  sens <- sum(pred == 1 & y == 1) / n1
  spec <- sum(pred == 0 & y == 0) / n0
  sens + spec - 1
}

.calibration <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  if (sd(lp) < 1e-8) return(c(NA_real_, NA_real_))   # no spread to calibrate
  fit <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial())),
                  error = function(e) NULL)
  if (is.null(fit) || any(is.na(coef(fit)))) return(c(NA_real_, NA_real_))
  slope <- unname(coef(fit)[2])
  # calibration intercept: intercept of a model with the linear predictor
  # entered as an offset (calibration in the large)
  fit0 <- tryCatch(suppressWarnings(glm(y ~ offset(lp), family = binomial())),
                   error = function(e) NULL)
  icpt <- if (!is.null(fit0)) unname(coef(fit0)[1]) else unname(coef(fit)[1])
  c(icpt, slope)                                       # unnamed: intercept, slope
}

compute_metrics <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  ok <- is.finite(p)
  y <- y[ok]; p <- p[ok]
  if (length(y) < 3 || length(unique(y)) < 2) {
    return(setNames(rep(NA_real_, 9),
      c("AUC","TSS","RMSE","MAE","R2","Spearman","Brier",
        "calib_intercept","calib_slope")))
  }
  resid <- y - p
  rmse <- sqrt(mean(resid^2))
  mae  <- mean(abs(resid))
  ss_tot <- sum((y - mean(y))^2)
  r2   <- if (ss_tot > 0) 1 - sum(resid^2) / ss_tot else NA_real_
  sp   <- suppressWarnings(cor(p, y, method = "spearman"))
  brier <- mean((p - y)^2)
  cal  <- .calibration(y, p)
  c(AUC = .auc(y, p), TSS = .tss(y, p), RMSE = rmse, MAE = mae,
    R2 = r2, Spearman = sp, Brier = brier,
    calib_intercept = cal[1], calib_slope = cal[2])
}
