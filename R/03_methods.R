# =====================================================================
# 03_methods.R
# Method wrappers. Each returns predicted occurrence probability on the
# test set, predictions on the training set (for optimism and overfitting
# diagnostics), and a convergence flag. Predictors are standardised on
# the training set and the same transform applied to the test set, so no
# information leaks from test to train.
#
# Design method            -> implementation in this environment
#   Random Forest          -> ranger probability forest (env only)
#   Spatial Random Forest  -> ranger + coordinates + buffer distances
#   BRT                    -> gbm bernoulli
#   MaxEnt (presence only) -> glmnet on MaxEnt style features (= maxnet)
#   Spatial GAM            -> mgcv binomial, thin plate spatial smooth
#   Geostatistical (GMRF)  -> mgcv binomial, Gaussian process smooth
#   JSDM (multi species)   -> mgcv joint model, shared field + species fx
# =====================================================================

suppressMessages({
  library(ranger); library(gbm); library(glmnet); library(mgcv)
})

# scale training predictors, return transform to apply to new data
.make_scaler <- function(df, vars) {
  mu <- sapply(vars, function(v) mean(df[[v]]))
  sdv <- sapply(vars, function(v) { s <- sd(df[[v]]); if (s == 0) 1 else s })
  list(mu = mu, sdv = sdv, vars = vars)
}
.apply_scaler <- function(df, sc) {
  for (v in sc$vars) df[[v]] <- (df[[v]] - sc$mu[v]) / sc$sdv[v]
  df
}
.clip01 <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)

# ---------------------------------------------------------------- RF ---
fit_RF <- function(train, test, env_used, ...) {
  preds <- c(env_used, "year")
  fm <- as.formula(paste("y ~", paste(preds, collapse = " + ")))
  tr <- train; tr$y <- factor(tr$y, levels = c(0, 1))
  m <- ranger(fm, data = tr, probability = TRUE, num.trees = 500,
              min.node.size = 5, num.threads = 1)
  p_te <- predict(m, test)$predictions[, "1"]
  p_tr <- predict(m, train)$predictions[, "1"]
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = m, predict = function(mod, nd) .clip01(predict(mod, nd)$predictions[, "1"]))
}

# ------------------------------------------------------- SpatialRF ---
# adds coordinates plus buffer distances to a set of spread reference
# anchors (the spatial features approach for random forests)
fit_SpatialRF <- function(train, test, env_used, ...) {
  set.seed(nrow(train) * 7 + 3)
  k <- min(8L, nrow(train))
  anchors <- train[sample.int(nrow(train), k), c("lon", "lat")]
  buf <- function(df) {
    B <- sapply(seq_len(k), function(j)
      sqrt((df$lon - anchors$lon[j])^2 + (df$lat - anchors$lat[j])^2))
    colnames(B) <- paste0("buf", seq_len(k)); as.data.frame(B)
  }
  tr <- cbind(train, buf(train)); te <- cbind(test, buf(test))
  preds <- c(env_used, "year", "lon", "lat", paste0("buf", seq_len(k)))
  fm <- as.formula(paste("y ~", paste(preds, collapse = " + ")))
  tr$y <- factor(tr$y, levels = c(0, 1))
  m <- ranger(fm, data = tr, probability = TRUE, num.trees = 500,
              min.node.size = 5, num.threads = 1)
  p_te <- predict(m, te)$predictions[, "1"]
  p_tr <- predict(m, tr)$predictions[, "1"]
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = m, predict = function(mod, nd) .clip01(predict(mod, cbind(nd, buf(nd)))$predictions[, "1"]))
}

# --------------------------------------------------------------- BRT ---
fit_BRT <- function(train, test, env_used, ...) {
  preds <- c(env_used, "year")
  fm <- as.formula(paste("y ~", paste(preds, collapse = " + ")))
  m <- suppressWarnings(gbm(fm, data = train, distribution = "bernoulli",
           n.trees = 600, interaction.depth = 3, shrinkage = 0.05,
           bag.fraction = 0.75, verbose = FALSE, n.minobsinnode = 5))
  p_te <- predict(m, test, n.trees = 600, type = "response")
  p_tr <- predict(m, train, n.trees = 600, type = "response")
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = m, predict = function(mod, nd) .clip01(predict(mod, nd, n.trees = 600, type = "response")))
}

# ---------------------------------------------------------- MaxEntPO ---
# maxnet style: glmnet on MaxEnt feature classes (linear + quadratic +
# product), presence versus available (the training hauls as background)
fit_MaxEntPO <- function(train, test, env_used, ...) {
  sc <- .make_scaler(train, env_used)
  trs <- .apply_scaler(train, sc); tes <- .apply_scaler(test, sc)
  make_feats <- function(df) {
    X <- as.matrix(df[, env_used, drop = FALSE])
    Q <- X^2; colnames(Q) <- paste0(env_used, "_sq")
    P <- NULL
    if (length(env_used) >= 2) {
      cmb <- combn(env_used, 2)
      P <- sapply(seq_len(ncol(cmb)), function(j) df[[cmb[1, j]]] * df[[cmb[2, j]]])
      colnames(P) <- apply(cmb, 2, function(z) paste0(z[1], "_x_", z[2]))
    }
    cbind(X, Q, P)
  }
  # presence points versus the sampled locations as background
  Xp <- make_feats(trs[trs$y == 1, , drop = FALSE])
  Xb <- make_feats(trs)                       # available = all sampled sites
  Xg <- rbind(Xp, Xb)
  yg <- c(rep(1, nrow(Xp)), rep(0, nrow(Xb)))
  # infinitely weighted logistic regression weights (Fithian and Hastie)
  wg <- c(rep(1, nrow(Xp)), rep(100, nrow(Xb)))
  fit <- tryCatch(
    glmnet(Xg, yg, family = "binomial", weights = wg, alpha = 1,
           lambda = exp(seq(log(0.05), log(1e-4), length.out = 20))),
    error = function(e) NULL)
  if (is.null(fit)) return(list(prob = rep(mean(train$y), nrow(test)),
                                prob_train = rep(mean(train$y), nrow(train)),
                                converged = FALSE))
  lam <- fit$lambda[length(fit$lambda)]       # least regularised
  raw_te <- as.numeric(predict(fit, make_feats(tes), s = lam))
  raw_tr <- as.numeric(predict(fit, make_feats(trs), s = lam))
  # cloglog transform of the relative suitability, standard for maxnet
  p_te <- 1 - exp(-exp(raw_te - max(raw_te)))
  p_tr <- 1 - exp(-exp(raw_tr - max(raw_tr)))
  mx_te <- max(raw_te)
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = fit, predict = function(mod, nd) {
         raw <- as.numeric(predict(mod, make_feats(.apply_scaler(nd, sc)), s = lam))
         .clip01(1 - exp(-exp(raw - mx_te))) })
}

# build an mgcv smooth formula from the available environmental predictors
.gam_formula <- function(env_used, spatial = c("none", "tp", "gp"), n = Inf) {
  spatial <- match.arg(spatial)
  # Basis dimensions are capped by the training sample size. The design
  # values (k = 6 per covariate, k = 30 for the spatial smooth) exceed the
  # data at the smallest sample sizes: 46 basis functions on 30 hauls.
  # mgcv's compiled REML optimiser corrupts memory on some such draws
  # (observed on mgcv 1.9-1 and 1.9-3; the process aborts rather than
  # erroring), and a smooth with more knots than observations is in any
  # case not identifiable. The cap leaves n >= 100 unchanged.
  ke <- max(3L, min(6L,  floor(n / 8)))
  ks <- max(5L, min(30L, floor(n / 3)))
  env_terms <- sapply(env_used, function(v) sprintf("s(%s, k=%d)", v, ke))
  terms <- env_terms
  if (spatial == "tp") terms <- c(terms, sprintf("s(lon, lat, k=%d)", ks))
  if (spatial == "gp") terms <- c(terms, sprintf("s(lon, lat, bs='gp', k=%d)", ks))
  as.formula(paste("y ~", paste(terms, collapse = " + ")))
}

# --------------------------------------------------------- SpatialGAM ---
fit_SpatialGAM <- function(train, test, env_used, ...) {
  fm <- .gam_formula(env_used, spatial = "tp", n = nrow(train))
  m <- tryCatch(gam(fm, data = train, family = binomial(),
                    method = "REML", select = TRUE),
                error = function(e) NULL)
  if (is.null(m)) return(list(prob = rep(mean(train$y), nrow(test)),
                              prob_train = rep(mean(train$y), nrow(train)),
                              converged = FALSE))
  p_te <- as.numeric(predict(m, test, type = "response"))
  p_tr <- as.numeric(predict(m, train, type = "response"))
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = m, predict = function(mod, nd) .clip01(as.numeric(predict(mod, nd, type = "response"))))
}

# --------------------------------------------------------- GeostatGP ---
# Gaussian process (Matern type) spatial field, the mgcv analogue of an
# sdmTMB GMRF. Fitted with bam and covariate discretisation so the dense
# GP basis stays affordable at the sample sizes in the data limited range.
fit_GeostatGP <- function(train, test, env_used, ...) {
  # basis dimensions capped by sample size, as for the spatial GAM
  n <- nrow(train)
  ke <- max(3L, min(6L,  floor(n / 8))); ks <- max(5L, min(20L, floor(n / 3)))
  env_terms <- sapply(env_used, function(v) sprintf("s(%s, k=%d)", v, ke))
  fm <- as.formula(paste("y ~",
    paste(c(env_terms, sprintf("s(lon, lat, bs='gp', k=%d)", ks)), collapse = " + ")))
  m <- tryCatch(bam(fm, data = train, family = binomial(),
                    discrete = TRUE, method = "fREML"),
                error = function(e) NULL)
  if (is.null(m)) return(list(prob = rep(mean(train$y), nrow(test)),
                              prob_train = rep(mean(train$y), nrow(train)),
                              converged = FALSE))
  p_te <- as.numeric(predict(m, test, type = "response"))
  p_tr <- as.numeric(predict(m, train, type = "response"))
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = m, predict = function(mod, nd) .clip01(as.numeric(predict(mod, nd, type = "response"))))
}

# ---------------------------------------------------------- JointSDM ---
# multi species joint model. Fits all co occurring selected species on the
# focal split's training hauls, sharing a global environmental response
# and spatial field with species specific deviations, then predicts the
# focal species. This is the borrowing strength mechanism a JSDM provides.
# co_occ_train / co_occ_test: data.frames of the other species' 0/1 labels
# aligned to train / test rows (columns = species names, focal excluded).
fit_JointSDM <- function(train, test, env_used, co_occ_train, co_occ_test,
                         focal, ...) {
  # stack focal + co species into long form
  sp_cols <- c(focal, setdiff(colnames(co_occ_train), focal))
  stack_one <- function(base_df, occ_df, sp) {
    d <- base_df[, c("lon", "lat", "year", env_used), drop = FALSE]
    d$y <- if (sp == focal) base_df$y else occ_df[[sp]]
    d$species <- sp
    d
  }
  long_tr <- do.call(rbind, lapply(sp_cols, function(sp)
    stack_one(train, co_occ_train, sp)))
  long_tr$species <- factor(long_tr$species, levels = sp_cols)

  env_terms <- sapply(env_used, function(v) sprintf("s(%s, k=6)", v))
  fm <- as.formula(paste(
    "y ~ s(species, bs='re') +", paste(env_terms, collapse = " + "),
    "+ s(lon, lat, k=30)"))
  m <- tryCatch(bam(fm, data = long_tr, family = binomial(),
                    discrete = TRUE, method = "fREML"),
                error = function(e) NULL)
  if (is.null(m)) return(list(prob = rep(mean(train$y), nrow(test)),
                              prob_train = rep(mean(train$y), nrow(train)),
                              converged = FALSE))
  newd_te <- test[, c("lon", "lat", "year", env_used), drop = FALSE]
  newd_te$species <- factor(focal, levels = sp_cols)
  newd_tr <- train[, c("lon", "lat", "year", env_used), drop = FALSE]
  newd_tr$species <- factor(focal, levels = sp_cols)
  p_te <- as.numeric(predict(m, newd_te, type = "response"))
  p_tr <- as.numeric(predict(m, newd_tr, type = "response"))
  list(prob = .clip01(p_te), prob_train = .clip01(p_tr), converged = TRUE,
       model = m, predict = function(mod, nd) {
         nd2 <- nd[, c("lon", "lat", "year", env_used), drop = FALSE]
         nd2$species <- factor(focal, levels = sp_cols)
         .clip01(as.numeric(predict(mod, nd2, type = "response"))) })
}

METHOD_FUNS <- list(
  RF = fit_RF, SpatialRF = fit_SpatialRF, BRT = fit_BRT,
  MaxEntPO = fit_MaxEntPO, SpatialGAM = fit_SpatialGAM,
  GeostatGP = fit_GeostatGP, JointSDM = fit_JointSDM
)
