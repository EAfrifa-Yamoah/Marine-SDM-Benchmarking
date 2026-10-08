# Numerical verification of the identities and inequalities stated in the
# theoretical foundation. Base R only. Run: Rscript verify_identities.R
#
# Sections 1 to 10 are the script of the theoretical foundation (16 September
# 2026), restored here because the manuscript (Appendix C) describes it as
# included with the analysis code; their output is unchanged. Sections 11 to
# 14 were added with the generalised spatial results (Lemma 6, Corollary 6)
# and check the results of Appendix A.5 without compact support, with a
# Matern kernel of the kind the benchmark's Gaussian process smooth uses.
# Section 14 needs mgcv and is skipped when it is not installed.
set.seed(20260916)
auc_fun <- function(p, y) {
  # Mann Whitney form with half credit for ties (identical to the population
  # functional in Definition 1)
  r <- rank(p)                       # midranks
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n0 * n1)
}
brier <- function(p, y) mean((y - p)^2)
mae   <- function(p, y) mean(abs(y - p))
r2ss  <- function(p, y) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
tss_at <- function(p, y, t) {
  sens <- mean(p[y == 1] > t); spec <- mean(p[y == 0] <= t)
  sens + spec - 1
}

# 1 -----------------------------------------------------------------------
cat("=== Theorem 1c: Spearman(p,y) = (AUC - 1/2) * sqrt(12 n0 n1 / (n^2 - 1)) (no ties) ===\n")
worst <- 0
for (rep in 1:200) {
  n <- sample(30:500, 1)
  pi <- runif(1, 0.03, 0.6)
  y <- rbinom(n, 1, pi); if (sum(y) %in% c(0, n)) next
  p <- plogis(rnorm(n, mean = 1.2 * y, sd = 1))   # continuous, no ties
  n1 <- sum(y); n0 <- n - n1
  lhs <- cor(p, y, method = "spearman")
  rhs <- (auc_fun(p, y) - 0.5) * sqrt(12 * n0 * n1 / (n^2 - 1))
  worst <- max(worst, abs(lhs - rhs))
}
cat("max |Spearman - affine(AUC)| over 200 random draws:", format(worst, digits = 3), "\n")
cat("large n approximation 2*sqrt(3*pi*(1-pi))*(AUC-1/2), pi = 0.25, AUC = 0.80 ->",
    round(2 * sqrt(3 * 0.25 * 0.75) * 0.3, 3), "\n\n")

# 2 -----------------------------------------------------------------------
cat("=== Theorem 1b: R2 (sum of squares form) = 1 - Brier / (pi (1 - pi)) ===\n")
worst <- 0
for (rep in 1:200) {
  n <- sample(30:500, 1); pi <- runif(1, 0.03, 0.6)
  y <- rbinom(n, 1, pi); if (sum(y) %in% c(0, n)) next
  p <- plogis(rnorm(n, mean = 0.8 * y, sd = 1.5))
  pih <- mean(y)
  worst <- max(worst, abs(r2ss(p, y) - (1 - brier(p, y) / (pih * (1 - pih)))))
}
cat("max |R2 - (1 - Brier/(pi(1-pi)))| over 200 draws:", format(worst, digits = 3), "\n")
cat("example: AUC unchanged and R2 reduced by an over confident recalibration\n")
n <- 200; y <- rbinom(n, 1, 0.2); lp <- -1.4 + 1.0 * y + rnorm(n, 0, 1)
p_cal <- plogis(lp); p_over <- plogis(4 * lp)       # slope 1/4 on recalibration
cat(" AUC calibrated:", round(auc_fun(p_cal, y), 3), " AUC over confident:", round(auc_fun(p_over, y), 3), "\n")
cat(" R2 calibrated:", round(r2ss(p_cal, y), 3), " R2 over confident:", round(r2ss(p_over, y), 3), "\n\n")

# 3 -----------------------------------------------------------------------
cat("=== Proposition 1: MAE^2 <= Brier <= MAE, and MAE = 2 Brier for a calibrated predictor ===\n")
ok <- TRUE
for (rep in 1:200) {
  n <- 300; y <- rbinom(n, 1, runif(1, 0.05, 0.6)); p <- runif(n)
  b <- brier(p, y); m <- mae(p, y)
  ok <- ok && (m^2 <= b + 1e-12) && (b <= m + 1e-12)
}
cat("MAE^2 <= Brier <= MAE held in all 200 draws:", ok, "\n")
n <- 2e5; p <- runif(n); y <- rbinom(n, 1, p)     # perfectly calibrated
cat("calibrated predictor, n = 2e5: MAE =", round(mae(p, y), 4), " 2*Brier =", round(2 * brier(p, y), 4), "\n\n")

# 4 -----------------------------------------------------------------------
cat("=== Proposition 2: AUC - 1/2 = integral of TSS(t) dF1(t); AUC - 1/2 <= max_t TSS(t) ===\n")
n <- 4000; y <- rbinom(n, 1, 0.3); p <- plogis(rnorm(n, 1.1 * y, 1))
ts <- sort(p[y == 1])
integral <- mean(sapply(ts, function(t) tss_at(p, y, t)))
cat("AUC - 1/2 =", round(auc_fun(p, y) - 0.5, 4), " mean_{t ~ F1} TSS(t) =", round(integral, 4),
    " max_t TSS =", round(max(sapply(seq(0.01, 0.99, by = 0.005), function(t) tss_at(p, y, t))), 4), "\n\n")

# 5 -----------------------------------------------------------------------
cat("=== Theorem 2: AUC invariant to monotone recalibration; Brier is not ===\n")
n <- 1000; y <- rbinom(n, 1, 0.3); lp <- -1 + 1.3 * y + rnorm(n)
slopes <- c(0.2, 0.5, 1, 2, 4)
res <- t(sapply(slopes, function(s) c(slope = s, AUC = auc_fun(plogis(s * lp), y), Brier = brier(plogis(s * lp), y))))
print(round(res, 4)); cat("\n")

# 6 -----------------------------------------------------------------------
cat("=== Theorem 3: between group variance >= ((K-1)/(2K)) * Delta^2 ===\n")
ok <- TRUE
for (rep in 1:1000) {
  K <- sample(2:6, 1); f <- rnorm(K)
  bv <- mean((f - mean(f))^2)
  D <- mean(abs(outer(f, f, "-"))[upper.tri(matrix(0, K, K))])
  ok <- ok && (bv >= (K - 1) / (2 * K) * D^2 - 1e-12)
}
cat("inequality held in 1000 random draws:", ok, "\n")
f <- c(0.3, 0.6); cat("K = 2: between variance", mean((f - mean(f))^2), " bound", 1 / 4 * 0.3^2, "(equality)\n\n")

# 7 -----------------------------------------------------------------------
cat("=== Lemmas 3 and 4: posterior variance beyond kernel range equals prior; monotone in the design ===\n")
wendland <- function(d, R) ifelse(d < R, (1 - d / R)^4 * (4 * d / R + 1), 0)
post_var <- function(S, s, R, sig2) {
  K <- wendland(as.matrix(dist(S)), R) + diag(sig2, nrow(S))
  k <- wendland(sqrt(rowSums((S - matrix(s, nrow(S), 2, byrow = TRUE))^2)), R)
  1 - drop(t(k) %*% solve(K, k))
}
S1 <- cbind(runif(60, 0, 1), runif(60, 0, 1)); S2 <- rbind(S1, cbind(runif(60, 0, 1), runif(60, 0, 1)))
s_far <- c(1.6, 0.5)
cat("prior variance 1; posterior at s_far, n = 60:", post_var(S1, s_far, 0.4, 0.05),
    " n = 120:", post_var(S2, s_far, 0.4, 0.05), "\n")
s_in <- c(0.5, 0.5)
cat("posterior at interior point, n = 60:", round(post_var(S1, s_in, 0.4, 0.05), 4),
    " n = 120:", round(post_var(S2, s_in, 0.4, 0.05), 4), "\n\n")

# 8 -----------------------------------------------------------------------
cat("=== Lemma 2: excess Brier equals L2 distance to the Bayes predictor ===\n")
n <- 2e5; x <- runif(n); pstar <- plogis(-2 + 4 * x); y <- rbinom(n, 1, pstar)
phat <- plogis(-1.5 + 3 * x)
cat("Brier(phat) - Brier(pstar) =", round(brier(phat, y) - brier(pstar, y), 5),
    " mean (phat - pstar)^2 =", round(mean((phat - pstar)^2), 5), "\n\n")

# 9 -----------------------------------------------------------------------
cat("=== Fit count ===\n")
cat("4 x 3 x 30 x 4 x 7 =", 4 * 3 * 30 * 4 * 7, "; minus 40 x 7 =", 4 * 3 * 30 * 4 * 7 - 40 * 7, "\n")

# 10 ----------------------------------------------------------------------
cat("full scale: 4 x 15 x 30 x 100 x 7 =", format(4 * 15 * 30 * 100 * 7, big.mark = ","),
    "; minus 11,747 degenerate splits x 7 =", format(4 * 15 * 30 * 100 * 7 - 11747 * 7, big.mark = ","), "\n\n")

# =========================================================================
# Added with Lemma 6 and Corollary 6: the spatial results without compact
# support. Domain D = [0, 1.5] x [0, 1], training locations in [0, 1] x [0, 1],
# held out edge block B = (1, 1.5] x [0, 1]. Matern kernel with smoothness
# 3/2 in the parameterisation of mgcv's Gaussian process smooth,
# k(d) = tau2 (1 + d / rho) exp(-d / rho), which has no compact support.
# =========================================================================
matern32 <- function(d, rho, tau2 = 1) tau2 * (1 + d / rho) * exp(-d / rho)
sqexp    <- function(d, rho, tau2 = 1) tau2 * exp(-(d / rho)^2)
dist_to  <- function(S, s) sqrt(rowSums((S - matrix(s, nrow(S), 2, byrow = TRUE))^2))
kbar     <- function(kern, h, rho) { d <- seq(h, h + 50 * rho, length.out = 5000); max(abs(kern(d, rho))) }
post_var_k <- function(kern, S, s, rho, sig2, tau2 = 1) {
  K <- kern(as.matrix(dist(S)), rho, tau2) + diag(sig2, nrow(S)); k <- kern(dist_to(S, s), rho, tau2)
  tau2 - drop(crossprod(k, solve(K, k)))
}

# 11 ----------------------------------------------------------------------
cat("=== Lemma 3 without compact support: |w_hat(s)| <= ||alpha||_1 kbar(h), v_n(s) >= tau2 - n kbar(h)^2 / sigma2 ===\n")
rho <- 0.15; sig2 <- 0.1; S <- cbind(runif(120, 0, 1), runif(120, 0, 1)); r <- rnorm(120)
alpha <- solve(matern32(as.matrix(dist(S)), rho) + diag(sig2, 120), r)
ok_w <- ok_v <- TRUE
for (s1 in seq(1.05, 1.5, by = 0.05)) {
  s <- c(s1, 0.5); h <- min(dist_to(S, s))
  w_hat <- sum(alpha * matern32(dist_to(S, s), rho)); kb <- kbar(matern32, h, rho)
  ok_w <- ok_w && abs(w_hat) <= sum(abs(alpha)) * kb + 1e-12
  ok_v <- ok_v && post_var_k(matern32, S, s, rho, sig2) >= 1 - 120 * kb^2 / sig2 - 1e-12
}
cat("field bound held at all ten block points:", ok_w, "| variance bound held:", ok_v, "\n\n")

# 12 ----------------------------------------------------------------------
cat("=== Corollary 6: a positive floor on the block variance for every n (Matern); the squared exponential for contrast ===\n")
s_b <- c(1.25, 0.5)                                   # depth 0.25 into the block
Smax <- cbind(runif(800, 0, 1), runif(800, 0, 1))      # nested designs: the first n rows
ns <- c(50, 100, 200, 400, 800)
vn_m <- sapply(ns, function(n) post_var_k(matern32, Smax[1:n, , drop = FALSE], s_b, 0.3, sig2))
vn_g <- sapply(ns, function(n) post_var_k(sqexp,    Smax[1:n, , drop = FALSE], s_b, 0.3, sig2))
# the floor: variance of w(s) given the noiseless field on ever finer grids of D minus B
floor_grid <- function(kern, m, jitter = 1e-9) {
  g <- seq(0, 1, length.out = m); G <- as.matrix(expand.grid(g, g))
  post_var_k(kern, G, s_b, 0.3, jitter)
}
fl_m <- sapply(c(21, 31, 41), function(m) floor_grid(matern32, m))
fl_g <- sapply(c(21, 31, 41), function(m) floor_grid(sqexp, m))
cat("Matern, posterior variance at s with n =", paste(ns, collapse = ", "), ":", paste(round(vn_m, 4), collapse = ", "), "\n")
cat("Matern, floor on grids of 21, 31, 41 points a side:", paste(round(fl_m, 4), collapse = ", "),
    "| every v_n above the finest floor:", all(vn_m >= min(fl_m)), "\n")
cat("squared exponential, posterior variance:", paste(signif(vn_g, 3), collapse = ", "), "\n")
cat("squared exponential, floor on the same grids:", paste(signif(fl_g, 3), collapse = ", "), "(falls as the grid is refined)\n\n")

# 13 ----------------------------------------------------------------------
cat("=== Lemma 6: decompositions agreeing on the training locations differ on the block by f - f' up to (||alpha||_1 + ||alpha'||_1) kbar(h_B) ===\n")
xcov <- function(S) sin(2 * S[, 1]) + 0.5 * S[, 2]           # a smooth covariate
S <- cbind(runif(150, 0, 1), runif(150, 0, 1)); Km <- matern32(as.matrix(dist(S)), 0.15)
alpha <- solve(Km + diag(sig2, 150), rnorm(150))               # field of the first decomposition
f  <- function(S) 0.5 + 1.0 * xcov(S); fp <- function(S) 0.5 + 1.6 * xcov(S)   # two environmental parts
eta_tr <- f(S) + drop(Km %*% alpha)                            # linear predictor on the training locations
alpha_p <- solve(Km, eta_tr - fp(S))                           # second field reproduces it exactly
cat("decompositions agree on the training locations:", isTRUE(all.equal(fp(S) + drop(Km %*% alpha_p), eta_tr)), "\n")
ok <- TRUE; worst <- 0
for (s1 in seq(1.05, 1.5, by = 0.05)) for (s2 in c(0.2, 0.5, 0.8)) {
  s <- c(s1, s2); ks <- matern32(dist_to(S, s), 0.15); h <- min(dist_to(S, s))
  gap <- (f(rbind(s)) + sum(alpha * ks)) - (fp(rbind(s)) + sum(alpha_p * ks)) - (f(rbind(s)) - fp(rbind(s)))
  bound <- (sum(abs(alpha)) + sum(abs(alpha_p))) * kbar(matern32, h, 0.15)
  ok <- ok && abs(gap) <= bound + 1e-10; worst <- max(worst, abs(gap))
}
cat("bound held at all 30 block points:", ok, "| largest departure from f - f':", signif(worst, 3), "\n\n")

# 14 ----------------------------------------------------------------------
cat("=== The fitted smooths: mgcv spatial bases used by the benchmark ===\n")
if (requireNamespace("mgcv", quietly = TRUE)) {
  d0 <- data.frame(lon = runif(200), lat = runif(200))
  gp <- mgcv::smoothCon(mgcv::s(lon, lat, bs = "gp", k = 20), data = d0)[[1]]
  tp <- mgcv::smoothCon(mgcv::s(lon, lat, k = 30), data = d0)[[1]]
  cat("Gaussian process smooth: correlation type", gp$gp.defn[1], "(3 = Matern, smoothness 3/2), range",
      round(gp$gp.defn[2], 4), "= largest distance between locations", round(max(dist(d0)), 4),
      "| unpenalised null space dimension", gp$null.space.dim, "(intercept and linear terms in the coordinates)\n")
  cat("thin plate regression spline: unpenalised null space dimension", tp$null.space.dim, "\n")
  cat("Matern 3/2 correlation at a distance equal to the range:", round(2 * exp(-1), 3), "\n")
} else cat("mgcv not installed; section skipped\n")
