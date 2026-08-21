library(testthat)
library(SimplexRegression)

# ==============================================================================
# test-simplexreg_ic.R
# Tests for: ic
# ==============================================================================

# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------
set.seed(42)
n   <- 100
x1  <- runif(n, 0, 1)
x2  <- runif(n, 0, 1)
mu  <- parametric_mean_link_inv(0.8 - 1.2*x1 - 1.5*x2, 0.25, "plogit2")
y   <- rsimplex(n, mu, 0.5)
dat <- data.frame(y = y, x1 = x1, x2 = x2)

fit_plogit2 <- simplexreg(y ~ x1 + x2 | 1, data = dat, link.mu = "plogit2")
fit_logit   <- simplexreg(y ~ x1 + x2 | 1, data = dat, link.mu = "logit")
fit_vardisp <- simplexreg(y ~ x1 + x2 | x1, data = dat, link.mu = "logit")

# ==============================================================================
# 1. Input validation
# ==============================================================================

test_that("ic rejects non-simplexregression objects", {
  expect_error(ic(list()),    "class 'simplexregression'")
  expect_error(ic(fit_logit, list()), "class 'simplexregression'")
})

# ==============================================================================
# 2. Single model — return structure
# ==============================================================================

test_that("ic returns a named numeric vector for a single model", {
  res <- ic(fit_logit)
  expect_true(is.numeric(res))
  expect_named(res, c("AIC", "AICc", "BIC", "BICc", "HQIC", "HQICc"))
})

test_that("ic returns a vector of length 6 for a single model", {
  expect_length(ic(fit_logit), 6L)
})

test_that("ic all components are finite for parametric link", {
  res <- ic(fit_plogit2)
  expect_true(all(is.finite(res)))
})

test_that("ic all components are finite for fixed link", {
  res <- ic(fit_logit)
  expect_true(all(is.finite(res)))
})

# ==============================================================================
# 3. Single model — consistency with fitted object
# ==============================================================================

test_that("ic AIC matches -2*loglik + 2k computed from model dimensions", {
  res <- ic(fit_logit)

  n <- fit_logit$nobs
  p <- ncol(fit_logit$mu.x)
  q <- ncol(fit_logit$sigma2.x)
  parametric <- !is.na(fit_logit$coefficients$lambda)
  k <- if (parametric) p + q + 1 else p + q

  ref <- -2 * fit_logit$loglik + 2 * k
  expect_equal(unname(res["AIC"]), ref, tolerance = 1e-10)
})

test_that("ic BIC matches -2*loglik + k*log(n) computed from model dimensions", {
  res <- ic(fit_logit)

  n <- fit_logit$nobs
  p <- ncol(fit_logit$mu.x)
  q <- ncol(fit_logit$sigma2.x)
  parametric <- !is.na(fit_logit$coefficients$lambda)
  k <- if (parametric) p + q + 1 else p + q

  ref <- -2 * fit_logit$loglik + k * log(n)
  expect_equal(unname(res["BIC"]), ref, tolerance = 1e-10)
})

test_that("ic HQIC matches -2*loglik + 2k*log(log(n)) computed from model dimensions", {
  res <- ic(fit_logit)

  n <- fit_logit$nobs
  p <- ncol(fit_logit$mu.x)
  q <- ncol(fit_logit$sigma2.x)
  parametric <- !is.na(fit_logit$coefficients$lambda)
  k <- if (parametric) p + q + 1 else p + q

  ref <- -2 * fit_logit$loglik + 2 * k * log(log(n))
  expect_equal(unname(res["HQIC"]), ref, tolerance = 1e-10)
})

test_that("ic parametric link has larger k (and thus different AIC) than fixed link", {
  res_par   <- ic(fit_plogit2)
  res_fixed <- ic(fit_logit)
  expect_false(isTRUE(all.equal(unname(res_par["AIC"]), unname(res_fixed["AIC"]))))
})

# ==============================================================================
# 4. Corrections — direction relative to uncorrected criteria
# ==============================================================================

test_that("ic AICc is larger than AIC (correction inflates the penalty)", {
  res <- ic(fit_logit)
  expect_true(unname(res["AICc"]) > unname(res["AIC"]))
})

test_that("ic BICc is larger than BIC (correction inflates the penalty)", {
  res <- ic(fit_logit)
  expect_true(unname(res["BICc"]) > unname(res["BIC"]))
})

test_that("ic HQICc is larger than HQIC (correction inflates the penalty)", {
  res <- ic(fit_logit)
  expect_true(unname(res["HQICc"]) > unname(res["HQIC"]))
})

test_that("ic corrections shrink toward the uncorrected criteria as k/n decreases", {
  # fit_vardisp has an extra sigma2 parameter (larger k) than fit_logit
  res_small_k <- ic(fit_logit)
  res_large_k <- ic(fit_vardisp)

  gap_small_k <- unname(res_small_k["AICc"] - res_small_k["AIC"])
  gap_large_k <- unname(res_large_k["AICc"] - res_large_k["AIC"])

  expect_true(gap_large_k > gap_small_k)
})

# ==============================================================================
# 5. Multiple models — return structure
# ==============================================================================

test_that("ic returns a data.frame for multiple models", {
  res <- ic(fit_plogit2, fit_logit)
  expect_s3_class(res, "data.frame")
})

test_that("ic data.frame has correct dimensions for two models", {
  res <- ic(fit_plogit2, fit_logit)
  expect_equal(nrow(res), 2L)
  expect_equal(ncol(res), 6L)
})

test_that("ic data.frame has correct column names", {
  res <- ic(fit_plogit2, fit_logit)
  expect_named(res, c("AIC", "AICc", "BIC", "BICc", "HQIC", "HQICc"))
})

test_that("ic data.frame row values match single-model calls", {
  res_multi  <- ic(fit_plogit2, fit_logit)
  res_single <- ic(fit_logit)
  expect_equal(as.numeric(res_multi[2, ]), as.numeric(res_single),
               tolerance = 1e-10)
})

test_that("ic data.frame all values are finite", {
  res <- ic(fit_plogit2, fit_logit)
  expect_true(all(is.finite(unlist(res))))
})

test_that("ic works for three models", {
  res <- ic(fit_plogit2, fit_logit, fit_vardisp)
  expect_equal(nrow(res), 3L)
})

test_that("ic data.frame row names use deparsed model names when unnamed", {
  res <- ic(fit_plogit2, fit_logit)
  expect_equal(rownames(res), c("fit_plogit2", "fit_logit"))
})

test_that("ic respects supplied names when models are named arguments", {
  res <- ic(m1 = fit_plogit2, m2 = fit_logit)
  expect_equal(rownames(res), c("m1", "m2"))
})
