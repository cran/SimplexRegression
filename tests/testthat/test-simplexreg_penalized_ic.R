library(testthat)
library(SimplexRegression)

# ---------------------------------------------------------------------------
# Fixture: simulate data and fit models once for all tests
# ---------------------------------------------------------------------------
set.seed(2026)
n <- 100
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
z1 <- runif(n, 0, 1)
mu <- parametric_mean_link_inv(0.6 - 2*x1 - 1.5*x2, 0.5, "plogit1")
sigma2 <- dispersion_link_inv(-2 - 2.5*z1, "log")
y <- rsimplex(n, mu, sigma2)
data <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1)

fit_plogit1 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit1")
fit_plogit2 <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "plogit2")
fit_logit   <- simplexreg(y ~ x1 + x2 | z1, data = data, link.mu = "logit")

# ---------------------------------------------------------------------------

test_that("ic.penalized requires objects of class simplexregression", {
  expect_error(ic.penalized(list()),          "class 'simplexregression'")
  expect_error(ic.penalized(fit_plogit1, "x"), "class 'simplexregression'")
})

test_that("ic.penalized returns a data.frame with correct columns (kappa > 0)", {
  result <- ic.penalized(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("df", "AICp", "BICp", "HQICp"))
  expect_equal(nrow(result), 2)
})

test_that("ic.penalized returns standard column names when kappa = 0", {
  result <- ic.penalized(fit_plogit1, fit_plogit2, kappa = 0, verbose = FALSE)
  expect_named(result, c("df", "AIC", "BIC", "HQIC"))
})

test_that("ic.penalized returns finite numeric values", {
  result <- ic.penalized(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  expect_true(all(is.finite(result$AICp)))
  expect_true(all(is.finite(result$BICp)))
  expect_true(all(is.finite(result$HQICp)))
})

test_that("ic.penalized uses named arguments as row names", {
  result <- ic.penalized(m1 = fit_plogit1, m2 = fit_plogit2,
                         kappa = 0.1, verbose = FALSE)
  expect_equal(rownames(result), c("m1", "m2"))
})

test_that("ic.penalized works with a single model", {
  result <- ic.penalized(fit_plogit1, kappa = 0.1, verbose = FALSE)
  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(unlist(result))))
})

test_that("ic.penalized warns and sets kappa = 0 with non-parametric models", {
  expect_warning(
    ic.penalized(fit_plogit1, fit_logit, kappa = 0.1, verbose = FALSE),
    "kappa > 0 is only valid"
  )
})

test_that("ic.penalized returns invisibly when verbose = TRUE", {
  result <- withVisible(
    ic.penalized(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = TRUE)
  )
  expect_false(result$visible)
})

test_that("ic.penalized returns visibly when verbose = FALSE", {
  result <- withVisible(
    ic.penalized(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  )
  expect_true(result$visible)
})

test_that("ic.penalized df column matches number of estimated parameters", {
  result <- ic.penalized(fit_plogit1, fit_plogit2, kappa = 0.1, verbose = FALSE)
  expected_df <- length(fit_plogit1$coefficients$mean) +
                 length(fit_plogit1$coefficients$dispersion) + 1L
  expect_equal(result$df[1], expected_df)
})

test_that("ic.penalized silently sets kappa = 0 ...", {
  result <- expect_no_warning(ic.penalized(fit_logit, verbose = TRUE))
  expect_named(result, c("df", "AIC", "BIC", "HQIC"))
  expect_true(all(is.finite(unlist(result))))
})
