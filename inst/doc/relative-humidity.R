## ----include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  fig.align  = "center",
  message    = FALSE,
  warning    = FALSE
)
old_opts <- options(width = 70, prompt = "R> ", continue = "+  ", digits = 5)

## ----setup----------------------------------------------------------
library(SimplexRegression)
data(RelativeHumidity, package = "SimplexRegression")
head(RelativeHumidity, 5)

## ----data-prep------------------------------------------------------
rh       <- RelativeHumidity
rh$hs    <- sin(2 * pi * seq_len(nrow(rh)) / 12)
rh$hc    <- cos(2 * pi * seq_len(nrow(rh)) / 12)
rh$dummy <- as.integer(as.integer(format(rh$Date, "%m")) %in% 10:12)

## ----summary-rh-----------------------------------------------------
summary(rh$RH)
cat(sprintf(
  "Std. dev.: %.4f  |  Skewness: %.4f\n",
  sd(rh$RH),
  mean(((rh$RH - mean(rh$RH)) / sd(rh$RH))^3)
))

## ----formula--------------------------------------------------------
formula <- RH ~ Ins2 + MT + WS + hs + hc + dummy + I(dummy * WS) | Pre2

## ----models-plogit--------------------------------------------------
fit_p1 <- simplexreg(formula, data = rh, link.mu = "plogit1")
fit_p2 <- simplexreg(formula, data = rh, link.mu = "plogit2")

## ----penalized-ss-param---------------------------------------------
penalized.ss(fit_p1, fit_p2, kappa = 0.1)

## ----penalized-ic---------------------------------------------------
penalized.ic(fit_p1, fit_p2, kappa = 0.1)

## ----models-fixed---------------------------------------------------
fit_loglog  <- simplexreg(formula, data = rh, link.mu = "loglog")
fit_logit   <- simplexreg(formula, data = rh, link.mu = "logit")
fit_probit  <- simplexreg(formula, data = rh, link.mu = "probit")
fit_cauchit <- simplexreg(formula, data = rh, link.mu = "cauchit")
fit_cloglog <- simplexreg(formula, data = rh, link.mu = "cloglog")

## ----penalized-ss-all-----------------------------------------------
penalized.ss(
  fit_loglog, fit_logit, fit_probit,
  fit_cauchit, fit_cloglog, fit_p1,
  kappa = 0
)

## ----summary-fit----------------------------------------------------
summary(fit_loglog)

## ----aic-bic--------------------------------------------------------
AIC(fit_loglog, fit_logit, fit_probit, fit_cauchit, fit_cloglog, fit_p1)
BIC(fit_loglog, fit_logit, fit_probit, fit_cauchit, fit_cloglog, fit_p1)
HQIC(fit_loglog, fit_logit, fit_probit, fit_cauchit, fit_cloglog, fit_p1)

## ----coef-vcov------------------------------------------------------
coef(fit_loglog)                          # full coefficient vector
coef(fit_loglog, model = "mean")          # mean submodel only
coef(fit_loglog, model = "dispersion")    # dispersion submodel only
round(vcov(fit_loglog, model = "mean"), 6)  # vcov of mean submodel

## ----loglik---------------------------------------------------------
logLik(fit_loglog)

## ----lrtest---------------------------------------------------------
fit_loglog_null <- update(fit_loglog, . ~ . | 1)
lmtest::lrtest(fit_loglog, fit_loglog_null)

## ----scoretest------------------------------------------------------
fit_logit_h0 <- simplexreg(formula, data = rh, link.mu = "logit")
scoretest(fit_logit_h0, link.mu = "plogit1")
scoretest(fit_logit_h0, link.mu = "plogit2")

## ----resettest------------------------------------------------------
resettest(fit_loglog)                      # both submodels augmented
resettest(fit_loglog, dispersion = FALSE)  # mean submodel only

## ----fitted-resid---------------------------------------------------
head(fitted(fit_loglog))
head(residuals(fit_loglog, type = "quantile"))  # approx. N(0,1)
head(residuals(fit_loglog, type = "pearson"))
head(residuals(fit_loglog, type = "weighted"))  # for halfnormal.plot

## ----press----------------------------------------------------------
press(fit_loglog)                            # single model
press(fit_loglog, fit_logit, fit_probit)     # comparing models

## ----predict--------------------------------------------------------
head(predict(fit_loglog, type = "response"))        # fitted means
head(predict(fit_loglog, type = "link")$mean)       # mean linear predictor
head(predict(fit_loglog, type = "link")$dispersion) # dispersion predictor
head(predict(fit_loglog, type = "dispersion"))      # fitted sigma^2

# Out-of-sample prediction
new_obs <- rh[1:3, ]
predict(fit_loglog, newdata = new_obs, type = "response")

## ----simulate-------------------------------------------------------
set.seed(2026)
sims <- simulate(fit_loglog, nsim = 3)
head(sims)

## ----influence------------------------------------------------------
hii  <- hatvalues(fit_loglog)
cook <- cooks.distance(fit_loglog, type = "pearson")
cat(sprintf("Leverages  — max: %.4f  mean: %.4f\n", max(hii),  mean(hii)))
cat(sprintf("Cook's D   — max: %.4f\n", max(cook)))

## ----gleverage------------------------------------------------------
gl <- gleverage(fit_loglog)
cat(sprintf("Generalized leverage — max: %.4f  mean: %.4f\n",
            max(gl), mean(gl)))

## ----plots-1-5, fig.height = 8, fig.cap = "Diagnostic plots (1–6) for the fitted simplex regression model with log-log link."----
oldpar <- par(mfrow = c(3, 2))
plot(fit_loglog, which = 1:5, reset.par = FALSE)
par(oldpar)

## ----plot-cook, fig.height = 4, fig.cap = "Cook's distances. Observations exceeding the threshold of 0.15 are labeled."----
plot(fit_loglog, which = 6, threshold = 0.15, label.pos = 4)

## ----plot-glev, fig.height = 4, fig.cap = "Generalized leverage values. Observations exceeding 0.08 are labeled."----
plot(fit_loglog, which = 7, threshold = 0.08, label.pos = 3)

## ----local-influence-cw, fig.height = 4.5, fig.cap = "Total local influence $C_i$ under case-weight perturbation for all parameters."----
local.influence(
  fit_loglog,
  scheme    = "case.weight",
  parameter = "theta",
  type      = "Ci",
  plot      = TRUE,
  threshold = 0.5,
  label.pos = c(3, 4, 3, 2, 2)
)

## ----local-influence-resp, fig.height = 4.5, fig.cap = "Total local influence $C_i$ under response perturbation for all parameters."----
local.influence(
  fit_loglog,
  scheme    = "response",
  parameter = "theta",
  type      = "Ci",
  plot      = TRUE,
  threshold = 0.4,
  label.pos = 2
)

## ----halfnormal, fig.height = 5, fig.cap = "Half-normal plot of absolute weighted residuals with 95% simulated envelope (100 replications)."----
halfnormal.plot(fit_loglog, nsim = 19, type = "weighted", seed = 2026)

## ----timeseries, fig.height = 4, fig.cap = "Observed (solid black) and fitted (dashed red) monthly relative humidity in Brasília, January 2000 to December 2025."----
plot(rh$Date, rh$RH,
     type = "l", col = "black", lwd = 1.2,
     xlab = "Date", ylab = "Relative humidity",
     main = "Observed vs Fitted RH — Brasília (2000–2025)")
lines(rh$Date, fitted(fit_loglog), col = "red", lwd = 1.5, lty = 2)
legend("bottomleft",
       legend = c("Observed", "Fitted"),
       col    = c("black", "red"),
       lty    = c(1, 2), lwd = c(1.2, 1.5),
       bty    = "n", cex = 0.85)

## ----restore-options, include = FALSE-----------------------------------------
options(old_opts)

## ----session------------------------------------------------------------------
sessionInfo()

