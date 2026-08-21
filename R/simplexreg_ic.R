################################################################################
#               SIMPLEX REGRESSION - INFORMATION CRITERIA                      #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2026-08-10                                                             #
# Description: Information criteria and their finite-sample corrections for    #
#              simplex regression models                                       #
################################################################################

#' @title Information Criteria for Simplex Regression
#' @description Extracts the Akaike information criterion (AIC), the Schwarz
#' (Bayesian) information criterion (BIC), and the Hannan-Quinn criterion
#' (HQIC), along with their finite-sample corrections AICc, BICc, and HQICc follow
#' Bayer and Cribari-Neto (2017).
#'
#' @param ... One or more objects of class \code{"simplexregression"} to be
#' evaluated.
#'
#' @details
#' Let \eqn{\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})}
#' denote the maximized log-likelihood of the fitted model and \eqn{k} the
#' total number of estimated parameters. The AIC (Akaike, 1973) is
#' \deqn{AIC = -2\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})
#' + 2k.}
#' Sugiura (1978) proposed an unbiased estimator of the Kullback-Leibler
#' distance in linear regressions, later generalized by Hurvich and Tsai
#' (1989) to nonlinear regressions and autoregressive models. The resulting
#' finite-sample-corrected AIC is
#' \deqn{AICc = -2\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})
#' + \frac{2nk}{n-k-1}.}
#'
#' The Schwarz (1978) information criterion, also known as the Bayesian
#' information criterion (BIC), is a consistent criterion given by
#' \deqn{BIC = -2\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})
#' + k\log(n).}
#' McQuarrie (1999) derived a finite-sample-corrected version, which remains
#' consistent:
#' \deqn{BICc = -2\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})
#' + \frac{nk\log(n)}{n-k-1}.}
#'
#' The Hannan and Quinn (1979) criterion, proposed for autoregressive model
#' selection, is
#' \deqn{HQIC = -2\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})
#' + 2k\log(\log(n)).}
#' McQuarrie and Tsai (1998) proposed a finite-sample-corrected variant:
#' \deqn{HQICc = -2\ell(\widehat{\boldsymbol\beta}, \widehat{\boldsymbol\gamma})
#' + \frac{2nk\log(\log(n))}{n-k-1}.}
#'
#' Smaller values indicate better model fit.
#'
#' @return When a single model is provided, a named numeric vector with
#' components \code{AIC}, \code{AICc}, \code{BIC}, \code{BICc}, \code{HQIC},
#' and \code{HQICc}. When multiple models are provided, a data frame with one
#' row per model.
#'
#' @references
#' Akaike, H. (1973). Information theory and an extension of the maximum
#' likelihood principle. In \emph{Proceedings of the 2nd International
#' Symposium on Information Theory}, 267--281.
#'
#' Akaike, H. (1978). A Bayesian analysis of the minimum AIC procedure.
#' \emph{Annals of the Institute of Statistical Mathematics}, \bold{30}(1),
#' 9--14. \doi{10.1007/BF02480194}
#'
#' Bayer, F. M. and Cribari-Neto, F. (2017). Model selection criteria in beta
#' regression with varying dispersion. \emph{Communications in Statistics -
#' Simulation and Computation}, \bold{46}(1), 729--746.
#' \doi{10.1080/03610918.2014.977918}
#'
#' Hannan, E. J. and Quinn, B. G. (1979). The determination of the order of
#' an autoregression. \emph{Journal of the Royal Statistical Society, Series
#' B}, \bold{41}(2), 190--195. \doi{10.1111/j.2517-6161.1979.tb01072.x}
#'
#' Hurvich, C. M. and Tsai, C.-L. (1989). Regression and time series model
#' selection in small samples. \emph{Biometrika}, \bold{76}(2), 297--307.
#' \doi{10.1093/biomet/76.2.297}
#'
#' McQuarrie, A. D. R. (1999). A small-sample correction for the Schwarz SIC
#' model selection criterion. \emph{Statistics & Probability Letters},
#' \bold{44}(1), 79--86. \doi{10.1016/S0167-7152(98)00294-6}
#'
#' McQuarrie, A. D. R. and Tsai, C.-L. (1998). \emph{Regression and Time
#' Series Model Selection}. Singapore: World Scientific.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals
#' of Statistics}, \bold{6}(2), 461--464. \doi{10.1214/aos/1176344136}
#'
#' Sugiura, N. (1978). Further analysis of the data by Akaike's information
#' criterion and the finite corrections. \emph{Communications in Statistics -
#' Theory and Methods}, \bold{7}(1), 13--26.
#' \doi{10.1080/03610927808827599}
#'
#' @examples
#' data(ReadingSkills, package = "SimplexRegression")
#' fit1 <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                    data = ReadingSkills)
#'
#' # Single model
#' ic(fit1)
#'
#' # Comparing multiple models
#' fit2 <- simplexreg(accuracy ~ dyslexia * iq | dyslexia + iq + I(iq^2),
#'                    data = ReadingSkills, link.mu = "loglog")
#' ic(fit1, fit2)
#'
#' @seealso \code{\link{r2}}, \code{\link{press}}, \code{\link{ic.penalized}}
#' @export
ic <- function(...) {

  models <- list(...)
  M <- length(models)

  if (!all(sapply(models, function(x) inherits(x, "simplexregression")))) {
    stop("All arguments must be objects of class 'simplexregression'.")
  }

  # Get model names (supports named arguments)
  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- vapply(as.list(substitute(list(...)))[-1L], deparse, character(1))
  }

  aic <- aic_c <- bic <- bic_c <- hq <- hq_c <- numeric(M)

  for (i in seq_len(M)) {
    model <- models[[i]]

    n <- model$nobs
    p <- ncol(model$mu.x)
    q <- ncol(model$sigma2.x)
    parametric <- !is.na(model$coefficients$lambda)
    k <- if (parametric) p + q + 1 else p + q

    ll <- model$loglik

    # --- AIC (Akaike, 1973) and its finite-sample correction ---
    # (Sugiura, 1978; Hurvich and Tsai, 1989)
    aic[i]   <- -2 * ll + 2 * k
    aic_c[i] <- -2 * ll + (2 * n * k) / (n - k - 1)

    # --- SIC/BIC (Schwarz, 1978) and its finite-sample correction ---
    # (McQuarrie, 1999)
    bic[i]   <- -2 * ll + k * log(n)
    bic_c[i] <- -2 * ll + (n * k * log(n)) / (n - k - 1)

    # --- HQIC (Hannan and Quinn, 1979) and its finite-sample correction ---
    # (McQuarrie and Tsai, 1998)
    hq[i]   <- -2 * ll + 2 * k * log(log(n))
    hq_c[i] <- -2 * ll + (2 * n * k * log(log(n))) / (n - k - 1)
  }

  if (M == 1) {
    return(c(AIC = aic, AICc = aic_c,
             BIC = bic, BICc = bic_c,
             HQIC = hq, HQICc = hq_c))
  } else {
    data.frame(AIC = aic, AICc = aic_c,
               BIC = bic, BICc = bic_c,
               HQIC = hq, HQICc = hq_c,
               row.names = model_names)
  }
}
