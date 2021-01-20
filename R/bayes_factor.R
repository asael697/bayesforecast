#' Log Marginal Likelihood via Bridge Sampling.
#'
#' Computes log marginal likelihood via bridge sampling,
#' which can be used in the computation of Bayes factors
#' and posterior model probabilities.
#'
#' The \code{varstan} class is just a thin wrapper that
#' contains the \code{stanfit} objects.
#'
#' @aliases bridge_sampler
#'
#' @param samples A \code{varstan} object.
#' @param ... Additional arguments passed to
#'   \code{\link[bridgesampling:bridge_sampler]{bridge_sampler.stanfit}}.
#'
#' @details
#' Computing the marginal likelihood  via the bridgesampler package
#' for stanfit objects.
#'
#' The computation of marginal likelihoods based on bridge sampling requires
#' a lot more posterior samples than usual. A good conservative rule of thump
#' is perhaps 10-fold more samples (read: the default of 4000 samples may not
#' be enough in many cases). If not enough posterior samples are provided, the
#' bridge sampling algorithm tends to be unstable leading to considerably different
#' results each time it is run. We thus recommend running \code{bridge_sampler}
#' multiple times to check the stability of the results.
#'
#' For  more details check the \pkg{bridgesampling} package.
#'
#' @return the model's marginals likelihood from the \code{bridge_sampler} package.
#'
#' @method bridge_sampler varstan
#' @importFrom bridgesampling bridge_sampler
#' @export bridge_sampler
#' @export
#'
#' @examples
#' \donttest{
#' library(astsa)
#' # Fitting a seasonal ARIMA model
#' mod1 = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#' fit1 = varstan(mod1,chains = 1)
#'
#' fit1
#' bridge_sampler(fit1)
#'
#' # Fitting a Dynamic harmonic regression
#' mod2  = Sarima(birth,order = c(0,1,2),xreg = fourier(birth,K=6))
#' fit2 = varstan(mod2,chains = 1)
#'
#' fit2
#' bridge_sampler(fit2)
#' }
#'
bridge_sampler.varstan <- function(samples, ...) {
  if(!is.varstan(samples))
    stop("The current object is not a varstan class")

  out = try(bridge_sampler(samples$stanfit, ...))
  return(out)
}
#' Bayes Factors from Marginal Likelihoods.
#'
#' Compute Bayes factors from marginal likelihoods.
#'
#' @aliases bayes_factor
#'
#' @param x1 A \code{varstan} object
#' @param x2 Another \code{varstan} object based on the same data.
#' @param log A boolean parameter for report the Bayes_factor in log scale.
#' The default value is FALSE.
#' @param ... Additional arguments passed to \code{bayes_factor}.
#'
#' @details
#' The computation of marginal likelihoods based on bridge sampling requires
#' a lot more posterior samples than usual. A good conservative rule of thump
#' is perhaps 10-fold more samples (read: the default of 4000 samples may not
#' be enough in many cases). If not enough posterior samples are provided, the
#' bridge sampling algorithm tends to be unstable leading to considerably different
#' results each time it is run. We thus recommend running \code{bridge_sampler}
#' multiple times to check the stability of the results.
#'
#' For  more details check the \pkg{bridgesampling} package.
#'
#' @return The bayes factors of two models.
#'
#' @method bayes_factor varstan
#' @importFrom bridgesampling bayes_factor
#' @export bayes_factor
#' @export
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  # Fitting a seasonal arima model
#'  mod1 = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(mod1,chains = 1)
#'
#'  # Fitting a Dynamic harmonic regression
#'  mod2  = Sarima(birth,order = c(0,1,2),xreg = fourier(birth,K=6))
#'  fit2 = varstan(mod2,chains = 1)
#'
#'  # compute the Bayes factor
#'  bayes_factor(fit1, fit2)
#' }
#'
bayes_factor.varstan <- function(x1, x2, log = FALSE, ...) {
  bridge1 = bridge_sampler(x1, ...)
  bridge2 = bridge_sampler(x2, ...)
  out = bayes_factor(bridge1, bridge2, log = log)
  attr(out, "model_names") = c("model1", "model2")
  return(out)
}
