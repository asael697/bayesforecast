#' Extract posterior sample of the point wise log-likelihood from a `varstan` object.
#'
#' Convenience function for extracting the point wise log-likelihood
#' matrix or array from a fitted Stan model.
#'
#' @param object a `varstan` object of the time series fitted model.
#' @param permuted a logical scalar indicating whether the draws after
#' the \code{warmup} period in each chain should be permuted and merged.
#' If FALSE, the original order is kept. For each \code{stanfit} object,
#' the permutation is fixed (i.e., extracting samples a second time
#' will give the same sequence of iterations).
#' @param ... additional values need in `log_lik` methods.
#'
#' @return
#' Usually, an S x N matrix containing the point wise log-likelihood
#' samples, where S is the number of samples and N is the number
#' of observations in the data. If \code{permuted} is \code{FALSE},
#' an S x N x R array is returned, where R is the number of fitted
#' chains.
#'
#' @aliases log_lik
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \code{doi:10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @importFrom rstantools log_lik
#' @importFrom loo  extract_log_lik
#' @method log_lik varstan
#' @export
#' @export log_lik
#'
#' @examples
#'
#' \donttest{
#'  model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(model,iter = 500,chains = 1)
#'
#'  log1 = log_lik(fit1)
#'  log1
#' }
#'
log_lik.varstan = function(object, permuted = TRUE, ...){

  if(!is.varstan(object))
    stop("The current object is not varstan class")

  log_lik = loo::extract_log_lik(object$stanfit, merge_chains = permuted)
  return(log_lik)
}
#' Extract posterior sample of the accumulated log-likelihood from a `varstan` object
#'
#' Convenience function for extracting the posterior sample of the accumulated
#' log-likelihood array from a fitted `varstan` object.
#'
#' @param object a `varstan` object of the time series fitted model.
#' @param permuted a logical scalar indicating whether the draws after
#' the `warmup`` period in each chain should be permuted and merged.
#' If `FALSE`, the original order is kept. For each `stanfit` object,
#' the permutation is fixed (i.e., extracting samples a second time
#' will give the same sequence of iterations).
#'
#' @return
#' A real value with the accumulated log likelihood.
#'
#' @export
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \code{doi:10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @examples
#' \donttest{
#'  model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(model,iter = 500,chains = 1)
#'
#'  log1 = loglik(fit1)
#'  log1
#' }
#'
loglik = function(object, permuted = TRUE){

  if(!is.varstan(object))
    stop("The current object is not varstan class")

  loglik = data.frame(extract_stan(object = object,
                                   pars = "loglik",
                                   permuted = permuted))

  if(permuted)
    loglik = as.numeric(loglik$loglik)
  else{
    colnames(loglik) = paste0("loglik.",1:ncol(loglik))
    loglik = as.matrix(loglik)
  }
  return(loglik)
}
#' Leave-one-out cross-validation
#'
#' The \code{loo} method for `varstan` objects. Computes approximate
#' leave-one-out cross-validation using Pareto smoothed importance
#' sampling (PSIS-LOO CV).
#'
#' @aliases loo
#'
#' @param x A `varstan` object
#' @param ... additional values need in loo methods
#'
#' @return an object from the loo class with the results of the Pareto-Smooth
#' Importance Sampling, leave one out cross validation for model selection.
#'
#' @seealso The \pkg{loo} package [vignettes](https://mc-stan.org/loo/articles/index.html)
#' for demonstrations. \code{psis()} for the underlying Pareto Smoothed Importance
#' Sampling (PSIS) procedure used in the LOO-CV approximation. \code{pareto-k-diagnostic}
#' for convenience functions for looking at diagnostics.\code{loo_compare()} for
#' model comparison.
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \code{doi:10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @importFrom rstan loo
#' @method loo varstan
#' @export loo
#' @export
#'
#' @examples
#' \donttest{
#' model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#' fit1 = varstan(model,iter = 500,chains = 1)
#'
#' loo1 = loo(fit1)
#' loo1
#' }
#'
loo.varstan = function(x, ...){
  if (!is.varstan(x))
    stop("The current object is not a varstan class")

  return(rstan::loo(x$stanfit) )
}
#' Widely Applicable Information Criterion (WAIC)
#'
#' Compute the widely applicable information criterion (WAIC)
#' based on the posterior likelihood using the \pkg{loo} package.
#' For more details see \code{\link[loo:waic]{waic}}.
#'
#' @param x A varstan object
#' @param ... additional values need in waic methods
#'
#' @aliases  waic
#'
#' @details See the \code{loo_compare} function of the \pkg{loo} package
#' for more details on model comparisons.
#'
#' @return An object of class \code{loo}. With the estimates of the
#' Watanabe-Akaike Information criteria.
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \code{doi:10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @importFrom loo waic
#' @method waic varstan
#' @export waic
#' @export
#'
#' @examples
#' \donttest{
#'  model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(model,iter = 500,chains = 1)
#'
#'  waic1 = waic(fit1)
#'  waic1
#' }
#'
waic.varstan = function(x,...){
  if (!is.varstan(x))
    stop("The current object is not a varstan class")

  return(loo::waic(log_lik.varstan(x)) )
}
#' Computes posterior sample of the point wise AIC method from a `varstan` object
#'
#' Convenience function for computing the point wise Akaike Information Criteria
#' method from a `varstan` object.
#'
#' @param x A `varstan` object of the time series fitted model.
#'
#' @return A numeric array  of size R, containing the posterior samples of the AICc
#' for a `varstan` object, where R is the number of iterations. If multiple chains are
#' fitted, then the array is of length M*R, where M is the number of chains
#'
#' @export
#'
#' @author Asael Alonzo Matamoros
#'
#' @examples
#' \donttest{
#'  model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(model,iter = 500,chains = 1)
#'
#'  aic1 = aic(fit1)
#'  mean(aic1)
#' }
#'
aic = function(x){
  if (!is.varstan(x))
    stop("The current object is not a varstan class")
  k = Total_order(x)
  aic = 2*k - 2*loglik(x)
  return(aic)
}
#' Computes posterior sample of the pointwise BIC method from a varstan object
#'
#' Convenience function for computing the pointwise Bayesian Information Criteria
#' method from a varstan object.
#'
#' @param x A varstan object of the time series fitted model.
#'
#' @return A numeric array  of size R, containing the posterior samples of the aic
#' for a varstan object, where R is the number of iterations. If multiple chains are
#' fitted, then the array is of length M*R, where M is the number of chains
#'
#' @export
#'
#' @author  Asael Alonzo Matamoros
#'
#' @examples
#'
#' \donttest{
#'  model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(model,iter = 500,chains = 1)
#'
#'  bic1 = bic(fit1)
#'  mean(bic1)
#' }
#'
bic = function(x){
  if (!is.varstan(x))
    stop("The current object is not a varstan class")
  k = Total_order(x)
  n = x$model$n
  bic = 2*k*log(n) - 2*loglik(x)
  return(bic)
}
#' Computes posterior sample of the point wise corrected AIC method from a `varstan` object
#'
#' Convenience function for computing the point wise corrected Akaike Information
#' Criterion method from a `varstan` object.
#'
#' @param x a `varstan` object of the time series fitted model.
#'
#' @return A numeric array  of size `R`, containing the posterior samples of the
#' AICc for a `varstan` object. where `R` is the number of iterations. If multiple
#' chains are fitted, then the array is of length `m*R`, where `m` is the number of
#' chains.
#'
#' @export
#'
#' @author  Asael Alonzo Matamoros
#'
#' @examples
#'
#' \donttest{
#'  model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(model,iter = 500,chains = 1)
#'
#'  aic1 = AICc(fit1)
#'  mean(aic1)
#' }
#'
AICc = function(x){
  if (!is.varstan(x))
    stop("The current object is not a varstan class")
  k = Total_order(x)
  n = x$model$n
  m = 2*(k^2 +k )/(n-k-1)
  aicc = 2*k - 2*loglik(x) +m
}
