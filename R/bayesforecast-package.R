#' Bayesian Time Series Modeling with \pkg{Stan}.
#'
#' @description
#' Fit univariate time series models using 'Stan' for full Bayesian inference.
#' A wide range of distributions and models are supported,  allowing users to
#' fit Seasonal ARIMA, ARIMAX, Dynamic Harmonic Regression,  GARCH, t-student
#' innovation GARCH models, asymmetric GARCH, Random Walks, and stochastic
#' volatility models. Prior specifications are flexible and explicitly encourage
#' users to apply prior distributions that actually reflect their beliefs. Model
#' fit can  easily  be assessed and compared with typical visualization methods,
#' information criteria such as loglik, AIC, BIC WAIC, Bayes factor and
#' leave-one-out cross-validation methods.
#'
#' @name bayesforecast-package
#' @aliases bayesforecast
#'
#' @useDynLib bayesforecast, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstan loo
#'
#' @references
#' Carpenter, B. and Gelman, A. and Hoffman, D. and  Lee, D. and Goodrich, B. and
#' Betancourt, M. and Brubaker, and Guo, L. and Riddell. 2017. Stan: A probabilistic
#' programming language. \emph{Journal of Statistical Software} 76(1).
#' \code{doi: 10.18637/jss.v076.i01}.
#'
#' Stan Development Team. (2018). Stan Modeling Language Users Guide and Reference Manual,
#' Version 2.18.0. \code{url: https://mc-stan.org}.
#'
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}.
#'
#' Tsay, R (2010). Analysis of Financial Time Series.
#' \emph{Wiley-Interscience}. 978-0470414354, second edition.
#'
#' Shumway, R.H. and Stoffer, D.S. (2010).Time Series Analysis and Its
#' Applications: With R Examples. \emph{Springer Texts in Statistics}.
#' isbn: 9781441978646. First edition.
#'
NULL
