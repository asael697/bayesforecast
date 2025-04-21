#' Constructor of an Stochastic volatility model object
#'
#' Constructor of the Stochastic Volatility model (SVM) for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#' \pkg{rstan} package.
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param arma Optionally, a specification of the  ARMA model,same
#' as order parameter: the two components `c(p, q)` are the AR, and
#' the  MA orders.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param series.name an optional string vector with the time series names.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return The function returns a list with the data for running \code{stan()}
#' function of \pkg{rstan} package.
#'
#' @export
#'
#' @references
#' Sangjoon,K. and Shephard, N. and Chib.S (1998). Stochastic Volatility: Likelihood
#' Inference and Comparison with ARCH Models. \emph{Review of Economic Studies}.
#' 65(1), 361-93. \code{url: https://www.jstor.org/stable/2566931}.
#'
#' Tsay, R (2010). Analysis of Financial Time Series.
#' \emph{Wiley-Interscience}. 978-0470414354, second edition.
#'
#' Shumway, R.H. and Stoffer, D.S. (2010).Time Series Analysis and Its
#' Applications: With R Examples. \emph{Springer Texts in Statistics}.
#' isbn: 9781441978646. First edition.
#'
#' @seealso \code{garch}, and  \code{et_prior}.
#'
#' @examples
#' # Declares a SVM model for the IPC data
#'
#' model = SVM(ipc, arma = c(1,1))
#' model
#'
SVM = function(ts, arma = c(0,0), xreg = NULL, series.name = NULL){

  m1 = garch(ts = ts,order = c(1,1,1),arma = arma,xreg = xreg,
             genT = FALSE,series.name = series.name)

  m1$prior_alpha = m1$prior_mgarch;m1$prior_garch = NULL;
  m1$prior_beta = m1$prior_arch; m1$prior_arch = NULL
  m1$h = 0;
  attr(m1,"class") = "SVM"
  return(m1)
}
#' Checks if is a SVM object.
#'
#' @param obj a SVM object.
#' @noRd
#'
is.SVM = function(obj){
  y = FALSE
  if( is(obj,"SVM")) y = TRUE
  return (y)
}
