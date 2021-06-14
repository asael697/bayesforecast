#' Constructor a Multiplicative Seasonal ARIMA model.
#'
#' Constructor of the SARIMA model  for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package
#'
#' @usage Sarima(ts,order = c(1,0,0),seasonal = c(0,0,0),xreg = NULL,period = 0,series.name = NULL)
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param order A specification of the non-seasonal part of the ARIMA model: the
#' three components (p, d, q) are the AR order, the number of differences, and the
#' MA order.
#' @param seasonal A specification of the seasonal part of the ARIMA model,same as
#' order parameter:  the three components (p, d, q) are the seasonal AR order,
#' the degree of seasonal differences, and the seasonal MA order.
#' @param xreg	Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param period an integer specifying the periodicity of the time series by
#' default the value frequency(ts) is used.
#' @param series.name an optional string vector with the series names.
#'
#' @details
#' If \code{xreg} option is used, the model by default will cancel the
#' seasonal differences adjusted (D = 0). If a value \code{d} > 0 is used, all
#' the regressor variables in \code{xreg} will be difference as well.
#'
#' The default priors used in Sarima are:
#'
#' \itemize{
#'  \item{ar ~ normal(0,0.5)}
#'  \item{ma ~ normal(0,0.5)}
#'  \item{mu0 ~ t-student(0,2.5,6)}
#'  \item{sigma0 ~ t-student(0,1,7)}
#'  \item{sar ~ normal(0,0.5)}
#'  \item{sma ~ normal(0,0.5)}
#'  \item{breg ~ t-student(0,2.5,6)}
#' }
#'
#' For changing the default prior use the function \code{set_prior}
#'
#' @return The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @author  Asael Alonzo Matamoros
#'
#' @importFrom stats as.ts time frequency
#' @importFrom utils tail
#' @export
#'
#' @references
#' Box, G. E. P. and Jenkins, G.M. (1978). Time series analysis: Forecasting and
#' control. San Francisco: Holden-Day. \emph{Biometrika}, 60(2), 297-303.
#' \code{doi:10.1093/biomet/65.2.297}.
#'
#' Kennedy, P. (1992). Forecasting with dynamic regression models: Alan Pankratz, 1991.
#' \emph{International Journal of Forecasting}. 8(4), 647-648.
#' \code{url: https://EconPapers.repec.org/RePEc:eee:intfor:v:8:y:1992:i:4:p:647-648}.
#'
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}
#'
#' @seealso \code{\link{garch}} \code{\link{set_prior}}
#'
#' @examples
#' # Declare a multiplicative seasonal ARIMA model for the birth data.
#'
#' library(astsa)
#' model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#' model
#'
#' #Declare an Dynamic Harmonic Regression model for the birth data.
#' model = Sarima(birth,order = c(1,0,1),xreg = fourier(birth,K = 2))
#' model
#'
Sarima = function(ts,order = c(1,0,0),seasonal = c(0,0,0),xreg = NULL,
                  period = 0,series.name = NULL){

  n = length(as.numeric(ts))

  # series name
  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  m1 = list(n = n,dimension = 1,time = as.numeric(stats::time(ts)),
            p = no_negative_check(order[1]),
            d = no_negative_check(order[2]),
            q = no_negative_check(order[3]),
            yreal = stats::as.ts(ts),y = as.numeric(ts),series.name = sn)

  m1$prior_mu0 = c(0,2.5,6,4)
  m1$prior_sigma0 = c(0,1,7,4)
  m1$prior_ar  = matrix(rep(c(0,0.5,1,1),m1$p),ncol = 4,byrow = TRUE)
  m1$prior_ma  = matrix(rep(c(0,0.5,1,1),m1$q),ncol = 4,byrow = TRUE)
  m1$n1 = m1$n

  if(period == 0) m1$period = stats::frequency(ts)
  else m1$period = period

  m1$P = seasonal[1]
  m1$D = seasonal[2]
  m1$Q = seasonal[3]

  if(m1$period <= 1){
    m1$P=0
    m1$D=0
    m1$Q=0
  }

  m1$prior_sar = matrix(rep(c(0,0.5,1,1),m1$P),ncol = 4,byrow = TRUE)
  m1$prior_sma=  matrix(rep(c(0,0.5,1,1),m1$Q),ncol = 4,byrow = TRUE)

  # arima regression model
  if( !is.null(xreg) ){

    if(!is.matrix(xreg))
      stop("xreg has to be a matrix with row dimension as same as the length of the time serie")

    if(nrow(xreg) != n)
      stop("The length of xreg don't match with the length of the time serie")

    # seasonal adjustment
    if(m1$D > 0) {
      warning("seasonal difference is not allowed in dynamic regressions D  = 0 \n")
      m1$D = 0
    }

    m1$d1 = ncol(xreg)
    m1$xreg = xreg
    m1$reg = xreg
    m1$prior_breg  = matrix(rep(c(0,2.5,6,4),m1$d1),ncol = 4,byrow = TRUE)

    if(m1$d > 0){
      m1$xreg = diff(m1$xreg,differences = m1$d)
      m1$xlast = matrix(0,nrow = m1$d,ncol = m1$d1)
      m1$xlast[1,] = utils::tail(xreg,n=1)
      if(m1$d > 1)
        for(i in 2:m1$d) m1$xlast[i,] = utils::tail( diff(xreg,differences = i-1),n=1)
    }
  }
  else{
    m1$d1 = 0
    n2 = m1$n1-m1$d -(m1$period*m1$D)
    m1$xreg = matrix(rep(0,m1$d1*n2 ),ncol = m1$d1,nrow = n2)
    m1$prior_breg  = matrix(rep(c(0,2.5,6,4),m1$d1),ncol = 4,byrow = TRUE)
  }

  sc = dif(ts = as.numeric(ts),d = m1$d,D = m1$D,period = m1$period)
  m1$y = sc$y
  m1$n1 = length(m1$y)
  m1$init = sc$init
  m1$inits = sc$inits

  attr(m1,"class") = "Sarima"

  return(m1)
}
#' Checks if is a Sarima object
#'
#' @param object an arima object
#' @noRd
#'
is.Sarima = function(object){
  y = FALSE
  if(is(object,"Sarima")) y = TRUE
  return (y)
}
#' Extracts all the order coefficients in a list
#'
#' @param dat a Sarima model
#' @noRd
#'
get_order_arima= function(dat){
  return(list(p = dat$p,d =dat$d,q=dat$q,
              P=dat$P,D = dat$D,Q=dat$Q,
              d1 = dat$d1,
              period = dat$period))

}
#' Max order  coefficients in an Sarima model
#'
#' @param dat a Sarima model
#' @noRd
#'
max_order_arima= function(dat){
  return(max(c(dat$p,dat$q,dat$period*dat$P,dat$period*dat$Q)))
}
