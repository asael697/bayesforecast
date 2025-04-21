#' Naive and Random Walk models.
#'
#' Naive is the model constructor for a random walk  model applied to \code{y}.
#' This is equivalent to an ARIMA(0,1,0) model. \code{naive()} is simply a wrapper
#' to  maintain forecast package similitude. \code{seasonal} returns the model
#' constructor for a seasonal random walk equivalent to an ARIMA(0,0,0)(0,1,0)m
#' model where `m `is the seasonal period.
#'
#' The random walk with drift model is
#' \deqn{Y_t = mu_0 + Y_{t-1} + epsilon_t}{Y[t]= mu_0 +Y[t-1] + epsilon[t]}
#' where  \eqn{epsilon_t}{epsilon[t]} is a normal i.i.d. error.
#'
#' The seasonal naive model is
#' \deqn{Y_t = mu_0 + Y_{t-m} + epsilon_t}{Y[t]= mu_0 +Y[t-m] + epsilon[t]}
#' where  \eqn{epsilon_t}{epsilon[t]} is a normal i.i.d. error.
#'
#' @aliases naive
#'
#' @param ts  a numeric or ts object with the univariate time series.
#' @param seasonal a Boolean value for select a seasonal random walk instead.
#' @param m  an optional integer value for the seasonal period.
#'
#' @return The function returns a list with the data for running \code{stan()}
#' function of \pkg{rstan} package.
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{\link{Sarima}}
#'
#' @references
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}.
#'
#' Box, G. E. P. and Jenkins, G.M. (1978). Time series analysis: Forecasting and
#' control. San Francisco: Holden-Day. \emph{Biometrika}, 60(2), 297-303.
#' \code{doi:10.1093/biomet/65.2.297}.
#'
#' Kennedy, P. (1992). Forecasting with dynamic regression models: Alan Pankratz, 1991.
#' \emph{International Journal of Forecasting}. 8(4), 647-648.
#' \code{url: https://EconPapers.repec.org/RePEc:eee:intfor:v:8:y:1992:i:4:p:647-648}.
#'
#' @examples
#' # A seasonal Random-walk model.
#' model = naive(birth,seasonal = TRUE)
#' model
#'
#' @export
#'
naive = function(ts, seasonal = FALSE, m = 0){
  if(seasonal == FALSE)
    dat = Sarima(ts = ts,order = c(0,1,0), xreg = NULL, period = 0)
  else
    dat = Sarima(ts = ts, order = c(0,0,0), seasonal = c(0,1,0), xreg = NULL, period = m)

  attr(dat,"class") = "naive"
  return(dat)
}
#' Checks if is a `naive` object
#'
#' @param obj a `naive` object
#' @noRd
#'
is.naive = function(obj){
  y = FALSE
  if( is(obj,"naive")) y = TRUE
  return (y)
}
