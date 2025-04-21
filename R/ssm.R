#' A constructor for a Additive linear State space model.
#'
#' Constructor of the \code{ets("Z","Z","Z")} object for Bayesian estimation in \pkg{Stan}.
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param trend a bool value to specify a trend local level model. By default,
#' \code{trend = FALSE}.
#' @param damped a bool value to specify a damped trend local level model. By default,
#' \code{damped = FALSE}. If \code{trend = FALSE} then \code{damped = FALSE}
#' automatically.
#' @param seasonal a bool value to specify a seasonal local level model. By default
#' \code{seasonal = FALSE}.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param period an integer specifying the periodicity of the time series by
#' default the value `frequency(ts)` is used.
#' @param genT a bool value to specify for a generalized t-student SSM model.
#' @param series.name an optional string vector with the time series names.
#'
#' @details
#' By  default  the  \code{ssm()}  function generates a local level
#' `ets("A","N","N")`, or exponential smoothing model. If \code{trend = TRUE},
#' then  the model transforms into a  local  trend, `ets("A","A","N")` or Holt model
#' from the n\pkg{forecast} package. For damped trend models set \code{damped = TRUE}.
#' When \code{seasonal = TRUE}, the model becomes a seasonal local level or
#' `ets("A","N","A")`` model from  the  \pkg{forecast} package. Finally, a
#' Holt-Winters method or `ets("A","A","A")`,is whenever both \code{Trend} and
#' \code{seasonal} options are  \code{TRUE}.
#'
#' The \code{genT = TRUE} defines a t-student innovations SSM model. Check, Ardia (2010))
#' and Fonseca, et. al (2019) for more details.
#'
#' The default priors used in a `ssm( )` model are:
#'
#' \itemize{
#'  \item{level ~ normal(0,0.5)}
#'  \item{Trend ~ normal(0,0.5)}
#'  \item{damped~ normal(0,0.5)}
#'  \item{Seasonal ~ normal(0,0.5)}
#'  \item{sigma0 ~ t-student(0,1,7)}
#'  \item{level1 ~ normal(0,1)}
#'  \item{trend1 ~ normal(0,1)}
#'  \item{seasonal1 ~ normal(0,1)}
#'  \item{dfv ~ gamma(2,0.1)}
#'  \item{breg ~ t-student(0,2.5,6)}
#' }
#'
#' For changing the default prior use the function \code{set_prior()}.
#'
#' @return The function returns a list with the data for running \code{stan()}
#' function of \pkg{rstan} package.
#'
#' @author Asael Alonzo Matamoros.
#'
#' @export
#' @importFrom stats as.ts time frequency
#'
#' @references
#' Fonseca, T. and Cequeira, V. and Migon, H. and Torres, C. (2019). The effects of
#' degrees of freedom estimation in the Asymmetric GARCH model with Student-t
#' Innovations. \emph{arXiv} \code{doi: arXiv: 1910.01398}.
#'
#' @seealso \code{Sarima}, \code{auto.arima}, \code{set_prior}, and \code{garch}.
#'
#' @examples
#  #Declaring a local level model model for the ipc data.
#' mod1 = ssm(ipc)
#'
#' # Declaring a Holt model for the ipc data.
#' mod2 = ssm(ipc,trend = TRUE,damped = TRUE)
#'
#' # Declaring an additive Holt-Winters model for the birth data
#' mod3 = ssm(birth,trend = TRUE,damped = TRUE,seasonal = TRUE)
#'
ssm = function(ts, trend = FALSE, damped = FALSE, seasonal = FALSE, xreg = NULL,
               period = 0, genT = FALSE, series.name = NULL){

    n = length(as.numeric(ts))
    y = as.numeric(ts)

    # series name
    if(is.null(series.name))
      sn = deparse(substitute(ts))
    else
      sn = as.character(series.name)

    # Check the damped trend
    is_dp = damped
    if(trend == FALSE) is_dp = FALSE

    m1 = list(n = n,time = as.numeric(stats::time(ts)),
              is_td = trend, is_dp = is_dp,is_ss = seasonal,
              y = y,yreal = stats::as.ts(ts),series.name = sn)

    # Period
    if(period == 0) m1$period = stats::frequency(ts)
    else m1$period = period

    # Priors
    m1$prior_sigma0    = c(0,1,7,4)
    m1$prior_level     = c(0,0.5,1,1)
    m1$prior_trend     = c(0,0.5,1,1)
    m1$prior_damped    = c(0,0.5,1,1)
    m1$prior_seasonal  = c(0,0.5,1,1)
    m1$prior_level1    = c(mean(y[1:m1$period]) ,2.5,6,4)
    m1$prior_trend1    = c(0,2.5,6,4)
    m1$prior_seasonal1 = c(mean(y[1:m1$period]) ,2.5,6,4)


    # Generalized t distribution
    m1$genT = genT
    m1$prior_dfv = c(2,0.1,1,9)

    # Simple regression model
    if( !is.null(xreg) ){

      if(!is.matrix(xreg))
        stop("xreg has to be a matrix with row dimension as same as the length of the time serie")

      if(nrow(xreg) != n)
        stop("The length of xreg matrix don't match with the length of the time serie")

      m1$d1 = ncol(xreg)
      m1$xreg = xreg
    }
    else{
      m1$d1 = 0
      m1$xreg = matrix(rep(0,m1$d1*n),ncol = m1$d1,nrow = n)
    }
    m1$prior_breg  = matrix(rep(c(0,2.5,6,4),m1$d1),ncol = 4,byrow = TRUE)


    attr(m1,"class") = "ssm"
    return(m1)
}
#' Checks if is a SSM object
#'
#' @param object a  SSM object.
#' @noRd
#'
is.ssm = function(object){
  y = FALSE
  if(is(object,"ssm")) y = TRUE
  return (y)
}
#' Extracts all the order coefficients in a list
#'
#' @param dat A SSM model.
#' @noRd
#'
get_order_ssm = function(dat){
  return(list(level = 1,
              trend = ifelse(dat$is_td, 1, 0),
              damped = ifelse(dat$is_dp, 1, 0),
              seasonal= ifelse(dat$is_ss, 1, 0)
              )
         )
}
#' Max order  coefficients in a SSM model
#'
#' @param dat A SSM model
#' @noRd
#'
max_order_ssm = function(dat){

  return(max(c(1,dat$d1)))
}
