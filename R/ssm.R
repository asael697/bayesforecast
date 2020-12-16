#' A  constructor for a Additive linear State space model.
#'
#' Constructor of the \code{ets("Z","Z","Z")} object for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @usage ssm(ts,trend = FALSE,damped = FALSE,seasonal = FALSE,xreg = NULL,
#'            period = 0,genT = FALSE,series.name = NULL)
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param trend a boolean value to specify a trend local level model. By default
#' is \code{FALSE}.
#' @param damped a boolean value to specify a damped trend local level model. By default
#' is \code{FALSE}. If \code{trend} option is \code{FALSE} then \code{damped} is set to
#' \code{FALSE} automatically.
#' @param seasonal a boolean value to specify a seasonal local level model. By default
#' is \code{FALSE}.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param period an integer specifying the periodicity of the time series by
#' default the value frequency(ts) is used.
#' @param genT a boolean value to specify for a generalized t-student SSM model.
#' @param series.name an optional string vector with the time series names.
#'
#' @details
#' By  default  the  \code{ssm()}  function generates a local level model (or a ets("A","N","N") or
#' exponential smoothing model from the \pkg{forecast} package). If \code{trend} is set \code{TRUE},
#' then  a  local  trend ssm model is defined (a equivalent ets("A","A","N") or Holt model from the
#' \pkg{forecast} package). For damped trend models set \code{damped} to \code{TRUE}. If \code{seasonal}
#' is  set  to  \code{TRUE} a seasonal local level model is defined (a equivalent ets("A","N","A") model
#' from  the  \pkg{forecast} package). For a Holt-Winters method (ets("A","A","A")) set \code{Trend} and
#' \code{seasonal} to \code{TRUE}.
#'
#' When \code{genT} option is \code{TRUE} a t-student innovations ssm model (see Ardia (2010)) is generated
#' see Fonseca, et. al (2019) for more details.
#'
#' The default priors used in a ssm( ) model are:
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
#' @seealso \code{\link{Sarima}} \code{\link{auto.arima}} \code{\link{set_prior}} \code{\link{garch}}
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
ssm = function(ts,trend = FALSE,damped = FALSE,seasonal = FALSE,xreg = NULL,
               period = 0,genT = FALSE,series.name = NULL){

    n = length(as.numeric(ts))
    y = as.numeric(ts)

    # series name
    if(is.null(series.name))
      sn = deparse(substitute(ts))
    else
      sn = as.character(series.name)

    # Check the damped trend
    is_dp = damped
    if(trend == FALSE)is_dp = FALSE

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
    m1$prior_seasonal1 = c(0,10,1,1)


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
#' @param object a  SSM object
#' @noRd
#'
is.ssm = function(object){
  y = FALSE
  if(is(object,"ssm")) y = TRUE
  return (y)
}
#' Extracts all the order coefficients in a list
#'
#' @param dat A SSM model
#' @noRd
#'
get_order_ssm = function(dat){
  return(list(level = q,trend = ifelse(dat$is_td,1,0),damped = ifelse(dat$is_dp,1,0),seasonal= ifelse(dat$is_ss,1,0)))
}
#' Max order  coefficients in a SSM model
#'
#' @param dat A SSM model
#' @noRd
#'
max_order_ssm = function(dat){

  return(max(c(1,dat$d1)))
}
