#' A  constructor for local level state-space model.
#'
#' Constructor of the \code{ets("A","N","N")} object for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @usage LocalLevel(ts,xreg = NULL,genT = FALSE,series.name = NULL)
#'
#' @param ts a numeric or ts object with the univariate time series.

#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
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
#'  \item{sigma0 ~ t-student(0,1,7)}
#'  \item{level1 ~ normal(0,1)}
#'  \item{dfv ~ gamma(2,0.1)}
#'  \item{breg ~ t-student(0,2.5,6)}
#' }
#'
#' For changing the default prior use the function \code{set_prior()}.
#'
#' @return The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
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
#' mod1 = LocalLevel(ipc)
#'
LocalLevel = function(ts,xreg = NULL,genT = FALSE,series.name = NULL){

  dat = ssm(ts = ts,trend = FALSE,damped = FALSE,seasonal = FALSE,
            xreg = xreg,period = 0,genT = genT, series.name = series.name)

  attr(dat,"class") = "LocalLevel"
  return(dat)
}
#' Checks if is a Local level object
#'
#' @param object a  Local level object
#' @noRd
#'
is.LocalLevel = function(obj){
  y = FALSE
  if( is(obj,"LocalLevel")) y = TRUE
  return (y)
}
#' A constructor for a Holt trend state-space model.
#'
#' Constructor of the \code{ets("A","A","Z")} object for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @usage Holt(ts,damped = FALSE,xreg = NULL,genT = FALSE,series.name = NULL)
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param damped a boolean value to specify a damped trend local level model. By default
#' is \code{FALSE}. If \code{trend} option is \code{FALSE} then \code{damped} is set to
#' \code{FALSE} automatically.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param genT a boolean value to specify for a generalized t-student SSM model.
#' @param series.name an optional string vector with the time series names.
#'
#' @details
#' When \code{genT} option is \code{TRUE} a t-student innovations ssm model (see Ardia (2010)) is generated
#' see Fonseca, et. al (2019) for more details.
#'
#' The default priors used in a ssm( ) model are:
#'
#' \itemize{
#'  \item{level ~ normal(0,0.5)}
#'  \item{trend ~ normal(0,0.5)}
#'  \item{damped~ normal(0,0.5)}
#'  \item{sigma0 ~ t-student(0,1,7)}
#'  \item{level1 ~ normal(0,1)}
#'  \item{trend1 ~ normal(0,1)}
#'  \item{dfv ~ gamma(2,0.1)}
#'  \item{breg ~ t-student(0,2.5,6)}
#' }
#'
#' For changing the default prior use the function \code{set_prior()}.
#'
#' @return The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
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
#  #Declaring a Holt model for the ipc data.
#' mod1 = Holt(ipc)
#'
#' # Declaring a Holt damped trend model for the ipc data.
#' mod2 = Holt(ipc,damped = TRUE)
#'
Holt = function(ts,damped = FALSE,xreg = NULL,genT = FALSE,series.name = NULL){

  dat = ssm(ts = ts,trend = TRUE,damped = damped,seasonal = FALSE,
            xreg = xreg,period = 0,genT = genT, series.name = series.name)

  attr(dat,"class") = "Holt"
  return(dat)
}
#' Checks if is a Holt object
#'
#' @param object a Holt object
#' @noRd
#'
is.Holt = function(obj){
  y = FALSE
  if( is(obj,"Holt")) y = TRUE
  return (y)
}
#' A  constructor for a Holt-Winters state-space model.
#'
#' Constructor of the \code{ets("A","A","A")} object for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @usage Hw(ts,damped = FALSE,xreg = NULL,period = 0,genT = FALSE,series.name = NULL)
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param damped a boolean value to specify a damped trend local level model. By default
#' is \code{FALSE}. If \code{trend} option is \code{FALSE} then \code{damped} is set to
#' \code{FALSE} automatically.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param period an integer specifying the periodicity of the time series by
#' default the value frequency(ts) is used.
#' @param genT a boolean value to specify for a generalized t-student SSM model.
#' @param series.name an optional string vector with the time series names.
#'
#' @details
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
#' @return The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
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
#' @seealso \code{\link{Sarima}} \code{\link{auto.arima}} \code{\link{set_prior}}
#' \code{\link{garch}}
#'
#' @examples
#  #Declaring a Holt Winters model for the ipc data.
#' mod1 = Hw(ipc)
#'
#' # Declaring a Holt Winters damped trend model for the ipc data.
#' mod2 = Hw(ipc,damped = TRUE)
#'
#' # Declaring an additive Holt-Winters model for the birth data
#' mod3 = Hw(birth,damped = FALSE)
#'
Hw = function(ts,damped = FALSE,xreg = NULL,
               period = 0,genT = FALSE,series.name = NULL){

  dat = ssm(ts = ts,trend = TRUE,damped = damped,seasonal = TRUE,
            xreg = xreg,period = 0,genT = genT,
            series.name = series.name)

  attr(dat,"class") = "Hw"
  return(dat)
}
#' Checks if is a Holt-Winters object
#'
#' @param object a Holt-Winters object
#' @noRd
#'
is.Hw = function(obj){
  y = FALSE
  if( is(obj,"Hw")) y = TRUE
  return (y)
}
