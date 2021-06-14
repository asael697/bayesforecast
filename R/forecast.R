#' Forecasting varstan objects
#'
#' \code{forecast} is a generic function for forecasting from time series or
#' varstan models. The function invokes particular \emph{methods} which
#' depend on the class of the first argument.
#'
#' If \code{model=NULL},the function \code{\link{forecast.ts}} makes forecasts
#' using \code{\link{ets}} models (if the data are non-seasonal or the seasonal
#' period is 12 or less) or \code{\link{stlf}} (if the seasonal period is 13 or
#' more).
#'
#' If \code{model} is not \code{NULL}, \code{forecast.ts} will apply the
#' \code{model} to the \code{object} time series, and then generate forecasts
#' accordingly.
#'
#' @aliases forecast
#'
#' @param object a time series or varstan model for which forecasts are
#' required.
#' @param h Number of periods for forecasting.
#' @param probs A numerical vector \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \code{90\%} and \code{80\%} intervals (\code{level=c(0.8,0.9)}).
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param robust A boolean for obtain the robust estimation. The default
#' @param draws An integer indicating the number of draws to return. The default
#'    number of draws is 1000.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Further arguments passed to  \code{posterior_predict}.
#'
#' @return
#' An object of class "\code{forecast}".
#'
#' The function \code{summary} is used to obtain and print a summary of the
#' results, while the function \code{plot} produces a plot of the forecasts and
#' prediction intervals.
#'
#' The generic accessors functions \code{fitted.values} and \code{residuals}
#' extract various useful features of the value returned by
#' \code{forecast$model}.
#'
#' An object of class \code{"forecast"} is a list usually containing at least
#' the following elements: \item{model}{A list containing information about the
#' fitted model} \item{method}{The name of the forecasting method as a
#' character string} \item{mean}{Point forecasts as a time series}
#' \item{lower}{Lower limits for prediction intervals} \item{upper}{Upper
#' limits for prediction intervals} \item{level}{The confidence values
#' associated with the prediction intervals} \item{x}{The original time series
#' (either \code{object} itself or the time series used to create the model
#' stored as \code{object}).} \item{residuals}{Residuals from the fitted model.
#' For models with additive errors, the residuals will be x minus the fitted
#' values.} \item{fitted}{Fitted values (one-step forecasts)}
#'
#' @author Asael Alonzo Matamoros.
#'
#' @method forecast varstan
#' @importFrom forecast forecast
#' @importFrom stats frequency ts tsp tsp<-
#' @export forecast
#' @export
#'
#' @seealso The \code{"forecast"} methods of the forecast package.
#'
#' @examples
#' \donttest{
#'  fit = auto.sarima(ts = birth,iter = 500,chains = 1)
#'  fc = forecast(fit,h = 12)
#' }
#'
forecast.varstan = function(object, h = 10,probs = c(0.80, 0.90),xreg = NULL,robust = FALSE,
                            draws = 1000,seed = NULL,...){

  if (! is.varstan(object))
    stop("The current object is not a varstan class",call. = FALSE)

  # Preparing the data
  level = sort(probs);tsp.x =  stats::tsp(object$ts)

  if (!is.null(tsp.x))
    start.f = stats::tsp(object$ts)[2] + 1 / stats::frequency(object$ts)
  else
    start.f =length(object$x) + 1

  # Posterior predict
  pp = posterior_predict(object = object,h = h,xreg = xreg,robust = robust,
                         draws = draws,seed = seed,...)
  pp = as.matrix(pp);

  # Forecast credible intervals
  ppi = lower = upper = NULL
  ppm = apply(pp, 2, mean)

  for (i in level) ppi = cbind(ppi,posterior_interval(mat = pp,prob = i))

  # Reorganice
  for (i in 1:(2*length(level))){
    if(i%%2 == 1)
      lower = cbind(lower,ppi[,i])
    if(i%%2 == 0)
      upper = cbind(upper,ppi[,i])
  }

  colnames(lower) = colnames(upper)  = paste(100*level, "%", sep = "")
  row.names(lower) = row.names(upper) = NULL
  names(ppm) = NULL

  out = list(model = object, mean = stats::ts(ppm, frequency = stats::frequency(object$ts), start = start.f),
             level = level, x = object$ts)

  out$lower = stats::ts(lower)
  out$upper = stats::ts(upper)
  stats::tsp(out$lower) = stats::tsp(out$upper) = stats::tsp(out$mean)

  out$fitted  =  fitted.varstan(object)

  if (!is.null(object$series.name))
    out$series = object$series.name
  else
    out$series = deparse(substitute(object$ts))

  out$residuals = residuals.varstan(object)

  return(structure(out, class = "forecast"))
}
#'
#' @aliases forecast
#' @importFrom prophet prophet make_future_dataframe
#' @importFrom lubridate date_decimal
#' @importFrom stats is.ts time tsp
#' @export
#'
forecast.ts = function(object, h = 10,probs = c(0.8,0.95),...){

  if (! stats::is.ts(object))
    stop("The current object is not a time series class",call. = FALSE)

  # Preparing the data
  level = sort(probs);tsp.x =  stats::tsp(object)

  if (!is.null(tsp.x))
    start.f = stats::tsp(object)[2] + 1 / stats::frequency(object)
  else
    start.f =length(object) + 1

  # transform to a data frame object
  mts = as.numeric(stats::time(object))
  b = data.frame(ds =lubridate::date_decimal(mts),y = as.numeric(object) )

  ppm = lower = upper = NULL

  for (i in level){
    # call prophet for fitting the GAM
    mod = suppressMessages(prophet::prophet(df = b,interval.width = i))

    # Forecast
    future = prophet::make_future_dataframe(mod, periods = h,
                                            freq = findfreq(object),include_history = FALSE)
    fct = predict(mod,future)

    # Forecast credible intervals
    ppm = fct$yhat
    lower = cbind(lower,fct$yhat_lower)
    upper = cbind(upper,fct$yhat_upper)
  }
  # Rename
  colnames(lower) = colnames(upper)  = paste(100*level, "%", sep = "")
  row.names(lower) = row.names(upper) = NULL
  names(ppm) = NULL

  out = list(model = mod, mean = stats::ts(ppm, frequency = stats::frequency(object), start = start.f),
             level = level, x = object)

  out$lower = stats::ts(lower)
  out$upper = stats::ts(upper)
  stats::tsp(out$lower) = stats::tsp(out$upper) = stats::tsp(out$mean)

  out$series = deparse(substitute(object))

  # Fitted values and residuals
  out$fitted   = out$residual = object
  out$method = "prophet GAMs"

  return(structure(out, class = "forecast"))
}
#' Find the frequency for prophet methods
#' @noRd
#'
findfreq = function(object){
  if (! stats::is.ts(object))
    stop("The current object is not a time series class",call. = FALSE)

  m = stats::frequency(object)
  freq = m
  if(identical(m,365))freq = 'day'
  if(identical(m,52))freq = 'week'
  if(identical(m,12))freq = 'month'
  if(identical(m, 4))freq = 'quarter'
  if(identical(m, 1))freq = 'year'
  return(freq)
}
