#' Automatic estimate of a Seasonal ARIMA model
#'
#' Returns the best seasonal ARIMA model using a \code{bic} value, this
#' function the\code{auto.arima} function of the \pkg{forecast} package
#' to select the seasonal ARIMA model and estimates the model using a
#' HMC sampler.
#'
#' @usage auto.sarima(ts,xreg = NULL,chains=1,iter=4000,warmup=floor(iter/2),
#'                 adapt.delta = 0.9,tree.depth =10,stepwise = TRUE, series.name = NULL,
#'                 prior_mu0 = NULL,prior_sigma0 = NULL,prior_ar = NULL, prior_ma = NULL,
#'                 prior_sar = NULL,prior_sma = NULL, prior_breg = NULL,...)
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param chains An integer of the number of Markov Chains chains to be run,
#' by default 4 chains are run.
#' @param iter An integer of total iterations per chain including the warm-up,
#' by default  the number of iterations are 2000.
#' @param warmup  A positive integer specifying number of warm-up (aka burn-in)
#'   iterations. This also specifies the number of iterations used for step-size
#'   adaptation, so warm-up samples should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param adapt.delta An optional real value between 0 and 1, the thin of the jumps
#' in a HMC method. By default is 0.9.
#' @param  tree.depth An integer of the maximum  depth of the trees  evaluated
#' during each iteration. By default is 10.
#' @param stepwise	If TRUE, will do stepwise selection (faster). Otherwise, it searches
#' over all models. Non-stepwise selection can be very slow, especially for seasonal models.
#' @param series.name an optional string vector with the series names.
#' @param prior_mu0 The prior distribution for the location parameter in an ARIMA model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param prior_sigma0 The prior distribution for the scale parameter in an ARIMA model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param prior_ar The prior distribution for the auto-regressive parameters in an ARIMA model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_ma The prior distribution for the moving average parameters in an ARIMA model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_sar The prior distribution for the seasonal auto-regressive parameters in a
#' SARIMA model. By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_sma The prior distribution for the seasonal moving average parameters in a
#' SARIMA model. By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_breg The prior distribution for the regression coefficient parameters in a
#' ARIMAX model. By default the value is set \code{NULL}, then the default student(7,0,1) priors
#' are used.
#' @param ... Further arguments passed to  \code{auto.arima} function.
#'
#' @details
#' Automatic ARIMA model fitting implemented by Rob Hyndman, this function finds the best
#' Seasonal ARIMA model using \code{bic}, and then proceeds to fit the model using
#' \code{varstan} function and the default priors of a \code{Sarima} model constructor.
#'
#' This function provides an initial model fit for beginning the Bayesian analysis
#' of the univariate time series. For better fit and model selection try different
#' models and other model selection criteria such as \code{loo} or \code{bayes_factor}.
#'
#' The default arguments are designed for rapid estimation of models for many time series.
#' If you are analyzing just one time series, and can afford to take some more time, it is
#' recommended that you set \code{stepwise}=\code{FALSE} and reduce the number of iterations
#' per chain (\code{iter}).
#'
#' For more information look at \code{auto.arima()} function of forecast package.
#'
#' @author Asael Alonzo Matamoros
#'
#' @importFrom forecast auto.arima
#' @export
#'
#' @return a varstan model
#'
#' @seealso \code{Sarima} \code{varstan}.
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
#' \dontrun{
#'  # Automatic Sarima model for the birth data
#'  auto.sarima(birth)
#'
#'  # Dynamic Harmonic regression
#'  auto.sarima(birth,xreg = fourier(birth,K= 6))
#'}
#'
auto.sarima = function(ts,xreg= NULL,chains = 4,iter = 4000,warmup = floor(iter/2),
                       adapt.delta = 0.9,tree.depth = 10,stepwise = TRUE,
                       series.name = NULL,prior_mu0 = NULL,prior_sigma0 = NULL,
                       prior_ar = NULL, prior_ma = NULL, prior_sar = NULL,
                       prior_sma = NULL, prior_breg = NULL,...){


  # Selecting the "best" Arima model using forecast package
  y  = auto.arima(y = ts,xreg = xreg,ic = "bic",stepwise = stepwise,...)
  arma = y$arma
  ord = c(arma[1],arma[6],arma[2])
  season = c(arma[3],arma[7],arma[4])

  # series name
  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  dat = Sarima(ts,order = ord,
               seasonal = season,
               period = arma[5],
               xreg = xreg,series.name = sn)

  # Priors selection
  if(!is.null(prior_mu0)) dat = set_prior(model = dat,par = "mu0",dist = prior_mu0)
  if(!is.null(prior_sigma0)) dat = set_prior(model = dat,par = "sigma0",dist = prior_sigma0)
  if(!is.null(prior_ar)) dat = set_prior(model = dat,par = "ar",dist = prior_ar)
  if(!is.null(prior_ma)) dat = set_prior(model = dat,par = "ma",dist = prior_ma)
  if(!is.null(prior_sar)) dat = set_prior(model = dat,par = "sar",dist = prior_sar)
  if(!is.null(prior_sma)) dat = set_prior(model = dat,par = "sma",dist = prior_sma)
  if(!is.null(prior_breg)) dat = set_prior(model = dat,par = "breg",dist = prior_breg)

  # Fitting the Sarima model.
  sf1 = varstan(model = dat,
                chains = chains,
                iter = iter,
                warmup = warmup,
                adapt.delta = adapt.delta,
                tree.depth = tree.depth)
  return(sf1)
}
#' Fourier terms for modeling seasonality.
#'
#' \code{fourier} returns a matrix containing terms from a Fourier series, up
#' to order \code{K}, suitable for use in \code{\link{Sarima}} or
#' \code{\link{auto.sarima}}.
#'
#' The period of the Fourier terms is determined from the time series
#' characteristics of \code{x}. When \code{h} is missing, the length of
#' \code{x} also determines the number of rows for the matrix returned by
#' \code{fourier}. Otherwise, the value of \code{h} determines the number of
#' rows for the matrix returned by \code{fourier}, typically used for
#' forecasting. The values within \code{x} are not used.
#'
#' Typical use would omit \code{h} when generating Fourier terms fitting a model
#' and include \code{h} when generating Fourier terms for forecasting.
#'
#' When \code{x} is a \code{ts} object, the value of \code{K} should be an
#' integer and specifies the number of sine and cosine terms to return. Thus,
#' the matrix returned has \code{2*K} columns.
#'
#' When \code{x} is a \code{msts} object, then \code{K} should be a vector of
#' integers specifying the number of sine and cosine terms for each of the
#' seasonal periods. Then the matrix returned will have \code{2*sum(K)}
#' columns.
#'
#' @param x Seasonal time series: a \code{ts} or a \code{msts} object
#' @param K Maximum order(s) of Fourier terms
#' @param h Number of periods ahead to forecast (optional)
#'
#' @return Numerical matrix.
#'
#' @author Rob J Hyndman
#'
#' @seealso \code{\link{seasonaldummy}}
#'
#' @keywords forecast
#'
#' @examples
#' \dontrun{
#'  library(astsa)
#'  # Dynaimc Harmonic regression
#'  sf1 = auto.sarima(birth,xreg = fourier(birth,K= 6))
#' }
#'
#' @importFrom forecast fourier
#' @export
#'
fourier <- function(x, K, h = NULL) {
  return(forecast::fourier(x = x,K = K,h = h))
}
