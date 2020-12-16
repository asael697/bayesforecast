#' Fitting a Multiplicative Seasonal ARIMA model.
#'
#' Fitting a SARIMA model  in \pkg{Stan}.
#'
#' The function returns a \code{varstan} object with the fitted model.
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param order A specification of the non-seasonal part of the ARIMA model: the
#' three components (p, d, q) are the AR order, the number of differences, and the
#' MA order.
#' @param seasonal A specification of the seasonal part of the ARIMA model,same as
#' order parameter:  the three components (p, d, q) are the seasonal AR order,
#' the degree of seasonal differences, and the seasonal MA order.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param period an integer specifying the periodicity of the time series by
#' default the value frequency(ts) is used.
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
#' @param stepwise If TRUE, will do stepwise selection (faster). Otherwise, it searches
#' over all models. Non-stepwise selection can be very slow, especially for seasonal models.
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
#' @param series.name an optional string vector with the series names.
#' @param ... Further arguments passed to  \code{varstan} function.
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
#' \dontrun{
#'  library(astsa)
#'  # Declare a multiplicative seasonal ARIMA model for the birth data.
#'  sf1 = stan_Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'
#'
#'  #Declare an Dynamic Harmonic Regression model for the birth data.
#'  sf2 = stan_Sarima(birth,order = c(1,0,1),xreg = fourier(birth,K = 2))
#' }
#'
stan_sarima = function(ts,order = c(1,0,0),seasonal = c(0,0,0),xreg = NULL,
                       period = 0,chains = 4,iter = 4000,warmup = floor(iter/2),
                       adapt.delta = 0.9,tree.depth = 10,stepwise = TRUE,
                       prior_mu0 = NULL,prior_sigma0 = NULL,prior_ar = NULL,
                       prior_ma = NULL, prior_sar = NULL,prior_sma = NULL,
                       prior_breg = NULL,series.name = NULL,...){

  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  dat = Sarima(ts = ts,order = order,seasonal = seasonal,xreg = xreg,period = period,series.name = sn)

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
                tree.depth = tree.depth,...)
  return(sf1)
}
#' Fitting for a GARCH(s,k,h) model.
#'
#' Fitting a \code{GARCH(s,k,h)} model in \pkg{Stan}.
#'
#' The function returns a \code{varstan} object with the fitted model.
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param order A specification of the garch  model: the three components (s, k, h)
#' are the arch order, the garch order, and the mgarch order.
#' @param arma A specification of the  ARMA model,same as order parameter:  the two
#' components (p, q) are the AR order,and the  MA order.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param genT a boolean value to specify for a generalized t-student garch model.
#' @param asym a string value for the asymmetric function for an asymmetric GARCH process. By default
#' the value \code{"none"} for standard GARCH process. If \code{"logit"} a logistic function
#' is used for asymmetry, and if \code{"exp"} an exponential function is used.
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
#' @param stepwise If TRUE, will do stepwise selection (faster). Otherwise, it searches
#' over all models. Non-stepwise selection can be very slow, especially for seasonal models.
#' @param prior_mu0 The prior distribution for the location parameter in an ARMA model. By default
#' the value is set \code{NULL}, then the default normal(0,1) prior is used.
#' @param prior_sigma0 The prior distribution for the scale parameter in an ARMA model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param prior_ar The prior distribution for the auto-regressive parameters in an ARMA model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_ma The prior distribution for the moving average parameters in an ARMA model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_arch The prior distribution for the arch parameters in a GARCH model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_garch The prior distribution for the GARCH parameters in a GARCH model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_mgarch The prior distribution for the mean GARCH parameters in a
#' GARCH model. By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_breg The prior distribution for the regression coefficient parameters in a
#' ARIMAX model. By default the value is set \code{NULL}, then the default student(7,0,1) priors
#' are used.
#' @param prior_df The prior distribution for the degree freedom parameters in a t-student innovations
#' GARCH model. By default the value is set \code{NULL}, then the default gamma(2,0.1) priors
#' are used.
#' @param prior_gamma The prior distribution for the asymmetric parameters in am Asymmetric
#' GARCH model. By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param series.name an optional string vector with the series names.
#' @param ... Further arguments passed to  \code{varstan} function.
#'
#' @details
#' By default the \code{garch()} function generates a GARCH(1,1) model, when
#' \code{genT} option is \code{TRUE} a t-student innovations GARCH model
#' (see Ardia (2010)) is generated, and for Asymmetric GARCH models use the
#' option \code{asym} for specify the asymmetric function, see Fonseca,
#' et. al (2019) for more details.
#'
#' The default priors used in a GARCH(s,k,h) model are:
#'
#' \itemize{
#'  \item{ar ~ normal(0,0.5)}
#'  \item{ma ~ normal(0,0.5)}
#'  \item{mu0 ~ t-student(0,2.5,6)}
#'  \item{sigma0 ~ t-student(0,1,7)}
#'  \item{arch ~ normal(0,0.5)}
#'  \item{garch ~ normal(0,0.5)}
#'  \item{mgarch ~ normal(0,0.5)}
#'  \item{dfv ~ gamma(2,0.1)}
#'  \item{breg ~ t-student(0,2.5,6)}
#' }
#'
#' For changing the default prior use the function \code{set_prior()}.
#'
#' @author Asael Alonzo Matamoros.
#'
#' @export
#' @importFrom stats as.ts time
#'
#' @references
#' Engle, R. (1982). Autoregressive Conditional Heteroscedasticity with Estimates of
#' the Variance of United Kingdom Inflation. \emph{Econometrica}, 50(4), 987-1007.
#' \code{url: http://www.jstor.org/stable/1912773}.
#'
#' Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity.
#' \emph{Journal of Econometrics}. 31(3), 307-327.
#' \code{doi: https://doi.org/10.1016/0304-4076(86)90063-1}.
#'
#' Fonseca, T. and Cequeira, V. and Migon, H. and Torres, C. (2019). The effects of
#' degrees of freedom estimation in the Asymmetric GARCH model with Student-t
#' Innovations. \emph{arXiv} \code{doi: arXiv: 1910.01398}.
#'
#' Ardia, D. and Hoogerheide, L. (2010). Bayesian Estimation of the GARCH(1,1) Model
#' with Student-t Innovations. \emph{The R Journal}. 2(7), 41-47.
#' \code{doi: 10.32614/RJ-2010-014}.
#'
#' @seealso \code{\link{Sarima}} \code{\link{auto.arima}} \code{\link{set_prior}}
#'
#' @examples
#' \dontrun{
#'  # Declaring a garch(1,1) model for the ipc data.
#'  sf1 = stan_garch(ipc,order = c(1,1,0))
#'
#'  # Declaring a t-student M-GARCH(2,3,1)-ARMA(1,1) process for the ipc data.
#'  sf2 = stan_garch(ipc,order = c(2,3,1),arma = c(1,1),genT = TRUE)
#'
#'  # Declaring a logistic Asymmetric GARCH(1,1) process.
#'  sf3 = stan_garch(ipc,order = c(1,1,0),asym = "logit")
#' }
#'
stan_garch = function(ts,order = c(1,1,0),arma = c(0,0),xreg = NULL,genT = FALSE,
                      asym = "none",chains = 4,iter = 4000,warmup = floor(iter/2),
                      adapt.delta = 0.9,tree.depth = 10,stepwise = TRUE,prior_mu0 = NULL,
                      prior_sigma0 = NULL,prior_ar = NULL, prior_ma = NULL,prior_mgarch = NULL,
                      prior_arch = NULL,prior_garch = NULL, prior_breg = NULL,prior_gamma = NULL,
                      prior_df = NULL,series.name = NULL,...){

  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  dat = garch(ts = ts,order = order,arma = arma,xreg = xreg,genT = genT,asym = asym,series.name = sn)

  # Priors selection
  if(!is.null(prior_mu0)) dat = set_prior(model = dat,par = "mu0",dist = prior_mu0)
  if(!is.null(prior_sigma0)) dat = set_prior(model = dat,par = "sigma0",dist = prior_sigma0)
  if(!is.null(prior_ar)) dat = set_prior(model = dat,par = "ar",dist = prior_ar)
  if(!is.null(prior_ma)) dat = set_prior(model = dat,par = "ma",dist = prior_ma)
  if(!is.null(prior_breg)) dat = set_prior(model = dat,par = "breg",dist = prior_breg)
  if(!is.null(prior_arch)) dat = set_prior(model = dat,par = "arch",dist = prior_arch)
  if(!is.null(prior_garch)) dat = set_prior(model = dat,par = "garch",dist = prior_garch)
  if(!is.null(prior_mgarch)) dat = set_prior(model = dat,par = "mgarch",dist = prior_mgarch)
  if(!is.null(prior_df)) dat = set_prior(model = dat,par = "df",dist = prior_df)
  if(!is.null(prior_gamma)) dat = set_prior(model = dat,par = "gamma",dist = prior_gamma)

  # Fitting the garch model.
  sf1 = varstan(model = dat,
                chains = chains,
                iter = iter,
                warmup = warmup,
                adapt.delta = adapt.delta,
                tree.depth = tree.depth,...)
  return(sf1)
}
#' Naive and Random Walk models.
#'
#' naive is the model constructor for a random walk  model applied to \code{y}.
#' This is equivalent to an ARIMA(0,1,0) model. \code{naive()} is simply a wrapper
#' to  maintain forecast package similitude. \code{seasonal} returns the model constructor
#' for a seasonal random walk equivalent to an ARIMA(0,0,0)(0,1,0)m model where m is the
#' seasonal period.
#'
#' The random walk with drift model is
#' \deqn{Y_t = mu_0 + Y_{t-1} + epsilon_t}{Y[t]= mu_0 +Y[t-1] + epsilon[t]}
#' where  \eqn{epsilon_t}{epsilon[t]} is a normal iid error.
#'
#' The seasonal naive model is
#' \deqn{Y_t = mu_0 + Y_{t-m} + epsilon_t}{Y[t]= mu_0 +Y[t-m] + epsilon[t]}
#' where  \eqn{epsilon_t}{epsilon[t]} is a normal iid error.
#'
#' @param ts  a numeric or ts object with the univariate time series.
#' @param seasonal a Boolean value for select a seasonal random walk instead.
#' @param m  an optional integer value for the seasonal period.
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
#' @param prior_mu0 The prior distribution for the location parameter in an ARIMA model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param prior_sigma0 The prior distribution for the scale parameter in an ARIMA model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param series.name an optional string vector with the series names.
#' @param ... Further arguments passed to  \code{varstan} function.
#'
#' @return The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @author Asael Alonzo Matamoros
#'
#' @seealso \code{\link{Sarima}}
#' @export
#'
#' @references
#'
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
#'  library(astsa)
#'  # A seasonal Random-walk model.
#'  sf1 = stan_naive(birth,seasonal = TRUE)
#' }
#'
stan_naive = function(ts,seasonal = FALSE,m = 0,chains = 4,iter = 4000,warmup = floor(iter/2),
                      adapt.delta = 0.9,tree.depth = 10,stepwise = TRUE,
                      prior_mu0 = NULL,prior_sigma0 = NULL,series.name = NULL,...){

  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  dat = naive(ts = ts,seasonal = seasonal,m = m)

  # Priors selection
  if(!is.null(prior_mu0)) dat = set_prior(model = dat,par = "mu0",dist = prior_mu0)
  if(!is.null(prior_sigma0)) dat = set_prior(model = dat,par = "sigma0",dist = prior_sigma0)

  # Fitting the naive model.
  sf1 = varstan(model = dat,
                chains = chains,
                iter = iter,
                warmup = warmup,
                adapt.delta = adapt.delta,
                tree.depth = tree.depth,...)
  return(sf1)
}
#' Fitting a Stochastic volatility model
#'
#' Fitting a Stochastic Volatility model (SVM) in \pkg{Stan}.
#'
#' The function returns a \code{varstan} object with the fitted model.
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param arma Optionally, a specification of the  ARMA model,same
#' as order parameter: the two components (p, q) are the AR order,and
#' the  MA order.
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
#' @param stepwise If TRUE, will do stepwise selection (faster). Otherwise, it searches
#' over all models. Non-stepwise selection can be very slow, especially for seasonal models.
#' @param prior_mu0 The prior distribution for the location parameter in an SVM model. By default
#' the value is set \code{NULL}, then the default normal(0,1) prior is used.
#' @param prior_sigma0 The prior distribution for the scale parameter in an SVM model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param prior_ar The prior distribution for the auto-regressive parameters in an ARMA model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_ma The prior distribution for the moving average parameters in an ARMA model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_alpha The prior distribution for the arch parameters in a GARCH model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_beta The prior distribution for the GARCH parameters in a GARCH model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_breg The prior distribution for the regression coefficient parameters in a
#' ARIMAX model. By default the value is set \code{NULL}, then the default student(7,0,1) priors
#' are used.
#' @param series.name an optional string vector with the series names.
#' @param ... Further arguments passed to  \code{varstan} function.
#'
#' @author Asael Alonzo Matamoros
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
#' @seealso \code{\link{garch}} \code{\link{set_prior}}
#'
#' @examples
#' \dontrun{
#'  # Declares a SVM model for the IPC data
#'  sf1 = stan_SVM(ipc,arma = c(1,1))
#' }
#'
stan_SVM = function(ts,arma = c(0,0),xreg = NULL,chains = 4,iter = 4000,
                    warmup = floor(iter/2),adapt.delta = 0.9,tree.depth = 10,
                    stepwise = TRUE,prior_mu0 = NULL,prior_sigma0 = NULL,prior_ar = NULL,
                       prior_ma = NULL, prior_alpha = NULL,prior_beta = NULL,
                       prior_breg = NULL,series.name = NULL,...){

  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  dat = SVM(ts = ts,arma = arma,xreg = xreg,series.name = sn)

  # Priors selection
  if(!is.null(prior_mu0)) dat = set_prior(model = dat,par = "mu0",dist = prior_mu0)
  if(!is.null(prior_sigma0)) dat = set_prior(model = dat,par = "sigma0",dist = prior_sigma0)
  if(!is.null(prior_ar)) dat = set_prior(model = dat,par = "ar",dist = prior_ar)
  if(!is.null(prior_ma)) dat = set_prior(model = dat,par = "ma",dist = prior_ma)
  if(!is.null(prior_alpha)) dat = set_prior(model = dat,par = "alpha",dist = prior_alpha)
  if(!is.null(prior_beta)) dat = set_prior(model = dat,par = "beta",dist = prior_beta)
  if(!is.null(prior_breg)) dat = set_prior(model = dat,par = "breg",dist = prior_breg)

  # Fitting the SVM model.
  sf1 = varstan(model = dat,
                chains = chains,
                iter = iter,
                warmup = warmup,
                adapt.delta = adapt.delta,
                tree.depth = tree.depth,...)
  return(sf1)
}
#' Fitting an Additive linear State space model.
#'
#' Fitting an Additive linear State space model in \pkg{Stan}.
#'
#' The function returns a \code{varstan} object with the fitted model.
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
#' @param stepwise If TRUE, will do stepwise selection (faster). Otherwise, it searches
#' over all models. Non-stepwise selection can be very slow, especially for seasonal models.
#' @param prior_sigma0 The prior distribution for the scale parameter in an SSM model. By default
#' the value is set \code{NULL}, then the default student(7,0,1) prior is used.
#' @param prior_level The prior distribution for the level parameter in a SSM model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_trend The prior distribution for the trend parameter in a SSM model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_damped The prior distribution for the damped trend parameter in a SSM model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_seasonal The prior distribution for the seasonal parameter in a SSM model.
#' By default the value is set \code{NULL}, then the default normal(0,0.5) priors are used.
#' @param prior_level1 The prior distribution for the initial level parameter in a SSM model.
#' By default the value is set \code{NULL}, then the default student(6,0,2.5) priors are used.
#' @param prior_trend1 The prior distribution for the initial trend parameter in a SSM model.
#' By default the value is set \code{NULL}, then the default student(6,0,2.5)  priors are used.
#' @param prior_seasonal1 The prior distribution for the initial seasonal parameters in a SSM model.
#' The prior is specified for the first m seasonal parameters, where m is the periodicity of the
#' defined time series. By default the value is set \code{NULL}, then the default normal(0,0.5) priors
#' are used.
#' @param prior_breg The prior distribution for the regression coefficient parameters in a
#' ARMAX model. By default the value is set \code{NULL}, then the default student(7,0,1) priors
#' are used.
#' @param prior_df The prior distribution for the degree freedom parameters in a t-student innovations
#' SSM model. By default the value is set \code{NULL}, then the default gamma(2,0.1) priors
#' are used.
#' @param series.name an optional string vector with the series names.
#' @param ... Further arguments passed to  \code{varstan} function.
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
#' \dontrun{
#'  # Declaring a local level model model for the ipc data.
#'  sf1 = stan_ssm(ipc)
#'
#'  # Declaring a Holt model for the ipc data.
#'  sf2 = stan_ssm(ipc,trend = TRUE,damped = TRUE)
#'
#'  # Declaring an additive Holt-Winters model for the birth data
#'  sf3 = stan_ssm(birth,trend = TRUE,damped = TRUE,seasonal = TRUE)
#' }
#'
stan_ssm = function(ts,trend = FALSE,damped = FALSE,seasonal = FALSE,xreg = NULL,
                    period = 0,genT = FALSE,chains = 4,iter = 4000,warmup = floor(iter/2),
                    adapt.delta = 0.9,tree.depth = 10,stepwise = TRUE,prior_sigma0 = NULL,
                    prior_level = NULL,prior_level1 = NULL, prior_trend = NULL,prior_trend1 = NULL,
                    prior_damped = NULL,prior_seasonal = NULL, prior_seasonal1 = NULL,prior_breg = NULL,
                    prior_df = NULL,series.name = NULL,...){

  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  dat = ssm(ts = ts,trend = trend,damped = damped,seasonal = seasonal,xreg = xreg,period = period,
            genT = genT,series.name = sn)

  # Priors selection
  if(!is.null(prior_sigma0)) dat = set_prior(model = dat,par = "sigma0",dist = prior_sigma0)
  if(!is.null(prior_level)) dat = set_prior(model = dat,par = "level",dist = prior_level)
  if(!is.null(prior_level1)) dat = set_prior(model = dat,par = "level1",dist = prior_level1)
  if(!is.null(prior_trend)) dat = set_prior(model = dat,par = "trend",dist = prior_trend)
  if(!is.null(prior_trend1)) dat = set_prior(model = dat,par = "trend1",dist = prior_trend1)
  if(!is.null(prior_damped)) dat = set_prior(model = dat,par = "damped",dist = prior_damped)
  if(!is.null(prior_seasonal)) dat = set_prior(model = dat,par = "seasonal",dist = prior_seasonal)
  if(!is.null(prior_seasonal1)) dat = set_prior(model = dat,par = "seasonal1",dist = prior_seasonal1)
  if(!is.null(prior_breg)) dat = set_prior(model = dat,par = "breg",dist = prior_breg)
  if(!is.null(prior_df)) dat = set_prior(model = dat,par = "df",dist = prior_df)

  # Fitting the SSM model.
  sf1 = varstan(model = dat,
                chains = chains,
                iter = iter,
                warmup = warmup,
                adapt.delta = adapt.delta,
                tree.depth = tree.depth,...)
  return(sf1)
}
