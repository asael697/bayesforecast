#' A  constructor for a GARCH(s,k,h) model.
#'
#' Constructor of the \code{GARCH(s,k,h)} object for Bayesian estimation in \pkg{Stan}.
#'
#' The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
#'
#' @usage garch(ts,order = c(1,1,0),arma = c(0,0),xreg = NULL,
#'              genT = FALSE,asym = "none",series.name = NULL)
#'
#' @param ts a numeric or ts object with the univariate time series.
#' @param order A specification of the garch  model: the three components (s, k, h)
#' are the arch order, the garch order, and the mgarch order.
#' @param arma A specification of the  ARMA model,same as order parameter:  the two
#' components (p, q) are the AR order,and the  MA order.
#' @param xreg	Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param genT a boolean value to specify for a generalized t-student garch model.
#' @param asym a string value for the asymmetric function for an asymmetric GARCH process. By default
#' the value \code{"none"} for standard GARCH process. If \code{"logit"} a logistic function
#' is used for asymmetry, and if \code{"exp"} an exponential function is used.
#' @param series.name an optional string vector with the time series names.
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
#' @return The function returns a list with the data for running \code{stan()} function of
#'  \pkg{rstan} package.
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
#' # Declaring a garch(1,1) model for the ipc data.
#' dat = garch(ipc,order = c(1,1,0))
#' dat
#'
#' # Declaring a t-student M-GARCH(2,3,1)-ARMA(1,1) process for the ipc data.
#' dat = garch(ipc,order = c(2,3,1),arma = c(1,1),genT = TRUE)
#' dat
#'
#' # Declaring a logistic Asymmetric GARCH(1,1) process.
#' dat = garch(ipc,order = c(1,1,0),asym = "logit")
#' dat
#'
garch = function(ts,order = c(1,1,0),arma = c(0,0),xreg = NULL,
                 genT = FALSE,asym = "none",series.name = NULL){

  n = length(as.numeric(ts))
  y = as.numeric(ts)

  # series name
  if(is.null(series.name))
    sn = deparse(substitute(ts))
  else
    sn = as.character(series.name)

  m1 = list(n = n,dimension = 1,time = as.numeric(stats::time(ts)),
            s = no_negative_check(order[1] ),
            k = no_negative_check(order[2]),
            h = no_negative_check(order[3]),
            y = y,yreal = stats::as.ts(ts),series.name = sn)

  m1$prior_mu0 = c(0,1,0,1)
  m1$prior_sigma0 = c(0,1,7,4)
  m1$prior_arch    = matrix(rep(c(0,0.5,1,1),order[1]),ncol = 4,byrow = TRUE)
  m1$prior_garch   = matrix(rep(c(0,0.5,1,1),order[2]),ncol = 4,byrow = TRUE)
  m1$prior_mgarch  = matrix(rep(c(0,0.5,1,1),order[3]),ncol = 4,byrow = TRUE)

  # arma representation
  m1$p = arma[1]
  m1$q = arma[2]
  m1$prior_ar =  matrix(rep(c(0,0.5,1,1),arma[1]),ncol = 4,byrow = TRUE)
  m1$prior_ma =  matrix(rep(c(0,0.5,1,1),arma[2]),ncol = 4,byrow = TRUE)

  # Generalized t distribution
  m1$genT = genT
  m1$prior_dfv = c(2,0.1,1,9)

  # arima regression model
  if( !is.null(xreg) ){

    if(!is.matrix(xreg))
      stop("xreg has to be a matrix with row dimension as same as the length of the time serie")

    if(nrow(xreg) != n)
      stop("The length of xreg don't match with the length of the time serie")

    m1$d1 = ncol(xreg)
    m1$xreg = xreg
  }
  else{
    m1$d1 = 0
    m1$xreg = matrix(rep(0,m1$d1*n),ncol = m1$d1,nrow = n)
  }
  m1$prior_breg  = matrix(rep(c(0,2.5,6,4),m1$d1),ncol = 4,byrow = TRUE)

  # asymmetric GARCH
  m1$asym = check_asym(asym)
  m1$asym1 = ifelse(m1$asym > 0,1,0)
  m1$prior_gamma = matrix(rep(c(0,0.5,1,1),2*m1$asym1),ncol = 4,byrow = TRUE)

  attr(m1,"class") = "garch"
  return(m1)
}
#' Checks if is a garch object
#'
#' @param object a  garch object
#' @noRd
#'
is.garch = function(object){
  y = FALSE
  if(is(object,"garch")) y = TRUE
  return (y)
}
#' Extracts all the order coefficients in a list
#'
#' @param dat A garch model
#' @noRd
#'
get_order_garch= function(dat){
    return(list(p = dat$p,q=dat$q,d1 = dat$d1,s=dat$s,k=dat$k,h=dat$h))
}
#' Max order  coefficients in a garch model
#'
#' @param dat A garch model
#' @noRd
#'
max_order_garch= function(dat){

  return(max(c(dat$p,dat$q,dat$s,dat$k,dat$h)))
}
