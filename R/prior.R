#' Set a prior distribution to a model parameter.
#'
#' setting a prior distribution to an specify model parameter.
#'
#' @param model a time series model class specified in \pkg{varstan}.
#' @param par a string value with the  desired parameter which a prior is defined could be:
#' "mu0", "sigma0", "ar", "ma", "arch", "garch", "mgarch", "dfv", "df", "LKJ" or "breg".
#' @param dist the distribution of the prior parameter. The only accepted is a prior_dist object.
#' @param lag an optional integer value, indicates the desired lag of the parameter which the prior
#' is defined if lag = 0, then the prior distribution will be applied for all lags.
#'
#' @return a time series model class specified in \pkg{varstan} with the changed prior.
#'
#' @details
#' varstan provides its own functions to manipulate the parameter prior, this functions return
#' a \code{prior_dist} class, the \code{dist} argument only accepts this objects.
#'
#' \code{lag} parameter is an optional value to change te prior distribution of one parameter in particular,
#' this argument is only valid for: "ar","ma", "arch", "garch", "mgarch", or "breg" par arguments. lag has to
#' be a integer lower than the total amount of lagged parameters of the model. For example, to  ONLY
#' change the prior of the second "arch" paramater in a garch(3,1) model, a lag = 2 values must be specified.
#'
#' For varma and Bekk models the covariance matrix Sigma is factorized as follows:
#'
#'                   Sigma = D' Omega D
#'
#' Where Omega is the correlation matrix that accepts an LKJ prior distribution D is a diagonal matrix with
#' the inverse std desviations
#'
#' For changing the degree freedom in a LKJ distribution for omega use par = "LKJ" and dist = LKJ(df),
#' where df are the desired degree freedom.
#'
#' For changing the the priors in the diagonal D use par = "sigma0" and select one of the available prior
#' distributions.
#'
#' For ar, ma garch, arch parameters in varma and Bekk models only normal distributions priors with different
#' mu and sd are accepted. Even if \code{get_prior} accepts its change, Stan will change it to a normal(0,1) prior.
#'
#' @author Asael Alonzo Matamoros
#' @export
#' @examples
#' library(astsa)
#' dat = Sarima(birth,order = c(1,1,2))
#' dat = set_prior(model = dat,par = "ar",dist = normal(0,2))
#' dat
#'
#' dat = set_prior(model = dat,par = "mu0",dist = student(mu=0,sd = 2.5,df = 7))
#' dat
#'
#' dat = set_prior(model = dat,par = "ma",dist= beta(shape1 = 2,shape2 = 2),lag = 2)
#' dat
#'
set_prior = function(model,par = "ar",dist = normal(),lag = 0){
  if(!is.model(model))
    stop("The defined model does not belong to varstan")

  if(!check_par(par))
    stop("par is not a defined parameter")

  if(!is(object = dist,class2 = "prior_dist"))
    stop("dist is not a prior distribution object")

  if(!check_dist1(par = par,dist = dist$x))
    stop("The parameter does not accepts the defined distribution")

  if(is.garch(model)){
    if( identical(model$asym1,0) && identical(par,"gamma"))
      stop("The model is not an asymmetric GARCH")
  }

  if(lag < 0 )
    stop("lag values lower than zero are not accepted")

  if(identical(lag,0)){
    if(identical(par,"ar"))    model$prior_ar         = matrix(rep(dist$x,model$p),ncol = 4,byrow = TRUE)
    if(identical(par,"ma"))    model$prior_ma         = matrix(rep(dist$x,model$q),ncol = 4,byrow = TRUE)
    if(identical(par,"sar"))   model$prior_sar        = matrix(rep(dist$x,model$P),ncol = 4,byrow = TRUE)
    if(identical(par,"sma"))   model$prior_sma        = matrix(rep(dist$x,model$Q),ncol = 4,byrow = TRUE)
    if(identical(par,"arch"))  model$prior_arch       = matrix(rep(dist$x,model$s),ncol = 4,byrow = TRUE)
    if(identical(par,"garch")) model$prior_garch      = matrix(rep(dist$x,model$k),ncol = 4,byrow = TRUE)
    if(identical(par,"mgarch"))model$prior_mgarch     = matrix(rep(dist$x,model$q),ncol = 4,byrow = TRUE)
    if(identical(par,"breg"))  model$prior_breg       = matrix(rep(dist$x,model$d1),ncol = 4,byrow = TRUE)
    if(identical(par,"gamma")) model$prior_gamma      = matrix(rep(dist$x,2),ncol = 4,byrow = TRUE)
    if(identical(par,"mu0"))   model$prior_mu0 = dist$x
    if(identical(par,"sigma0"))model$prior_sigma0 = dist$x
    if(identical(par,"dfv"))   model$prior_dfv = dist$x
    if(identical(par,"df"))    model$prior_dfv = dist$x
    if(identical(par,"LKJ"))   model$prior_lkj = dist$x
    if(identical(par,"alpha")) model$prior_alpha = dist$x
    if(identical(par,"beta"))  model$prior_beta = dist$x
    if(identical(par,"level")) model$prior_level = dist$x
    if(identical(par,"trend")) model$prior_trend = dist$x
    if(identical(par,"damped"))   model$prior_damped = dist$x
    if(identical(par,"seasonal")) model$prior_seasonal = dist$x
    if(identical(par,"level1"))   model$prior_level1 = dist$x
    if(identical(par,"trend1"))   model$prior_trend1 = dist$x
    if(identical(par,"seasonal1"))model$prior_seasonal1 = dist$x
  }
  else{
    if(identical(par,"ar"))    if(lag <= model$p)  model$prior_ar[lag,]= dist$x
    if(identical(par,"ma"))    if(lag <= model$q)  model$prior_ma[lag,] = dist$x
    if(identical(par,"sar"))   if(lag <= model$P)  model$prior_sar[lag,]= dist$x
    if(identical(par,"sma"))   if(lag <= model$Q)  model$prior_sma[lag,]= dist$x
    if(identical(par,"arch"))  if(lag <= model$s)  model$prior_arch[lag,]= dist$x
    if(identical(par,"garch")) if(lag <= model$k)  model$prior_garch[lag,]= dist$x
    if(identical(par,"mgarch"))if(lag <= model$h)  model$prior_mgarch[lag,]= dist$x
    if(identical(par,"breg"))  if(lag <= model$d1) model$prior_breg[lag,]= dist$x
    if(identical(par,"gamma")) if(lag <= 2)        model$prior_gamma[lag,]= dist$x
    if(identical(par,"mu0"))   model$prior_mu0 = dist$x
    if(identical(par,"sigma0"))model$prior_sigma0 = dist$x
    if(identical(par,"dfv"))   model$prior_dfv = dist$x
    if(identical(par,"df"))    model$prior_dfv = dist$x
    if(identical(par,"LKJ"))   model$prior_lkj = dist$x
    if(identical(par,"alpha")) model$prior_alpha = dist$x
    if(identical(par,"beta"))  model$prior_beta = dist$x
    if(identical(par,"level")) model$prior_level = dist$x
    if(identical(par,"trend")) model$prior_trend = dist$x
    if(identical(par,"damped"))   model$prior_damped = dist$x
    if(identical(par,"seasonal")) model$prior_seasonal = dist$x
    if(identical(par,"level1"))   model$prior_level1 = dist$x
    if(identical(par,"trend1"))   model$prior_trend1 = dist$x
    if(identical(par,"seasonal1"))model$prior_seasonal1 = dist$x
  }
  return(model)
}
#' Get the prior distribution of a model parameter
#'
#' The functions gets the defined distribution of a defined model parameter
#'
#' @usage get_prior(model,par,lag = 0)
#'
#' @param model a time series model class specified in varstan.
#' @param par a string value with the  desired parameter which a prior is defined could be:
#' "mu0", "sigma0", "ar", "ma", "arch", "garch", "mgarch", "dfv", "df", "LKJ" or "breg".
#' @param lag an optional integer value, indicates the desired lag of the parameter which the prior
#' is defined if lag = 0, then the prior distribution will be applied for all lags
#'
#' @author  Asael Alonzo Matamoros
#'
#' @export
#'
#' @examples
#' library(astsa)
#' # get all the ar parameters
#' dat = Sarima(birth,order = c(2,1,2))
#' get_prior(model = dat,par = "ar")
#'
#' # change the mean constant parameter
#' dat = set_prior(model = dat,par = "mu0",dist = student(0,2.5,7))
#' get_prior(dat,par = "mu0")
#'
#' # change and print only the second ma parameter
#' dat = set_prior(model = dat,par = "ma",dist = beta(2,2),lag = 2)
#' get_prior(dat,par = "ma")
#'
get_prior = function(model,par,lag = 0){
  if(!is.model(model))
    stop("The defined model does not belong to varstan")

  if(!check_par(par))
    stop("par is not a defined parameter")

  if(is.garch(model)){
    if( identical(model$asym1,0) && identical(par,"gamma"))
      stop("The model is not an asymmetric GARCH")
  }

  if(lag < 0 )
    stop("lag values lower than zero are not accepted")

  y = NULL

  if( identical(lag,0) ){
    if( identical(par,"ar"))     if(model$p > 0) y = print_dist_matrix(model$prior_ar,par = "ar")
    if( identical(par,"ma"))     if(model$q > 0) y = print_dist_matrix(model$prior_ma,par = "ma")
    if( identical(par,"sar"))    if(model$P > 0) y = print_dist_matrix(model$prior_sar,par = "sar")
    if( identical(par,"sma"))    if(model$Q > 0) y = print_dist_matrix(model$prior_sma,par = "sma")
    if( identical(par,"arch"))   if(model$s > 0) y = print_dist_matrix(model$prior_arch,par = "arch")
    if( identical(par,"garch"))  if(model$k > 0) y = print_dist_matrix(model$prior_garch,par = "garch")
    if( identical(par,"mgarch")) if(model$h > 0) y = print_dist_matrix(model$prior_mgarch,par = "mgarch")
    if( identical(par,"breg"))   if(model$d1 > 0) y = print_dist_matrix(model$prior_breg,par = "breg")
    if( identical(par,"gamma"))  if(model$asym1) y = print_dist_matrix(model$prior_gamma,par = "gamma")
    if( identical(par,"mu0"))    y = print_dist_numeric(model$prior_mu0,"mu0")
    if( identical(par,"sigma0")) y = print_dist_numeric(model$prior_sigma0,"sigma0")
    if( identical(par,"dfv"))    y = print_dist_numeric(model$prior_dfv,"df")
    if( identical(par,"df"))     y = print_dist_numeric(model$prior_dfv,"df")
    if( identical(par,"LKJ"))    y = print_dist_numeric(model$prior_lkj,"LKJ")
    if( identical(par,"alpha"))  y = print_dist_numeric(model$prior_alpha,"alpha")
    if( identical(par,"beta"))   y = print_dist_numeric(model$prior_beta,"beta")
    if( identical(par,"level"))  y = print_dist_numeric(model$prior_level,"level")
    if( identical(par,"level1")) y = print_dist_numeric(model$prior_level1,"level1")
    if( identical(par,"trend"))  if(model$is_td) y = print_dist_numeric(model$prior_trend,"trend")
    if( identical(par,"trend1")) if(model$is_td) y = print_dist_numeric(model$prior_trend1,"trend1")
    if( identical(par,"damped")) if(model$is_dp) y = print_dist_numeric(model$prior_damped,"damped")
    if( identical(par,"seasonal")) if(model$is_ss) y = print_dist_numeric(model$prior_seasonal,"seasonal")
    if( identical(par,"seasonal1"))if(model$is_ss) y = print_dist_numeric(model$prior_seasonal1,"seasonal1")
  }
  else{
    if( identical(par,"ar") )    if(lag <= model$p & model$p > 0) y = print_dist_numeric(model$prior_ar[lag,],par = "ar",lag = lag)
    if(identical(par,"ma") )     if(lag <= model$q & model$q > 0) y = print_dist_numeric(model$prior_ma[lag,],par = "ma",lag = lag)
    if(identical(par,"sar") )    if(lag <= model$P & model$P > 0) y = print_dist_numeric(model$prior_sar[lag,],par = "sar",lag = lag)
    if(identical(par,"sma") )    if(lag <= model$Q & model$Q > 0) y = print_dist_numeric(model$prior_sma[lag,],par = "sma",lag = lag)
    if(identical(par,"arch") )   if(lag <= model$s & model$s > 0) y = print_dist_numeric(model$prior_arch[lag,],par = "arch",lag = lag)
    if(identical(par,"garch") )  if(lag <= model$k & model$k > 0) y = print_dist_numeric(model$prior_garch[lag,],par = "garch",lag = lag)
    if(identical(par,"mgarch") ) if(lag <= model$h & model$h > 0) y = print_dist_numeric(model$prior_mgarch[lag,],par = "mgarch",lag = lag)
    if(identical(par,"breg") )   if(lag <= model$d1 & model$d1 > 0) y = print_dist_numeric(model$prior_breg[lag,],par = "breg",lag = lag)
    if(identical(par,"gamma") )  if(lag <= 2 & model$asym1 > 0)     y = print_dist_numeric(model$prior_gamma[lag,],par = "gamma",lag = lag)
    if(identical(par,"mu0") )    y = print_dist_numeric(model$prior_mu0,"mu0")
    if(identical(par,"sigma0") ) y = print_dist_numeric(model$prior_sigma0,"sigma0")
    if(identical(par,"dfv") )    y = print_dist_numeric(model$prior_dfv,"df")
    if(identical(par,"df") )     y = print_dist_numeric(model$prior_dfv,"df")
    if(identical(par,"LKJ") )    y = print_dist_numeric(model$prior_lkj,"LKJ")
    if( identical(par,"alpha"))  y = print_dist_numeric(model$prior_alpha,"alpha")
    if( identical(par,"beta"))   y = print_dist_numeric(model$prior_beta,"beta")
    if( identical(par,"level"))  y = print_dist_numeric(model$prior_level,"level")
    if( identical(par,"level1")) y = print_dist_numeric(model$prior_level1,"level1")
    if( identical(par,"trend"))  if(model$is_td) y = print_dist_numeric(model$prior_trend,"trend")
    if( identical(par,"trend1")) if(model$is_td) y = print_dist_numeric(model$prior_trend1,"trend1")
    if( identical(par,"damped")) if(model$is_dp) y = print_dist_numeric(model$prior_damped,"damped")
    if( identical(par,"seasonal")) if(model$is_ss) y = print_dist_numeric(model$prior_seasonal,"seasonal")
    if( identical(par,"seasonal1"))if(model$is_ss) y = print_dist_numeric(model$prior_seasonal1,"seasonal1")
  }
  print(y)
}

######################### Distributions #######################################

#' Define a normal prior distribution
#'
#' normal(mu,sd)
#'
#' @param mu the location parameter mu
#' @param sd the standard desviation parameter sigma
#'
#' @details
#' Define a normal prior distribution using the hyper parameters
#' mu and sigma, by default a standard normal distribution is
#' return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
normal = function(mu = 0,sd = 1){
  #Normal == 1
  m = list(x = c(check_loc(mu),check_scl(sd),1,1))
  attr(m,"class") = c("prior_dist","normal")
  return(m)
}
#' @method print normal
#' @export
#' @noRd
#'
print.normal = function(x,...){
  if(!inherits(x,what = "normal"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( mu = ",x$x[1],",sd = ",x$x[2],")" )
}
#' Define a beta prior distribution
#'
#' beta(shape1,shape2)
#'
#' @param shape1 the first form parameter
#' @param shape2 the second form parameter
#'
#' @details
#' Define a beta prior distribution using the hyper parameters
#' shape1 and shape2, by default a beta(2,2) distribution is
#' return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
beta= function(shape1 = 2,shape2 = 2){
  #Beta == 1
  m = list(x = c(check_form(shape1),check_form(shape2),1,2))
  attr(m,"class") = c("prior_dist","beta")
  return(m)
}
#' @method print beta
#' @export
#' @noRd
#'
print.beta = function(x,...){
  if(!inherits(x,what = "beta"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( shape1 = ",x$x[1],",shape2 = ",x$x[2],")" )
}
#' Define a uniform prior distribution
#'
#' uniform(shape1,shape2)
#'
#' @param min the first form parameter
#' @param max the second form parameter
#'
#' @details
#' Define a beta prior distribution using the hyper parameters
#' min and max, by default a uniform(0,1) distribution is
#' return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
uniform= function(min = 0,max = 1){
  #unif == 3
  if(min > max) x = c(check_form(max),check_form(min),1,3)
  else if(min==max) x = c(0,1,1,3)
  else x = c(check_form(min),check_form(max),1,3)

  m = list(x = x)
  attr(m,"class") = c("prior_dist","uniform")
  return(m)
}
#' @method print uniform
#' @export
#' @noRd
#'
print.uniform = function(x,...){
  if(!inherits(x,what = "uniform"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( min = ",x$x[1],",max = ",x$x[2],")" )
}
#' Define a t student prior distribution
#'
#' student(mu,sd)
#'
#' @param mu the location parameter mu
#' @param sd the standard desviation parameter sigma
#' @param df the degree freedom parameter df
#'
#' @details
#' Define a t student prior distribution using the hyper parameters
#' mu, sigma and df as degree freedom, by default a standard t-sutdent(0,1,5)
#' distribution with 5 degree freedom is return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
student = function(mu= 0,sd = 1,df = 5){
  #student == 4
  m = list(x = c(check_loc(mu),check_scl(sd),check_df(df),4))
  attr(m,"class") = c("prior_dist","student")
  return(m)
}
#' @method print student
#' @export
#' @noRd
#'
print.student = function(x,...){
  if(!inherits(x,what = "student"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( mu = ",x$x[1],",sd = ",x$x[2],",df = ",x$x[3],")" )
}
#' Define a Cauchy prior distribution
#'
#' cauchy(mu,sd)
#'
#' @param mu the location parameter mu
#' @param sd the standard desviation parameter sigma
#'
#' @details
#' Define a cauchy prior distribution using the hyper parameters
#' mu and sigma, by default a standard cauchy(0,1) distribution is
#' return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
cauchy = function(mu= 0,sd = 1){
  #cauchy == 5
  m = list(x = c(check_loc(mu),check_scl(sd),1,5))
  attr(m,"class") = c("prior_dist","cauchy")
  return(m)
}
#' @method print cauchy
#' @export
#' @noRd
#'
print.cauchy = function(x,...){
  if(!inherits(x,what = "cauchy"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( mu = ",x$x[1],",sd = ",x$x[2],")" )
}
#' Define an inverse gamma prior distribution
#'
#' inverse.gamma(shape,rate)
#'
#' @param shape the form parameter alpha in gamma distribution
#' @param rate  the rate parameter beta in gamma distribution
#'
#' @details
#' Define a inverse.gamma prior distribution using the hyper parameters
#' shape and rate, by default an inverse.gamma(2,1) distribution is
#' return.
#'
#' If sigma has a gamma distribution then 1/sigma has n  inverse gamma
#' distribution. The rate parameter is the inverse of an scale parameter.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
inverse.gamma = function(shape= 2,rate = 1){
  #inverse_gamma == 6
  m = list(x = c(check_form(shape),check_form(rate),1,6))
  attr(m,"class") = c("prior_dist","inverse.gamma")
  return(m)
}
#' @method print inverse.gamma
#' @export
#' @noRd
#'
print.inverse.gamma = function(x,...){
  if(!inherits(x,what = "inverse.gamma"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( shape = ",x$x[1],",rate = ",x$x[2],")" )
}
#' Define an inverse gamma prior distribution
#'
#' inverse.chisq(df)
#'
#' @param df the degree freedom parameter df
#'
#' @details
#' Define a inverse chi square prior distribution using the hyper parameter
#' df, by default an inverse.chisq(df = 2) distribution is return.
#'
#' If sigma has a chi square distribution then 1/sigma has n  inverse chi square
#' distribution.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
inverse.chisq = function(df = 7){
  #inverse_chi_square == 7
  m = list(x = c(0,0,check_df(df),7))
  attr(m,"class") = c("prior_dist","inverse.chisq")
  return(m)
}
#' @method print inverse.chisq
#' @export
#' @noRd
#'
print.inverse.chisq = function(x,...){
  if(!inherits(x,what = "inverse.chisq"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( df = ",x$x[3],")" )
}
#' Define a non informative Jeffrey's prior for the degree freedom hyper parameter
#'
#' jeffey.df()
#'
#' @details
#' Define a non informative Jeffrey's prior distribution, by default an jeffrey.df( )
#' distribution is return.
#'
#' This prior can only be used in garch models with t-student innovations, or Bekk
#' models with generalized t-student distribution.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
jeffrey = function(){
  #Jeffrey == 8
  warning("This prior could cause ill estimations")
  m = list(x = c(0,1,1,8))
  attr(m,"class") = c("prior_dist","jeffrey")
  return(m)
}
#' @method print jeffrey
#' @export
#' @noRd
#'
print.jeffrey = function(x,...){
  if(!inherits(x,what = "jeffrey"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( )" )
}
#' Define a gamma prior distribution
#'
#' gamma(shape,rate)
#'
#' @param shape the form parameter alpha in gamma distribution
#' @param rate  the rate parameter beta in gamma distribution
#'
#' @details
#' Define a gamma prior distribution using the hyper parameters
#' shape and rate, by default an gamma(2,1) distribution is
#' return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
gamma = function(shape= 2,rate = 1){
  #gamma == 9
  m = list(x = c(check_form(shape),check_form(rate),1,9))
  attr(m,"class") = c("prior_dist","gamma")
  return(m)
}
#' @method print gamma
#' @export
#' @noRd
#'
print.gamma = function(x,...){
  if(!inherits(x,what = "gamma"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( shape = ",x$x[1],",rate = ",x$x[2],")" )
}
#' Define an exponetial prior distribution
#'
#' exponetial(rate)
#'
#' @param rate  the rate parameter lambda in exponential distribution
#'
#' @details
#' Define a gamma prior distribution using the rate hyper parameter,
#' by default an exponential(1) distribution is return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
exponential= function(rate = 1){
  #exponential == 10
  m = list(x = c(1,check_form(rate),1,10))
  attr(m,"class") = c("prior_dist","exponential")
  return(m)
}
#' @method print exponential
#' @export
#' @noRd
#'
print.exponential = function(x,...){
  if(!inherits(x,what = "exponential"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( rate = ",x$x[2],")" )
}
#' Define a chi square prior distribution
#'
#' chisq(df)
#'
#' @param df the degree freedom parameter df
#'
#' @details
#' Define a gamma prior distribution using the degree freedom df hyper parameter,
#' by default an chisq(7) distribution is return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
chisq = function(df = 7){
  #chisq == 11
  m = list(x = c(0,1,check_df(df),11))
  attr(m,"class") = c("prior_dist","chisq")
  return(m)
}
#' @method print chisq
#' @export
#' @noRd
#'
print.chisq = function(x,...){
  if(!inherits(x,what = "chisq"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( df = ",x$x[3],")" )
}
#' Define a Laplace prior distribution
#'
#' laplace(mu,sd)
#'
#' @param mu the location parameter mu
#' @param sd the standard desviation parameter sigma
#'
#' @details
#' Define a Laplace prior distribution using the hyper parameters
#' mu and sigma, by default a standard Laplace distribution is
#' return.
#'
#' The laplace distribution is exactly the same as the double exponential distribution
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
laplace = function(mu = 0,sd = 1){
  # laplace == 12
  m = list(x = c(check_loc(mu),check_scl(sd),1,12))
  attr(m,"class") = c("prior_dist","laplace")
  return(m)
}
#' @method print laplace
#' @export
#' @noRd
#'
print.laplace = function(x,...){
  if(!inherits(x,what = "laplace"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( mu = ",x$x[1],",sd = ",x$x[2],")" )
}
#' Define a LKJ matrix prior distribution
#'
#' LKJ(df)
#'
#' @param df the degree freedom parameter df
#'
#' @details
#' Define a Lewandowski Kurowicka and Joe (LKJ) matrix correlation prior
#' distribution using the degree freedom df hyper parameter,by default
#' a LKJ(2) distribution is return.
#'
#' @return a numerical vector interpreted as a prior in Stan
#'
#' @export
#'
#' @references
#' Lewandowski D, Kurowicka D, Joe H (2009). "Generating random correlation matrices based
#' on vines and extended onion method." Journal of Multivariate Analysis, 100(9), 1989 2001.
#' ISSN 0047-259X. doi:https://doi.org/10.1016/j.jmva.2009.04.008.
#' URL: http://www.sciencedirect.com/science/article/pii/S0047259X09000876.
#'
LKJ = function(df = 2){
  # LKJ == 13
  m = list(x = c(check_df(df),check_df(df),check_df(df),13))
  attr(m,"class") = c("prior_dist","LKJ")
  return(m)
}
#' @method print LKJ
#' @export
#' @noRd
#'
print.LKJ = function(x,...){
  if(!inherits(x,what = "LKJ"))
    stop("The current class is not a prior_dist object")
  cat(get.dist(x$x[4]),"( df = ",x$x[1],")" )
}

###############################################################################################
#   Prior internals
###############################################################################################


#'
#' get the distribution character  string of vector
#' @noRd
#'
print_dist_numeric = function(x,par,lag = NULL){
  if(!is.numeric(x))
    stop("the value must be a numeric class")

  y = NULL
  if(identical(x[4],1)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( mu = ",x[1],",sd = ",x[2],")")
  if(identical(x[4],2)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( shape1 = ",x[1],",shape2 = ",x[2],")")
  if(identical(x[4],3)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( min = ",x[1],",max = ",x[2],")")
  if(identical(x[4],4)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( mu = ",x[1],",sd = ",x[2],",df = ",x[3],")")
  if(identical(x[4],5)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( mu = ",x[1],",sd = ",x[2],")")
  if(identical(x[4],6)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( shape = ",x[1],",rate = ",x[2],")")
  if(identical(x[4],7)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( df = ",x[3],")")
  if(identical(x[4],8)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( )")
  if(identical(x[4],9)) y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( shape = ",x[1],",rate = ",x[2],")")
  if(identical(x[4],10))y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( rate = ",x[2],")")
  if(identical(x[4],11))y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( df = ",x[3],")")
  if(identical(x[4],12))y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( mu = ",x[1],",sd = ",x[2],")")
  if(identical(x[4],13))y = paste(par,"[",lag,"] ~",get.dist(x[4]),"( df = ",x[1],")")

  return(y)
}
#'
#' get the distribution character  string of matrix
#' @noRd
#'
print_dist_matrix = function(x,par){
  if(!is.matrix(x))
    stop("the value must be a matrix class")
  n = nrow(x);y = NULL
  if(n <= 0)
    return(y)

  for (i in 1:n){
    y = c(y,print_dist_numeric(x[i,],par,i))
  }
  return(matrix(y,ncol = 1))
}

#'
#' get the distribution name
#' @noRd
#'
get.dist = function(code = 1){
  if( !(code %in% 1:13) )
    stop("code value must be a value between 1 and 13")

  dist= c("normal","beta","uniform","student","cauchy","inverse.gamma","inverse.chisq",
          "jeffrey.df","gamma","exponential","chisq","laplace","LKJ")

  return(dist[code])
}

#'
#' get the distribution code
#' @noRd
#'
get.code = function(dist = "normal"){
  dist1= c("normal","beta","uniform","student","cauchy","inverse.gamma","inverse.chisq",
           "jeffrey.df","gamma","exponential","chisq","laplace","LKJ")

  code = match(x = dist,table = dist1)

  if(is.na(code))
    stop("The current distribution is not available")

  return(code)
}

#' Check if the value is in the domain of a scale parameter
#'
#' @param  dist the distribution
#' @param  par the parameter
#'
#' @noRd
#'
check_dist1 = function(par,dist){
  y = FALSE
  if(!is.na(match(table = c("ma","ar","sma","sar","arch","garch","beta","level","trend","damped","seasonal"),x=par))){
    if(dist[4] %in% c(1,2,3)) y = TRUE
  }
  else if(!is.na(match(table = c("mu0","sigma0","mgarch","breg","dfv","df","alpha","gamma","level1","trend1","seasonal1"),x = par))){
    if(dist[4] %in% 1:12) y = TRUE
  }
  else if(identical(par,"LKJ")) if(identical(dist[4],13)) y = 13

  return(y)
}

#' Check if the value is a type of parameter
#'
#' @param x a parameter
#' @noRd
#'
check_par <- function(x) {
  y = FALSE
  if(identical(x,"sma"))    y = TRUE
  if(identical(x,"ma"))     y = TRUE
  if(identical(x,"ar"))     y = TRUE
  if(identical(x,"sar"))    y = TRUE
  if(identical(x,"arch"))   y = TRUE
  if(identical(x,"garch"))  y = TRUE
  if(identical(x,"mgarch")) y = TRUE
  if(identical(x,"mu0"))    y = TRUE
  if(identical(x,"breg"))   y = TRUE
  if(identical(x,"sigma0")) y = TRUE
  if(identical(x,"dfv"))    y = TRUE
  if(identical(x,"df"))     y = TRUE
  if(identical(x,"alpha"))  y = TRUE
  if(identical(x,"beta"))   y = TRUE
  if(identical(x,"LKJ"))    y = TRUE
  if(identical(x,"gamma"))  y = TRUE
  if(identical(x,"level"))  y = TRUE
  if(identical(x,"level1")) y = TRUE
  if(identical(x,"trend"))  y = TRUE
  if(identical(x,"trend1")) y = TRUE
  if(identical(x,"damped"))   y = TRUE
  if(identical(x,"seasonal")) y = TRUE
  if(identical(x,"seasonal1"))y = TRUE
  return(y)
}

