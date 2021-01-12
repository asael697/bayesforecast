#' Constructor of a varstan object.
#'
#' Constructor of the varstan object for Bayesian estimation in \pkg{Stan}.
#'
#' The function estimates one of  the defined models in \pkg{Stan} using
#' the \code{stan()} function for sampling.
#'
#' @usage varstan(model,chains=4,iter=2000,warmup=floor(iter/2),
#'                adapt.delta = 0.90,tree.depth =10,...)
#'
#' @param model One of the \code{varstan} model classes defined in the package.
#' @param chains An integer of the number of Markov Chains chains to be run,
#' by default 4 chains are run.
#' @param iter An integer of total iterations per chain including the warm-up,
#' by default  the number of iterations are 2000.
#' @param warmup  A positive integer specifying number of warm-up (aka burn-in)
#'   iterations. This also specifies the number of iterations used for stepsize
#'   adaptation, so warmup samples should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param adapt.delta An optional real value between 0 and 1, the thin of the jumps
#' in a HMC method. By default is 0.9.
#' @param  tree.depth An integer of the maximum  depth of the trees  evaluated
#' during each iteration. By default is 10.
#' @param ... Further arguments passed to  \code{rstan()} function.
#'
#' @details
#' This is the principal package's function and the link with \pkg{Stan}, this function
#' fits the posterior distribution of every parameter for a defined model using a
#' HMC method.
#'
#' Every estimated model become a \code{varstan} object, with different methods
#' for summary, diagnostic, forecast and plotting.
#'
#' \bold{Defining priors}
#'
#' Default priors are chosen to be non or very weakly informative so that their
#' influence on the results will. However, after getting more familiar with Bayesian
#' statistics, I recommend you to start thinking about reasonable informative priors
#' for your model parameters.
#'
#' Those can be changed using the function \code{set_prior()} before estimating the
#' model with the \code{varstan()} function. For checking the defined priors use
#' \code{get_prior()} and \code{report()} functions.
#'
#' \bold{Adjusting the sampling behavior of \pkg{Stan}}
#'
#' In addition to choosing the number of iterations, warmup samples, and
#' chains, users can control the behavior of the NUTS sampler, by using the
#' \code{control} argument. The most important reason to use \code{control}
#' is to decrease (or eliminate at best) the number of divergent transitions
#' that cause a bias in the obtained posterior samples. Whenever you see the
#' warning "There were x divergent transitions after warmup." you should really
#' think about increasing \code{adapt_delta}.  Increasing \code{adapt_delta} will
#' slow down the sampler but will decrease the number of divergent transitions
#' threatening the validity of your posterior samples.
#'
#' Another problem arises when the depth of the tree being evaluated in each iteration
#' is exceeded. This is less common than having divergent transitions, but may also
#' bias the posterior samples. When it happens, \pkg{Stan} will throw out a warning
#' suggesting to increase \code{max_treedepth}. For more details on the \code{control}
#' argument see \code{\link[rstan:stan]{stan}}.
#'
#' @author Asael Alonzo Matamoros
#'
#' @export
#' @importFrom stats frequency
#'
#' @return  a \code{varstan} object with the estimated time series model.
#'
#' @seealso \code{\link[rstan:stan]{rstan:stan}}.
#'
#' @references
#' Carpenter, B. and Gelman, A. and Hoffman, D. and  Lee, D. and Goodrich, B. and
#' Betancourt, M. and Brubaker, and Guo, L. and Riddell. 2017. Stan: A probabilistic
#' programming language. \emph{Journal of Statistical Software} 76(1).
#' \code{doi: 10.18637/jss.v076.i01}.
#'
#' Stan Development Team. (2018). Stan Modeling Language Users Guide and Reference Manual,
#' Version 2.18.0. \code{url: http://mc-stan.org}.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel
#' Models Using Stan. \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \code{doi:10.18637/jss.v080.i01}
#'
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}.
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  # Fitting a seasonal ARIMA model
#'  mod1 = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#'  fit1 = varstan(mod1,chains = 1)
#'  fit1
#'
#'  # Fitting a GARCH(1,1) model
#'  dat = garch(ipc,order = c(1,1,0))
#'  fit2 = varstan(dat,chains = 1)
#'  fit2
#' }
#'
varstan = function(model,chains = 4,iter = 2000,warmup = floor(iter/2),
                   adapt.delta = 0.90,tree.depth = 10,...){

  if(!is.model(model))
    stop(class(model),"is not an available current model in varstan")

  m = list()

  if(is.Sarima(model)) sft = fit_Sarima(model,chains,iter,warmup,adapt.delta,tree.depth)
  if(is.naive(model))  sft = fit_Sarima(model,chains,iter,warmup,adapt.delta,tree.depth)
  if(is.garch(model))  sft = fit_garch(model,chains,iter,warmup,adapt.delta,tree.depth)
  if(is.SVM(model))    sft = fit_SVM(model,chains,iter,warmup,adapt.delta,tree.depth)
  if(is.ssm(model))    sft = fit_ssm(model,chains,iter,warmup,adapt.delta,tree.depth)

  sp = list(Algorithm = "HMC NUTS",chains = chains,iter = iter,warmup = warmup,
            adapt.delta =adapt.delta,max_treedepth = tree.depth)

  m = list(stanfit = sft,
           stan.parmaters = sp,
           model = model,
           series.name = model$series.name,
           ts = model$yreal
           )

  attr(m,"class") = "varstan"

  return(m)
}
#' Checks if is a varstan object
#'
#' @method is varstan
#'
#' @param object a varstan object
#'
#' @noRd
#'
is.varstan = function(object){
  y = FALSE
  if(is(object,"varstan") ) y = TRUE
  return (y)
}
#' Extract chains of an stanfit object implemented in rstan package
#'
#' @usage extract_stan(object,pars,permuted = TRUE,inc_warmup = FALSE,include = TRUE, ...)
#'
#' @param object a varstan object
#' @param pars n optional character vector providing the parameter names
#' (or other quantity names) of interest. If not specified, all parameters
#' and other quantities are used. The log-posterior with name lp__ is also
#' included by default.
#' @param permuted A logical scalar indicating whether the draws
#' after the warmup period in each chain should be permuted and merged.
#' If FALSE, the original order is kept. For each stanfit object,
#' the permutation is fixed (i.e., extracting samples a second time will
#' give the same sequence of iterations).
#' @param inc_warmup A logical scalar indicating whether to include the
#' warmup draws. This argument is only relevant if permuted is FALSE.
#' @param include A logical scalar indicating whether the parameters named
#' in pars should be included (TRUE) or excluded (FALSE).
#' @param ... Further arguments passed to  \code{extract} function.
#'
#' @author Asael Alonzo Matamoros
#'
#' @importFrom rstan extract
#' @export
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  # Fitting a GARCH(1,1) model
#'  dat = garch(ipc,order = c(1,1,0))
#'  fit2 = varstan(dat,chains = 1)
#'
#'  # Extracting the mean parameter
#'  mu0 = extract_stan(fit2,pars = "mu0")
#' }
#'
extract_stan = function(object,pars,permuted = TRUE,
                        inc_warmup = FALSE,include = TRUE,...){

  if(!is.varstan(object) )
    stop("The current object is not a varstan class")

  chains = rstan::extract(object$stanfit,pars,permuted,inc_warmup,include,...)

  return(chains)
}
#' Convert to a stanfit object.
#'
#' Convert a \code{varstan} object to a \code{stanfit} object of the
#' \pkg{rstan} package.
#'
#' @usage  as.stan(object)
#'
#' @param object a varstan object.
#' @return  a stanfit object.
#'
#' @author  Asael Alonzo Matamoros
#'
#' @export
#'
#' @examples
#' \donttest{
#'  # Fitting a GARCH(1,1) model
#'  dat = garch(ipc,order = c(1,1,0))
#'  fit1 = varstan(dat,chains = 1)
#'
#'  # Converting to a Stanfit object
#'  stanfit1 = as.stan(fit1)
#' }
as.stan = function(object){
  if( !is.varstan(object) )
    stop("The current object is not a varstan class")

  stanfit = object$stanfit

  return(stanfit)
}
#' Extracts all the order coefficients in a list
#'
#' @param object A varstan object
#' @noRd
#'
get_order = function(object){
  if(!is.varstan(object))
    stop("The object is not a varstan class")

  if(is.Sarima(object$model)) return(get_order_arima(object$model))
  if(is.garch(object$model))  return(get_order_garch(object$model))
  if(is.SVM(object$model))    return(get_order_garch(object$model))
  if(is.naive(object$model))  return(get_order_arima(object$model))
}
#' Max order  coefficients in a varma model
#'
#' @param object A varstan object
#' @noRd
#'
max_order = function(object){
  if(!is.varstan(object))
    stop("The object is not a varstan class")

  if(is.Sarima(object$model)) return(max_order_arima(object$model))
  if(is.garch(object$model))  return(max_order_garch(object$model))
  if(is.SVM(object$model))    return(max_order_garch(object$model))
  if(is.naive(object$model))  return(max_order_arima(object$model))
}
#' Extracts all the order coefficients in a list
#'
#' @param object A varstan object
#' @noRd
#'
Total_order = function(object){
  if(!is.varstan(object))
    stop("The object is not a varstan class")

  order = as.numeric(get_order(object))
  if(is.Sarima(object$model)){
    n = length(order)-1
    s1 = sum(order[1:n])+2
  }
  else s1 = sum(order)+2

  return(s1)
}
