#' Fit an arima model
#'
#' Fit a arima(p,d,q) model  in STAN
#'
#' The function returns a list with the fitted model in stan
#'
#' @usage  fit.arima(model,chains,iter,warmup,...)
#'
#' @param model A time series object for the varstan models
#' @param chains the number of chains to be run
#' @param iter the number of iteration per chain
#' @param warmup the number of initial iteration to be burned
#' @param adapt.delta the thin of the jumps in a HMC method
#' @param  tree.depth maximum  depth of the trees  evaluated
#' during each iteration
#'
#' @import rstan
#'
#' @author  Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return  a stanfit object
#'
fit_Sarima = function(model,chains=4,iter=2000,warmup=floor(iter/2),adapt.delta = 0.90,tree.depth,...){
  stanfit = rstan::sampling(stanmodels$Sarima,
                            data = model,
                            chains = chains,
                            iter = iter,
                            warmup = warmup,
                            control = list(adapt_delta = adapt.delta,max_treedepth = tree.depth))

  return(stanfit)
}
#' Fit a mgarch model
#'
#' Fit a garch(s,k,h) model  in STAN
#'
#' The function returns a list with the fitted model in stan
#'
#' @usage  fit.garch(model,chains,iter,warmup,...)
#'
#' @param model A time series object for the varstan models
#' @param chains the number of chains to be run
#' @param iter the number of iteration per chain
#' @param warmup the number of initial iteration to be burned
#' @param adapt.delta the thin of the jumps in a HMC method
#' @param  tree.depth maximum  depth of the trees  evaluated
#' during each iteration
#'
#' @import rstan
#'
#' @author  Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return  a stanfit object
#'
fit_garch = function(model,chains=4,iter=2000,warmup=floor(iter/2),adapt.delta = 0.90,tree.depth,...){

  pars = get_params_garch(model)$exclude
  stanfit = rstan::sampling(stanmodels$tgarch,
                            data = model,
                            chains = chains,
                            iter = iter,
                            warmup = warmup,
                            control = list(adapt_delta = adapt.delta,max_treedepth = tree.depth))

  return(stanfit)
}
#' Fit a SVM model
#'
#' Fit a SVM(ts, arma = c(p,q),xreg = NULL) model  in STAN
#'
#' The function returns a list with the fitted model in stan
#'
#' @usage  fit_SVM(model,chains,iter,warmup,...)
#'
#' @param model A time series object for the varstan models
#' @param chains the number of chains to be run
#' @param iter the number of iteration per chain
#' @param warmup the number of initial iteration to be burned
#' @param adapt.delta the thin of the jumps in a HMC method
#' @param  tree.depth maximum  depth of the trees  evaluated
#' during each iteration
#'
#' @import rstan
#'
#' @author Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return  a stanfit object
#'
fit_SVM = function(model,chains=4,iter=2000,warmup=floor(iter/2),adapt.delta = 0.90,tree.depth,...){

  pars = get_params_garch(model)$exclude


  stanfit = rstan::sampling(stanmodels$SVM,
                              data = model,
                              chains = chains,
                              iter = iter,
                              warmup = warmup,
                              control = list(adapt_delta = adapt.delta,max_treedepth = tree.depth))

  return(stanfit)
}
#' Fit a ssm model
#'
#' Fit a ssm(ts, arma = c(p,q),xreg = NULL) model  in STAN
#'
#' The function returns a list with the fitted model in stan
#'
#' @param model A time series object for the varstan models
#' @param chains the number of chains to be run
#' @param iter the number of iteration per chain
#' @param warmup the number of initial iteration to be burned
#' @param adapt.delta the thin of the jumps in a HMC method
#' @param  tree.depth maximum  depth of the trees  evaluated
#' during each iteration
#'
#' @import rstan
#'
#' @author Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return  a stanfit object
#'
fit_ssm = function(model,chains=4,iter=2000,warmup=floor(iter/2),adapt.delta = 0.90,tree.depth,...){

  pars = get_params_ssm(model)$exclude


  stanfit = rstan::sampling(stanmodels$ets,
                            data = model,
                            chains = chains,
                            iter = iter,
                            warmup = warmup,
                            control = list(adapt_delta = adapt.delta,max_treedepth = tree.depth))

  return(stanfit)
}
