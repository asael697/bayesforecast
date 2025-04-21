#' Fit an ARIMA model.
#'
#' Fit an `arima(p,d,q)` model in STAN.
#'
#' The function returns a list with the fitted model in Stan.
#'
#' @param model a time series object for the `varstan` models.
#' @param chains the number of chains to be run.
#' @param iter the number of iteration per chain.
#' @param warmup the number of initial iteration to be burned.
#' @param adapt.delta the thin of the jumps in a HMC method.
#' @param tree.depth maximum tree depth per iteration.
#'
#' @import rstan
#'
#' @author Asael Alonzo Matamoros.
#'
#' @noRd
#'
#' @return a stanfit object
#'
fit_Sarima = function(model, chains = 4, iter = 2000, warmup = floor(iter/2),
                      adapt.delta = 0.90, tree.depth, ...){

  stanfit = rstan::sampling(stanmodels$Sarima,
                            data = model,
                            chains = chains,
                            iter = iter,
                            warmup = warmup,
                            control = list(adapt_delta = adapt.delta,
                                           max_treedepth = tree.depth))

  return(stanfit)
}
#' Fit a MGARCH model.
#'
#' Fit a `garch(s,k,h)` model in Stan.
#'
#' The function returns a list with the fitted model in Stan.
#'
#' @param model a time series object for the `varstan` models.
#' @param chains the number of chains to be run.
#' @param iter the number of iteration per chain.
#' @param warmup the number of initial iteration to be burned.
#' @param adapt.delta the thin of the jumps in a HMC method.
#' @param  tree.depth maximum  tree depth per iteration.
#'
#' @import rstan
#'
#' @author Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return a stanfit object
#'
fit_garch = function(model, chains = 4, iter = 2000, warmup = floor(iter/2),
                     adapt.delta = 0.90, tree.depth, ...){

  pars = get_params_garch(model)$exclude
  stanfit = rstan::sampling(stanmodels$tgarch,
                            data = model,
                            chains = chains,
                            iter = iter,
                            warmup = warmup,
                            control = list(adapt_delta = adapt.delta,
                                           max_treedepth = tree.depth))

  return(stanfit)
}

#' Fit a SVM model
#'
#' Fit a `SVM(ts, arma = c(p,q), xreg = NULL)` model  in `Stan`.
#'
#' The function returns a list with the fitted model in `Stan`.
#'
#' @param model A time series object for the `varstan` models.
#' @param chains the number of chains to be run.
#' @param iter the number of iteration per chain.
#' @param warmup the number of initial iteration to be burned.
#' @param adapt.delta the thin of the jumps in a HMC method.
#' @param  tree.depth maximum tree depth per iteration.
#'
#' @import rstan
#'
#' @author Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return a `stanfit` object
#'
fit_SVM = function(model, chains = 4, iter = 2000, warmup = floor(iter/2),
                   adapt.delta = 0.90, tree.depth, ...){

  pars = get_params_garch(model)$exclude


  stanfit = rstan::sampling(stanmodels$SVM,
                              data = model,
                              chains = chains,
                              iter = iter,
                              warmup = warmup,
                              control = list(adapt_delta = adapt.delta,
                                             max_treedepth = tree.depth))

  return(stanfit)
}
#' Fit a SSM model.
#'
#' Fits a `ssm(ts, arma = c(p,q), xreg = NULL)` model  in `Stan`.
#'
#' The function returns a list with the fitted model in `Stan`.
#'
#' @param model A time series object for the `varstan` models.
#' @param chains the number of chains to be run.
#' @param iter the number of iteration per chain.
#' @param warmup the number of initial iteration to be burned.
#' @param adapt.delta the thin of the jumps in a HMC method.
#' @param  tree.depth maximum tree depth per iteraiton.
#'
#' @import rstan
#'
#' @author Asael Alonzo Matamoros
#'
#' @noRd
#'
#' @return a `stanfit` object.
#'
fit_ssm = function(model, chains = 4, iter = 2000, warmup = floor(iter/2),
                   adapt.delta = 0.90, tree.depth, ...){

  pars = get_params_ssm(model)$exclude

  stanfit = rstan::sampling(stanmodels$ets,
                            data = model,
                            chains = chains,
                            iter = iter,
                            warmup = warmup,
                            control = list(adapt_delta = adapt.delta,
                                           max_treedepth = tree.depth))

  return(stanfit)
}
