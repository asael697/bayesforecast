#' Get parameters of a varstan object
#'
#' Get the sampled parameters of a varstan object.
#'
#' @param object a varstan object
#'
#' @return a vector with the sampled parameters
#'
#' @author Asael Alonzo Matamoros
#' @export
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  sf1 = auto.sarima(birth)
#'  get_parameters(sf1)
#' }
#'
get_parameters = function(object){
  if(!is.varstan(object))
    stop("The current object is not a varstan class")

  if(is.Sarima(object$model)) gp = get_params_arima(object$model)
  if(is.naive(object$model))  gp = get_params_arima(object$model)
  if(is.garch(object$model))  gp = get_params_garch(object$model)
  if(is.SVM(object$model))    gp = get_params_svm(object$model)
  if(is.ssm(object$model))    gp = get_params_ssm(object$model)
  if(is.LocalLevel(object$model))gp = get_params_ssm(object$model)
  if(is.Holt(object$model))   gp = get_params_ssm(object$model)
  if(is.Hw(object$model))     gp = get_params_ssm(object$model)

  return(gp$include)
}
#'
#' Excluded parameters in a Sarima model
#'
#' @param dat a Sarima model
#' @noRd
#'
get_params_arima = function(dat){
  include = c("mu0","sigma0")
  if(dat$p > 0) include = c(include,"ar")
  if(dat$q > 0) include = c(include,"ma")
  if(dat$P > 0) include = c(include,"sar")
  if(dat$Q > 0) include = c(include,"sma")
  if(dat$d1 >0) include = c(include,"breg")
  exclude = c("ar0","ma0","sar0","sma0")
  pars = list(include = c(include,"loglik"),exclude = exclude)
  return(pars)
}
#' Excluded parameters in a  garch model
#'
#' @param dat a garch model
#' @noRd
#'
get_params_garch = function(dat){
  include = c("mu0","sigma0")
  if(dat$s > 0) include = c(include,"arch")
  if(dat$k > 0) include = c(include,"garch")
  if(dat$h > 0) include = c(include,"mgarch")
  if(dat$p > 0) include = c(include,"ar")
  if(dat$q > 0) include = c(include,"ma")
  if(dat$d1 >0) include = c(include,"breg")
  if(dat$genT)  include = c(include,"v")
  if(dat$asym1) include = c(include,"gamma")

  exclude = c("ar0","ma0")
  pars = list(include = c(include,"loglik"),exclude = exclude)
  return(pars)
}
#' Excluded parameters in a SVM model
#'
#' @param dat a SVM model
#' @noRd
#'
get_params_svm = function(dat){
  include = c("mu0","sigma0")
  if(dat$s > 0) include = c(include,"alpha")
  if(dat$k > 0) include = c(include,"beta")
  if(dat$h > 0) include = c(include,"mgarch")
  if(dat$p > 0) include = c(include,"ar")
  if(dat$q > 0) include = c(include,"ma")
  if(dat$d1 >0) include = c(include,"breg")
  if(dat$genT)  include = c(include,"v")
  if(dat$asym1) include = c(include,"gamma")

  exclude = c("ar0","ma0")
  pars = list(include = c(include,"loglik"),exclude = exclude)
  return(pars)
}
#' Excluded parameters in a SSM model
#'
#' @param dat a SSM model
#' @noRd
#'
get_params_ssm = function(dat){
  include = c("sigma0","level")
  if(dat$is_td) include = c(include,"trend")
  if(dat$is_dp) include = c(include,"damped")
  if(dat$is_ss) include = c(include,"seasonal")
  include = c(include,"level1")
  if(dat$is_td) include = c(include,"trend1")
  if(dat$is_ss) include = c(include,"seasonal1")
  if(dat$d1 >0) include = c(include,"breg")
  if(dat$genT)  include = c(include,"v")

  exclude = c("ar0","ma0")
  pars = list(include = c(include,"loglik"),exclude = exclude)
  return(pars)
}
