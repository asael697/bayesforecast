#' Print a full report of the time series model in a  varstan object.
#'
#' The function returns a report with the users defined model for the given time series data
#' and all the current defined priors of the model.
#'
#' @usage  report(object,...)
#'
#' @aliases report report.varstan report.Sarima report.garch report.varma report.Bekk report.naive
#'
#' @param object an object varstan object or one of the defined current defined reports in varstan package
#' @param ... additional values need in print methods
#'
#' @details if \code{object} is a varstan object the function will print the information of the
#' defined model inside of the object. If \code{object} is one of the model classes (like Sarima or garch)
#' then it will print the report information as well.
#'
#' @author Asael Alonzo Matamoros
#'
#' @export
#'
#' @return none. prints a string with the defined time series model report
#'
#' @examples
#' library(astsa)
#' dat2 = garch(birth,order = c(1,1,0))
#' report(dat2)
#'
report <- function(object,...) {
  UseMethod("report")
}
#' @aliases report
#' @method report varstan
#' @export
#'
report.varstan = function(object,...){
  if(!is.varstan(object))
    stop("The current object is not a varstan class")

  if( is.Sarima(object$model)) report.Sarima(object$model)
  if( is.naive(object$model))  report.naive(object$model)
  if( is.garch(object$model))  report.garch(object$model)
  if( is.SVM(object$model))    report.SVM(object$model)
  if( is.ssm(object$model))    report.ssm(object$model)
  if( is.LocalLevel(object$model))report.LocalLevel(object$model)
  if( is.Holt(object$model))   report.Holt(object$model)
  if( is.Hw(object$model))     report.Hw(object$model)
}
#' @aliases report
#' @method report Sarima
#' @export
#'
report.Sarima = function(object,...){
  if( !is.Sarima(object))
    stop("The object is not a Sarima model \n")

  model.Sarima(object)
  cat("Priors: \n Intercept:\n")
  get_prior(object,par = "mu0")
  cat("\n Scale Parameter: \n")
  get_prior(object,par = "sigma0")
  cat("\n")
  if(object$p  > 0 )get_prior(model = object,par = "ar")
  if(object$q  > 0 )get_prior(model = object,par = "ma")
  if(object$P != 0 || object$Q != 0 || object$D != 0){
    cat("\n Seasonal Parameters: \n")
    if(object$P  > 0 )get_prior(model = object,par = "sar")
    if(object$Q  > 0 )get_prior(model = object,par = "sma")
  }
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
}
#' @aliases report
#' @method report naive
#' @export
#'
report.naive = function(object,...){
  if( !is.naive(object))
    stop("The object is not a naive model \n")

  model.naive(object)
  cat("Priors: \n Intercept:\n")
  get_prior(object,par = "mu0")
  cat("\n Scale Parameter: \n")
  get_prior(object,par = "sigma0")
  cat("\n")
  if(object$D > 0 ){
    cat("\n period:",object$period,"\n")
  }
}
#' @aliases report
#' @method report garch
#' @export
#'
report.garch = function(object,...){
  if(!is.garch(object))
    stop("The object is not a garch model \n")

  model.garch(object)
  cat("Priors: \n Intercept:\n")
  get_prior(object,par = "mu0")
  cat("\n Scale Parameter: \n")
  get_prior(object,par = "sigma0")
  cat("\n")
  if(object$s  > 0 )get_prior(model = object,par = "arch")
  if(object$k  > 0 )get_prior(model = object,par = "garch")
  if(object$h  > 0 )get_prior(model = object,par = "mgarch")
  if(object$p != 0 || object$q != 0 ){
    cat("\n mean Parameters: \n")
    if(object$p  > 0 )get_prior(model = object,par = "ar")
    if(object$q  > 0 )get_prior(model = object,par = "ma")
  }
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
  if(object$genT){
    cat("\n Generalized t-student \n")
    cat("\n lambda ~ G(v/2,v/2) \n")
    get_prior(model = object,par = "dfv")
  }
  if(object$asym1){
    if(identical(object$asym,1)) cat("\n logistic asymmetric \n")
    if(identical(object$asym,2)) cat("\n exponential asymmetric \n")
    get_prior(model = object,par = "gamma")
  }
}
#' @aliases report
#' @method report SVM
#' @export
#'
report.SVM = function(object,...){
  if(!is.SVM(object))
    stop("The object is not a stochastic volatility model \n")

  model.SVM(object)
  cat("Priors: \n Intercept:\n")
  get_prior(object,par = "mu0")
  cat("\n log Scale Parameter: \n")
  get_prior(object,par = "sigma0")
  cat("\n")
  get_prior(model = object,par = "alpha")
  get_prior(model = object,par = "beta")
  if(object$p != 0 || object$q != 0 ){
    cat("\n mean Parameters: \n")
    if(object$p  > 0 )get_prior(model = object,par = "ar")
    if(object$q  > 0 )get_prior(model = object,par = "ma")
  }
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
}
#' @aliases report
#' @method report ssm
#' @export
#'
report.ssm = function(object,...){
  if(!is.ssm(object))
    stop("The object is not a ssm model \n")

  model.ssm(object)
  cat("Priors: \n Scale Parameter:\n")
  get_prior(object,par = "sigma0")
  cat("\n")
  get_prior(object,par = "level")
  if(object$is_td)get_prior(model = object,par = "trend")
  if(object$is_dp)get_prior(model = object,par = "damped")
  if(object$is_ss)get_prior(model = object,par = "seasonal")
  cat("Initial values \n")
  get_prior(object,par = "level1")
  if(object$is_td)get_prior(model = object,par = "trend1")
  if(object$is_ss)get_prior(model = object,par = "seasonal1")
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
  if(object$genT){
    cat("\n Generalized t-student \n")
    cat("\n lambda ~ G(v/2,v/2) \n")
    get_prior(model = object,par = "dfv")
  }
}
#' @aliases report
#' @method report LocalLevel
#' @export
#'
report.LocalLevel = function(object,...){
  if(!is.LocalLevel(object))
    stop("The object is not a Local level model \n")

  model.LocalLevel(object)
  cat("Priors: \n Scale Parameter:\n")
  get_prior(object,par = "sigma0")
  cat("\n")
  get_prior(object,par = "level")
  cat("Initial values \n")
  get_prior(object,par = "level1")
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
  if(object$genT){
    cat("\n Generalized t-student \n")
    cat("\n lambda ~ G(v/2,v/2) \n")
    get_prior(model = object,par = "dfv")
  }
}
#' @aliases report
#' @method report Holt
#' @export
#'
report.Holt = function(object,...){
  if(!is.Holt(object))
    stop("The object is not a Holt model \n")

  model.Holt(object)
  cat("Priors: \n Scale Parameter:\n")
  get_prior(object,par = "sigma0")
  cat("\n")
  get_prior(object,par = "level")
  if(object$is_td)get_prior(model = object,par = "trend")
  if(object$is_dp)get_prior(model = object,par = "damped")
  cat("Initial values \n")
  get_prior(object,par = "level1")
  if(object$is_td)get_prior(model = object,par = "trend1")
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
  if(object$genT){
    cat("\n Generalized t-student \n")
    cat("\n lambda ~ G(v/2,v/2) \n")
    get_prior(model = object,par = "dfv")
  }
}
#' @aliases report
#' @method report Hw
#' @export
#'
report.Hw = function(object,...){
  if(!is.Hw(object))
    stop("The object is not a Holt Winters model \n")

  model.Hw(object)
  cat("Priors: \n Scale Parameter:\n")
  get_prior(object,par = "sigma0")
  cat("\n")
  get_prior(object,par = "level")
  if(object$is_td)get_prior(model = object,par = "trend")
  if(object$is_dp)get_prior(model = object,par = "damped")
  if(object$is_ss)get_prior(model = object,par = "seasonal")
  cat("Initial values \n")
  get_prior(object,par = "level1")
  if(object$is_td)get_prior(model = object,par = "trend1")
  if(object$is_ss)get_prior(model = object,par = "seasonal1")
  if(object$d1 > 0 ){
    cat("\n Regression Parameters: \n")
    get_prior(model = object,par = "breg")
  }
  if(object$genT){
    cat("\n Generalized t-student \n")
    cat("\n lambda ~ G(v/2,v/2) \n")
    get_prior(model = object,par = "dfv")
  }
}
#' Generic function for extracting information about prior distributions
#'
#' The function returns a report with the users defined model for the given time series data
#' and all the current defined priors of the model.
#'
#' @aliases prior_summary
#'
#' @param object a varstan object or one of the defined current defined reports in varstan package
#' @param ... additional values need in print methods
#'
#' @details if \code{object} is a varstan object the function will print the information of the
#' defined model inside of the object. If \code{object} is one of the model classes (like Sarima or garch)
#' then it will print the report information as well.
#'
#' @return none. prints a string with the defined time series model report
#'
#' @author Asael Alonzo Matamoros
#'
#' @method prior_summary varstan
#' @importFrom rstantools prior_summary
#' @export
#' @export prior_summary
#'
#' @examples
#' library(astsa)
#' dat2 = garch(birth,order = c(1,1,0))
#' prior_summary(dat2)
#'
prior_summary.varstan = function(object, ...){
  if(!is.varstan(object))
    stop("The current object is not a varstan class")
  report(object)
}
