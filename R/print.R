#' Print a Sarima model
#'
#' @param  x a Sarima model from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object
#' @method print Sarima
#' @export
#'
print.Sarima = function(x,...){
  if(!is.Sarima(x))
    stop("The current object is not an arima model")
  report(x)
}
#' Print a naive model
#'
#' @param  x a naive model from the varstan package
#'
#' @return None. prints the object
#'
#' @method print naive
#' @param ... additional values need in print methods
#' @export
#'
print.naive = function(x,...){
  if(!is.naive(x))
    stop("The current object is not a naive model")
  report(x)
}
#' Print a garch model
#'
#' @param x a garch model from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object
#'
#' @method print garch
#' @export
#'
print.garch = function(x,...){
  if(!is.garch(x))
    stop("The current object is not a garch model")
  report(x)
}
#' Print a Stochastic Volatility model
#'
#' @param x a SVM model from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object
#'
#' @method print SVM
#' @export
#'
print.SVM = function(x,...){
  if(!is.SVM(x))
    stop("The current object is not a SVM model")
  report(x)
}
#' Print a state-space model
#'
#' @param x a ssm model from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object.
#'
#' @method print ssm
#' @export
#'
print.ssm = function(x,...){
  if(!is.ssm(x))
    stop("The current object is not a ssm model")
  report(x)
}
#' Print a Local Level model
#'
#' @param x a LocalLevel model from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object.
#'
#' @method print LocalLevel
#' @export
#'
print.LocalLevel = function(x,...){
  if(!is.LocalLevel(x))
    stop("The current object is not a ssm model")
  report(x)
}
#' Print a Holt model
#'
#' @param x a Holt model from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object.
#'
#' @method print Holt
#' @export
#'
print.Holt = function(x,...){
  if(!is.Holt(x))
    stop("The current object is not a ssm model")
  report(x)
}
#' Print a Holt-Winter model
#'
#' @param x a Hw from the varstan package
#' @param ... additional values need in print methods
#'
#' @return None. prints the object.
#'
#' @method print Hw
#' @export
#'
print.Hw = function(x,...){
  if(!is.Hw(x))
    stop("The current object is not a ssm model")
  report(x)
}
#' Print a varstan object
#'
#' @param x a varstan object
#' @param ... additional values need in print methods
#'
#' @return None. prints the object.
#'
#' @method print varstan
#' @export
#'
print.varstan = function(x,...){
  if(!is.varstan(x))
    stop("The current object is not a varstan model")

  model(x)
  print(summary(x))
  cat("\n Samples were drawn using sampling(NUTS). For each parameter, ess")
  cat("\n is the effective sample size, and Rhat is the potential")
  cat("\n scale reduction factor on split chains (at convergence, Rhat = 1). \n")
}
