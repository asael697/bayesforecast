#' Print the defined model of a `varstan` object.
#'
#' The function returns a string with the users defined model for the given time series data.
#'
#' @param object a `varstan` object or one of the defined current defined models.
#' @param ... additional values need in print methods.
#'
#' @details
#' if \code{object} is a `varstan` object the function will print the information of the
#' defined model inside of the object. If \code{object} is one of the model classes (like
#' \code{Sarima}, \code{garch}, \code{SVM} or \code{varma}), then it will print the model
#' information as well.
#'
#' For full information of the model with the used priors use the function report or just
#' print the object.
#'
#' @author Asael Alonzo Matamoros.
#'
#' @seealso \code{report} \code{print}
#'
#' @aliases model model.varstan model.Sarima model.garch model.varma model.Bekk model.SVM
#'
#' @export
#' @return a  string with the defined time series model.
#'
#' @examples
#' model1 = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
#' model(model1)
#'
model = function(object,...) {
  UseMethod("model")
}
#'
#' @method  model varstan
#' @export model
#' @export
#'
model.varstan = function(object,...){
  if(!is.varstan(object))
    stop("The current object is not a varstan class")

  if( is.Sarima(object$model)) model.Sarima(object$model)
  if( is.naive(object$model))  model.naive(object$model)
  if( is.garch(object$model))  model.garch(object$model)
  if( is.SVM(object$model))    model.SVM(object$model)
  if( is.ssm(object$model))    model.ssm(object$model)
  if( is.LocalLevel(object$model))model.LocalLevel(object$model)
  if( is.Holt(object$model))   model.Holt(object$model)
  if( is.Hw(object$model))     model.Hw(object$model)
}
#'
#' @method model Sarima
#' @export model
#' @export
#'
model.Sarima = function(object,...){
  if( !is.Sarima(object))
    stop("The object is not a Sarima model \n")

  cat("\n")
  log = paste0("y ~ Sarima(",object$p,",",object$d,",",object$q,")")

  if(object$P != 0 || object$Q != 0 || object$D != 0)
    log = paste0(log,"(",object$P,",",object$D,",",object$Q,")[",object$period,"]")

  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")

  cat(log,"\n")
  cat(object$n,"observations and 1 dimension \n")
  cat("Differences:",object$d,"seasonal Differences:",object$D,"\n")
  cat("Current observations:",object$n1,"\n \n")
}
#'
#' @method  model naive
#' @export model
#' @export
#'
model.naive = function(object,...){
  if( !is.naive(object))
    stop("The object is not a naive model \n")

  cat("\n")
  log = paste0("y ~ Sarima(",object$p,",",object$d,",",object$q,")")

  if(object$d == 0)
    log = paste0("y ~ Random Walk()")
  else
    log = paste0("y ~ Random Walk(",object$period,")")

  cat(log,"\n")
  cat(object$n,"observations and 1 dimension \n")
  cat("Differences:",object$d,"seasonal Diferences:",object$D,"\n")
  cat("Current observations:",object$n1,"\n \n")
}
#'
#' @method  model garch
#' @export model
#' @export
#'
model.garch = function(object,...){
  if(!is.garch(object))
    stop("The object is not a garch model \n")

  cat("\n")
  if(object$p != 0 || object$q != 0){
    log = paste0("y ~ arma(",object$p,",",object$q,")")
    if(object$asym1)
      log = paste0(log,"+asymmetric garch(",object$s,",",object$k,",",object$h,")")
    else
      log = paste0(log,"+garch(",object$s,",",object$k,",",object$h,")")
  }
  else{
    if(object$asym1)
      log = paste0("y ~ asymmetric garch(",object$s,",",object$k,",",object$h,")")
    else
      log = paste0("y ~ garch(",object$s,",",object$k,",",object$h,")")
  }
  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")
  cat(log,"\n")
  if(object$genT) cat("Generalized t-student model \n")
  if(object$asym1){
    if(identical(object$asym,1)) cat("asymmetric logistic model \n")
    if(identical(object$asym,2)) cat("asymmetric exponential model \n")
  }
  cat(object$n,"observations and 1 dimension \n \n")
}
#' @method  model SVM
#' @export model
#' @export
#'
model.SVM = function(object,...){
  if(!is.SVM(object))
    stop("The object is not a garch model \n")

  cat("\n")
  if(object$p != 0 || object$q != 0){
    log = paste0("y ~ arma(",object$p,",",object$q,")")
    log = paste0(log,"+SVM")
  }
  else log = paste0("y ~ SVM")

  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")
  cat(log,"\n")
  cat(object$n,"observations and 1 dimension \n \n")
}
#'
#' @method  model ssm
#' @export model
#' @export
#'
model.ssm = function(object,...){
  if(!is.ssm(object))
    stop("The object is not a ssm model \n")

  cat("\n")
  log = "Local level"
  if(object$is_td){
    if(object$is_dp)log = paste0(log,"-damped trend")
    else log = paste0(log,"-trend")
  }
  if(object$is_ss)log = paste0(log,"-seasonal")
  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")
  cat(log,"\n")
  cat("state-space model \n")
  if(object$genT) cat("Generalized t-student model \n")
  cat(object$n,"observations and 1 dimension \n \n")
}
#'
#' @method model LocalLevel
#' @export model
#' @export
#'
model.LocalLevel = function(object,...){
  if(!is.LocalLevel(object))
    stop("The object is not a Local level model \n")

  cat("\n")
  log = "Local level"

  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")
  cat(log,"\n")
  cat("state-space model \n")
  if(object$genT) cat("Generalized t-student model \n")
  cat(object$n,"observations and 1 dimension \n \n")
}
#'
#' @method model Holt
#' @export model
#' @export
#'
model.Holt = function(object,...){
  if(!is.Holt(object))
    stop("The object is not a Holt model \n")

  cat("\n")
  log = "Holt"
  if(object$is_td){
    if(object$is_dp)log = paste0(log,"-damped trend")
  }
  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")
  cat(log,"\n")
  cat("state-space model \n")
  if(object$genT) cat("Generalized t-student model \n")
  cat(object$n,"observations and 1 dimension \n \n")
}
#'
#' @method model Hw
#' @export model
#' @export
#'
model.Hw = function(object,...){
  if(!is.Hw(object))
    stop("The object is not a Holt-Winters model \n")

  cat("\n")
  log = "Holt-Winters"
  if(object$is_td){
    if(object$is_dp)log = paste0(log,"-damped trend")
  }
  if(object$d1 > 0) log = paste0(log,".reg[",object$d1,"]")
  cat(log,"\n")
  cat("state-space model \n")
  if(object$genT) cat("Generalized t-student model \n")
  cat(object$n,"observations and 1 dimension \n \n")
}
