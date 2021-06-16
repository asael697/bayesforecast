#' Out-of-sample predictive errors
#'
#' This is a convenience function for computing \eqn{y - y_{h}}{y - yh}
#' The method for stanreg objects  calls \code{\link{posterior_predict}}
#' internally, where as the method accepts the data.frame returned by
#' \code{posterior_predict} as input and can be used to avoid multiple calls to
#' \code{posterior_predict}.
#'
#' @aliases predictive_error
#'
#' @param object Either a fitted model object returned by one of the \pkg{rstanarm}
#' modeling functions (a stanreg object) or, for the \code{"ppd"} method, a matrix
#' of draws from the posterior predictive distribution returned by \code{\link{posterior_predict}}.
#' @param xreg Optional, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param newdata An array with the newdata vector.
#' @param draws,seed Optional arguments passed to \code{\link{posterior_predict}}.
#' Please see the \strong{Note} section below if \code{newdata} will be specified.
#' @param ... Further arguments passed to  \code{predictive_error}.
#'
#'
#' @return
#' A \code{draws} by \code{nrow(newdata)} data.frame.
#'
#' @note
#' If \code{object} is a \strong{varstan} object of a varma model then newdata has to be a matrix
#' with number of \strong{cols} as the dimension of the time series and number of \strong{rows}
#' as the number new elements.
#'
#' @note If \code{object} is a \code{posterior_predict} data.frame, then the
#' length of \code{newdata} has to be equal to the \code{ncol} of \code{object}.
#'
#' @note If \code{object} is a \code{posterior_predict} data.frame, for a \strong{varma} model,
#' then the dimension product of \code{newdata} matrix has to be equal to
#' the \code{ncol} of \code{object}.
#'
#' @seealso \code{posterior_predict} function from rstanarm package, to draw
#'   from the posterior predictive distribution without computing predictive
#'   errors.
#'
#' @importFrom rstantools predictive_error
#' @method predictive_error varstan
#' @export
#' @export predictive_error
#'
predictive_error.varstan = function(object,newdata = NULL,xreg = NULL,draws = 1000,
                                    seed = NULL,...){

  if( !is.varstan(object) )
     stop("The current object is not a varstan object")

  if(is.null(newdata))
    return(as.data.frame(extract_stan(object,"residuals", permuted = TRUE) ))

  if (!is.numeric(newdata) | !stats::is.ts(newdata))
    stop("The current object is not a time series object")

  yh = as.matrix(posterior_predict(object = object,h = length(newdata),xreg = xreg,draws = draws,seed = seed) )


  return( rstantools::predictive_error(object = yh,y = newdata) )
}
#'
#'@aliases predictive_error.varstan
#'@export
#'
predictive_error.data.frame = function(object,newdata = NULL,xreg = NULL,draws = 1000,
                                       seed = NULL,...){

  if(!is.data.frame(object))
    stop("The current object is not a data.frame object")

  if (!is.numeric(newdata) | !stats::is.ts(newdata))
    stop("The current newdata is not a time series object")

  obj = as.matrix(object)

  if(ncol(obj) != length(newdata))
    stop("The dimensions not match")

  return( rstantools::predictive_error(object = obj,y = newdata) )
}
#' Generic function and method for extract the residual of a varstan object
#'
#' The function returns the posterior estimate of the residuals
#' of a varstan model, similar to the residual functions of other
#' packages.
#'
#' @param object A varstan object, \code{\link[=varstan]{varstan}}
#' @param robust A boolean value, if its \code{TRUE} it returns the median of the posterior distribution,
#' and if its \code{FALSE} it returns the mean, by default is the \code{FALSE} value
#' @param ... Further arguments passed to  \code{posterior_residual}.
#'
#' @details
#' This function only extracts the point-wise estimate of the time series residuals
#' for extracting all the data use \code{extract_stan()} or posterior_intervals function
#'
#' @author Asael Alonzo Matamoros
#'
#' @export
#' @importFrom stats median ts
#'
#' @return
#' An array with the posterior estimate of the residuals of the time series model,
#'
#' @method residuals varstan
#' @export
#'
residuals.varstan = function(object,robust = FALSE,...){
  if( !is.varstan(object) )
    stop("The current object is not a varstan class")

  post = as.data.frame(extract_stan(object,"residuals", permuted = TRUE) )
  if(robust) sum1 = t(matrix(apply(post,2,stats::median),nrow = 1,byrow = TRUE))
  else sum1 = t(matrix(apply(post,2,mean),nrow = 1,byrow = TRUE))

  resd = stats::ts(sum1,start = stats::start(object$ts) ,frequency = stats::frequency(object$ts))

  return(resd)
}
