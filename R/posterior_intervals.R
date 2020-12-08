#' Posterior uncertainty intervals
#'
#' The \code{posterior_interval} function computes Bayesian posterior uncertainty
#' intervals. These intervals are often referred to as \emph{credible}
#' intervals, for more details see \pkg{rstanarm}
#'
#' @aliases posterior_interval
#'
#' @author Asael Alonzo Matamoros
#'
#' @param mat a matrix containing the posterior samples of a fitted parameter
#' @param ... Further arguments passed to  \code{posterior_intervals}.
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \code{90\%} intervals (\code{prob=0.9}) rather than the traditionally used
#'   \code{95\%}.
#'
#' @return A matrix with two columns and as many rows as model parameters (or
#'   the subset of parameters specified by \code{pars} and/or
#'   \code{regex_pars}). For a given value of \code{prob}, \eqn{p}, the columns
#'   correspond to the lower and upper \code{100*p\%} interval limits and have the
#'   names \eqn{100\alpha/2} and \eqn{100(1 - \alpha/2)}\code{\%}, where \eqn{\alpha
#'   = 1-p}. For example, if \code{prob=0.9} is specified (a \code{90\%}
#'   interval), then the column names will be \code{"5\%"} and \code{"95\%"},
#'   respectively.
#'
#' @export
#'
posterior_interval = function(mat, prob = 0.90, ...){
  if(is.data.frame(mat))
    mat1 = as.matrix(mat)
  else mat1 = mat
  rstantools::posterior_interval(object = mat1,prob = prob)
}
