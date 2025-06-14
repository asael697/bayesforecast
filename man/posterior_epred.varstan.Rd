% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_predict.R
\name{posterior_epred.varstan}
\alias{posterior_epred.varstan}
\alias{posterior_epred}
\title{Expected Values of the Posterior Predictive Distribution}
\usage{
\method{posterior_epred}{varstan}(
  object,
  h = 0,
  xreg = NULL,
  robust = FALSE,
  draws = 1000,
  seed = NULL,
  ...
)
}
\arguments{
\item{object}{a \code{varstan} object.}

\item{h}{An integer indicating the number of predictions. The default number
of predictions is 12.}

\item{xreg}{Optionally, a numerical matrix of external regressors,
which must have the same number of rows as ts. It should not be a data frame.}

\item{robust}{a bool for obtain the robust estimation.}

\item{draws}{a integer indicating the number of draws to return. The default
number of draws is 1000.}

\item{seed}{An optional \code{seed} to use.}

\item{...}{Further arguments passed to  \code{posterior_predict}.}
}
\value{
An \code{array} of predicted \emph{mean} response values. For categorical and
ordinal models, the output is an S x N x C array. Otherwise, the output is an
\code{S x N} matrix, where S is the number of posterior samples, N is the number
of observations, and C is the number of categories. In multivariate models, an
additional  dimension is added to the output which indexes along the different
response variables.
}
\description{
Compute posterior samples of the expected value/mean of the posterior
predictive distribution. Can be performed for the data used to fit the model
(posterior predictive checks) or for new data. By definition, these
predictions have smaller variance than the posterior predictions performed by
the \code{posterior_predict.varstan} method. This is because only the
uncertainty in the mean is incorporated in the samples computed by
\code{posterior_epred} while any residual error is ignored. However, the
estimated means of both methods averaged across samples should be very
similar.
}
