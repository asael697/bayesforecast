% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SVM.R
\name{SVM}
\alias{SVM}
\title{Constructor of an Stochastic volatility model object}
\usage{
SVM(ts, arma = c(0, 0), xreg = NULL, series.name = NULL)
}
\arguments{
\item{ts}{a numeric or ts object with the univariate time series.}

\item{arma}{Optionally, a specification of the  ARMA model,same
as order parameter: the two components \code{c(p, q)} are the AR, and
the  MA orders.}

\item{xreg}{Optionally, a numerical matrix of external regressors,
which must have the same number of rows as ts. It should not be a data frame.}

\item{series.name}{an optional string vector with the time series names.}
}
\value{
The function returns a list with the data for running \code{stan()}
function of \pkg{rstan} package.
}
\description{
Constructor of the Stochastic Volatility model (SVM) for Bayesian estimation in \pkg{Stan}.
}
\details{
The function returns a list with the data for running \code{stan()} function of
\pkg{rstan} package.
}
\examples{
# Declares a SVM model for the IPC data

model = SVM(ipc, arma = c(1,1))
model

}
\references{
Sangjoon,K. and Shephard, N. and Chib.S (1998). Stochastic Volatility: Likelihood
Inference and Comparison with ARCH Models. \emph{Review of Economic Studies}.
65(1), 361-93. \code{url: https://www.jstor.org/stable/2566931}.

Tsay, R (2010). Analysis of Financial Time Series.
\emph{Wiley-Interscience}. 978-0470414354, second edition.

Shumway, R.H. and Stoffer, D.S. (2010).Time Series Analysis and Its
Applications: With R Examples. \emph{Springer Texts in Statistics}.
isbn: 9781441978646. First edition.
}
\seealso{
\code{garch}, and  \code{et_prior}.
}
\author{
Asael Alonzo Matamoros
}
