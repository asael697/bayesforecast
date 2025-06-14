% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sarima.R
\name{Sarima}
\alias{Sarima}
\title{Constructor a Multiplicative Seasonal ARIMA model.}
\usage{
Sarima(
  ts,
  order = c(1, 0, 0),
  seasonal = c(0, 0, 0),
  xreg = NULL,
  period = 0,
  series.name = NULL
)
}
\arguments{
\item{ts}{a numeric or ts object with the univariate time series.}

\item{order}{a three length vector with the specification of the non-seasonal
part of the ARIMA model: the three components \code{c(p, d, q)} are the AR,number of
differences, and the MA orders respectively.}

\item{seasonal}{a vector of length three with the  specification of the seasonal
part of the SARIMA model. The three components \code{c(P, D, Q)} are the seasonal AR,
the degree of seasonal differences, and the seasonal MA orders respectively.}

\item{xreg}{Optionally, a numerical matrix of external regressors,
which must have the same number of rows as ts. It should not be a data frame.}

\item{period}{an integer specifying the periodicity of the time series by
default the value frequency(ts) is used.}

\item{series.name}{an optional string vector with the series names.}
}
\value{
The function returns a list with the data for running \code{stan()}
function of \pkg{rstan} package.
}
\description{
Constructor of the SARIMA model for Bayesian estimation in \pkg{Stan}.
}
\details{
The function returns a list with the data for running \code{stan} function of
\pkg{rstan} package

If \code{xreg} option is used, the model by default will cancel the
seasonal differences adjusted (D = 0). If a value \code{d} > 0 is used, all
the regressor variables in \code{xreg} will be difference as well.

The default priors used in Sarima are:

\itemize{
\item{ar ~ normal(0,0.5)}
\item{ma ~ normal(0,0.5)}
\item{mu0 ~ t-student(0,2.5,6)}
\item{sigma0 ~ t-student(0,1,7)}
\item{sar ~ normal(0,0.5)}
\item{sma ~ normal(0,0.5)}
\item{breg ~ t-student(0,2.5,6)}
}

For changing the default prior use the function \code{set_prior}.
}
\examples{
# Declare a multiplicative seasonal ARIMA model for the birth data.

model = Sarima(birth,order = c(0,1,2),seasonal = c(1,1,1))
model

#Declare an Dynamic Harmonic Regression model for the birth data.
model = Sarima(birth,order = c(1,0,1),xreg = fourier(birth,K = 2))
model

}
\references{
Box, G. E. P. and Jenkins, G.M. (1978). Time series analysis: Forecasting and
control. San Francisco: Holden-Day. \emph{Biometrika}, 60(2), 297-303.
\code{doi:10.1093/biomet/65.2.297}.

Kennedy, P. (1992). Forecasting with dynamic regression models: Alan Pankratz, 1991.
\emph{International Journal of Forecasting}. 8(4), 647-648.
\code{url: https://EconPapers.repec.org/RePEc:eee:intfor:v:8:y:1992:i:4:p:647-648}.

Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
1-22.\code{doi:	10.18637/jss.v027.i03}
}
\seealso{
\code{garch}, and \code{set_prior} functions.
}
\author{
Asael Alonzo Matamoros
}
