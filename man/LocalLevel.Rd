% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ets.R
\name{LocalLevel}
\alias{LocalLevel}
\title{A  constructor for local level state-space model.}
\usage{
LocalLevel(ts, xreg = NULL, genT = FALSE, series.name = NULL)
}
\arguments{
\item{ts}{a numeric or ts object with the univariate time series.}

\item{xreg}{Optionally, a numerical matrix of external regressors,
which must have the same number of rows as ts. It should not be a data frame.}

\item{genT}{a boolean value to specify for a generalized t-student SSM model.}

\item{series.name}{an optional string vector with the time series names.}
}
\value{
The function returns a list with the data for running \code{stan()}
function of \pkg{rstan} package.
}
\description{
Constructor of the \code{ets("A","N","N")} object for Bayesian estimation in
\pkg{Stan}.
}
\details{
By  default  the  \code{ssm()} function generates a local-level, \code{ets("A","N","N")},
or exponential smoothing model from the \pkg{forecast} package. When
\code{trend = TRUE} the SSM transforms into a local-trend, \code{ets("A","A","N")},
or the equivalent Holt model. For damped trend models set \code{damped = TRUE}.
If \code{seasonal = TRUE}, the model is a seasonal local level model, or
\code{ets("A","N","A")} model. Finally, the Holt-Winters method (\code{ets("A","A","A")})
is obtained by setting both \code{Trend = TRUE} and \code{seasonal = TRUE}.

The \code{genT = TRUE} option generates a t-student innovations SSM model. For
a detailed explanation, check Ardia (2010); or Fonseca, et. al (2019).

The default priors used in a \code{ssm( )} model are:

\itemize{
\item{level ~ normal(0,0.5)}
\item{sigma0 ~ t-student(0,1,7)}
\item{level1 ~ normal(0,1)}
\item{dfv ~ gamma(2,0.1)}
\item{breg ~ t-student(0,2.5,6)}
}

For changing the default prior use the function \code{set_prior()}.
}
\examples{
mod1 = LocalLevel(ipc)

}
\references{
Fonseca, T. and Cequeira, V. and Migon, H. and Torres, C. (2019). The effects of
degrees of freedom estimation in the Asymmetric GARCH model with Student-t
Innovations. \emph{arXiv} \code{doi: arXiv: 1910.01398}.
}
\seealso{
\code{Sarima}, \code{auto.arima}, \code{set_prior}, and \code{garch}.
}
\author{
Asael Alonzo Matamoros.
}
