#' U.S. Monthly Live Births.
#'
#' Monthly live births (adjusted) in thousands for the United States, 1948-1979.
#'
#' @docType data
#'
#' @format the format is: Time-Series (1:373) from 1948 to 1979:
#'
#' @source \pkg{astsa}
#'
#' @references
#' http://www.stat.pitt.edu/stoffer/tsa4/
#' http://www.stat.pitt.edu/stoffer/tsda/
#'
"birth"

#' DEM/GBP exchange rate log-returns
#'
#' The vector dem2gbp contains daily observations of the Deutschmark vs British
#' Pound foreign exchange rate log-returns. This data set has been promoted as
#' an informal benchmark for GARCH time-series software validation. See
#' McCullough and Renfro (1999), and Brooks, Burke, and Persand (2001) for details.
#' The nominal returns are expressed in percent as in Bollerslev and Ghysels (1996).
#' The sample period is from January 3, 1984, to December 31, 1991, for a total
#' of 1974 observations.
#'
#' @docType data
#'
#' @format the format is: Time-Series (1:350) from 1984 to 1985:
#'
#' @source \pkg{bayesGARCH}
#'
#' @references
#' Engle, R. (1982). Autoregressive Conditional Heteroscedasticity with Estimates of
#' the Variance of United Kingdom Inflation. \emph{Econometrica}, 50(4), 987-1007.
#' \code{url: http://www.jstor.org/stable/1912773}.
#'
#' Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity.
#' \emph{Journal of Econometrics}. 31(3), 307-327.
#' \code{doi: https://doi.org/10.1016/0304-4076(86)90063-1}.
#'
#' Ardia, D. and Hoogerheide, L. (2010). Bayesian Estimation of the GARCH(1,1) Model
#' with Student-t Innovations. \emph{The R Journal}. 2(7), 41-47.
#' \code{doi: 10.32614/RJ-2010-014}.
"demgbp"

#' Monthly inflation coefficients from 1980-2018.
#'
#' Monthly return coefficients for the inflation. An economic indicator of
#' a country's economy.
#'
#' @docType data
#'
#' @format A time series of monthly data from 1980 to 2018.
#'
#' @source
#' https://www.bch.hn/series_estadisticas.php
#'
#'
#'
"ipc"

#' Annual oil production in Saudi Arabia
#'
#' Annual oil production (millions of tonnes), Saudi Arabia, 1965-2013.
#'
#' @docType data
#'
#' @format the format is: Annual Time-Series (1:18) from 1996 to 2013:
#'
#' @source \pkg{fpp2}
#'
#' @references
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}
"oildata"

#' Air Transport Passengers Australia
#'
#' Total annual air passengers (in millions) including domestic and international
#' aircraft passengers of air carriers registered in Australia. 1970-2016.
#'
#' @docType data
#'
#' @format the format is: Annual Time-Series (1:27) from 1990 to 2016:
#'
#' @source \pkg{fpp2}
#'
#' @references
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}
"air"

#' International Tourists to Australia: Total visitor nights.
#'
#' Quarterly visitor nights (in millions) spent by international tourists
#' to Australia. 1999-2015
#'
#' @docType data
#'
#' @format the format is: Quarterly Time-Series (1:44) from 1999 to 2015:
#'
#' @source \pkg{fpp2}
#'
#' @references
#' Hyndman, R. & Khandakar, Y. (2008). Automatic time series forecasting: the
#' forecast package for \code{R}. \emph{Journal of Statistical Software}. 26(3),
#' 1-22.\code{doi:	10.18637/jss.v027.i03}
"aust"
