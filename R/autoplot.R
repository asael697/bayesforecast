#' Automatically create a ggplot for time series objects.
#'
#' \code{autoplot} takes an object of type \code{ts} or \code{mts} and creates
#' a ggplot object suitable for usage with \code{stat_forecast}.
#'
#' \code{fortify.ts} takes a \code{ts} object and converts it into a data frame
#' (for usage with ggplot2).
#'
#' @param object Object of class \dQuote{\code{ts}} or \dQuote{\code{mts}}.
#' @param series Identifies the time series with a colour, which integrates well
#' with the functionality of \link{geom_forecast}.
#' @param xlab a string with the plot's x axis label. By default a NUll value.
#' @param ylab a string with the plot's y axis label. By default a counts" value.
#' @param main a string with the plot's title.
#' @param facets If TRUE, multiple time series will be faceted (and unless
#' specified, colour is set to FALSE). If FALSE, each series will be assigned a
#' colour.
#' @param colour If TRUE, the time series will be assigned a colour aesthetic
#' @param model Object of class \dQuote{\code{ts}} to be converted to
#' \dQuote{\code{data.frame}}.
#' @param data Not used (required for \link{fortify} method)
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @author Mitchell O'Hara-Wild.
#'
#' @seealso \code{\link[stats]{plot.ts}}, \code{\link[ggplot2]{fortify}}
#'
#' @examples
#'
#' library(ggplot2)
#' autoplot(USAccDeaths)
#'
#' lungDeaths <- cbind(mdeaths, fdeaths)
#' autoplot(lungDeaths)
#' autoplot(lungDeaths, facets=TRUE)
#'
#' @importFrom forecast autoplot
#' @importFrom stats is.ts
#' @export
#'
autoplot.ts <- function(object, series=NULL, xlab = "Time", ylab = deparse(substitute(object)),
                        main = NULL,facets = FALSE,colour = TRUE, ...) {
  if (!stats::is.ts(object))
    stop("autoplot.ts requires a ts object, use object=object")

  # Create data frame with time as a column labelled x
  # and time series as a column labelled y.
  data <- data.frame(y = as.numeric(object), x = as.numeric(time(object)))
  if (!is.null(series)) {
    data <- transform(data, series = series)
  }

  # Initialise ggplot object
  p <- ggplot2::ggplot(ggplot2::aes_(y = ~y, x = ~x), data = data)

  # Add data
  if (!is.null(series)) {
    p <- p + ggplot2::geom_line(ggplot2::aes_(group = ~series, colour = ~series), na.rm = TRUE, ...)
  }
  else {
    p <- p + ggplot2::geom_line(na.rm = TRUE, ...)
  }

  # Add labels
  p <- p + ggplot2::labs(x = xlab, y = ylab, title = main)

  # Make x axis contain only whole numbers (e.g., years)
  p <- p + ggplot2::scale_x_continuous(breaks = ggtsbreaks)
  return(p)
}
#'
#' @rdname autoplot.ts
#' @importFrom zoo as.Date
#' @export
#'
fortify.ts = function(model, data, ...) {
  # Use ggfortify version if it is loaded
  # to prevent cran errors
  if (exists("ggfreqplot")) {
    tsp <- attr(model, which = "tsp")
    dtindex <- time(model)
    if (any(tsp[3] == c(4, 12))) {
      dtindex <- zoo::as.Date.yearmon(dtindex)
    }
    model <- data.frame(Index = dtindex, Data = as.numeric(model))
    return(ggplot2::fortify(model))
  }
  else {
    model <- cbind(x = as.numeric(time(model)), y = as.numeric(model))
    as.data.frame(model)
  }
}
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#'
#'
ggAddExtras = function(xlab=NA, ylab=NA, main=NA) {
  dots <- eval.parent(quote(list(...)))
  extras <- list()
  if ("xlab" %in% names(dots) || is.null(xlab) || any(!is.na(xlab))) {
    if ("xlab" %in% names(dots)) {
      extras[[length(extras) + 1]] <- ggplot2::xlab(dots$xlab)
    }
    else {
      extras[[length(extras) + 1]] <- ggplot2::xlab(paste0(xlab[!is.na(xlab)], collapse = " "))
    }
  }
  if ("ylab" %in% names(dots) || is.null(ylab) || any(!is.na(ylab))) {
    if ("ylab" %in% names(dots)) {
      extras[[length(extras) + 1]] <- ggplot2::ylab(dots$ylab)
    }
    else {
      extras[[length(extras) + 1]] <- ggplot2::ylab(paste0(ylab[!is.na(ylab)], collapse = " "))
    }
  }
  if ("main" %in% names(dots) || is.null(main) || any(!is.na(main))) {
    if ("main" %in% names(dots)) {
      extras[[length(extras) + 1]] <- ggplot2::ggtitle(dots$main)
    }
    else {
      extras[[length(extras) + 1]] <- ggplot2::ggtitle(paste0(main[!is.na(main)], collapse = " "))
    }
  }
  if ("xlim" %in% names(dots)) {
    extras[[length(extras) + 1]] <- ggplot2::xlim(dots$xlim)
  }
  if ("ylim" %in% names(dots)) {
    extras[[length(extras) + 1]] <- ggplot2::ylim(dots$ylim)
  }
  return(extras)
}
#'
#' @noRd
#'
ggtsbreaks = function(x) {
  # Make x axis contain only whole numbers (e.g., years)
  return(unique(round(pretty(floor(x[1]):ceiling(x[2])))))
}
#' Histogram with optional normal density functions
#'
#' Plots a histogram and density estimates using ggplot.
#'
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param title a string with the plot's title.
#' @param xlab a string with the plot's x axis label. By default a NUll value
#' @param ylab a string with the plot's y axis label. By default a "counts" value
#' @param add.normal A boolean value. Add a normal density function for comparison,
#' by default \code{add.normal = TRUE}.
#' @param bins The number of bins to use for the histogram. Selected by default
#' using the Friedman-Diaconis rule.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @author Rob J Hyndman
#'
#' @importFrom grDevices nclass.FD
#' @importFrom stats dnorm is.ts na.exclude
#' @export
#'
#' @examples
#' x = rnorm(100)
#' gghist(x,add.normal = TRUE)
#'
gghist = function(y,title = NULL,xlab = NULL,ylab = "counts",bins,add.normal = TRUE){

  if (!stats::is.ts(y) && !is.numeric(y))
    stop("gghist requires a ts or numeric object")

  if (missing(bins))
    bins = min(grDevices::nclass.FD(stats::na.exclude(y)),500)

  xlab1 = xlab
  if(is.null(xlab))
    xlab1 =  deparse(substitute(y))

  data = data.frame(y = as.numeric(c(y)))
  # Initialise ggplot object and plot histogram

  boundary=0
  binwidth = (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)) / bins

  p = ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(y),data = data,
                            binwidth = binwidth, boundary = boundary)

  if (add.normal) {
    xmin  = min(y, na.rm = TRUE)
    xmax  = max(y, na.rm = TRUE)
    xmean = mean(y, na.rm = TRUE)
    xsd   = sd(y, na.rm = TRUE)
    xmin  = min(xmin, xmean - 3 * xsd)
    xmax  = max(xmax, xmean + 3 * xsd)
  }
  xgrid  =  seq(xmin, xmax, length.out = 512)

  if (add.normal) {
    df = data.frame(x = xgrid, y = length(y)*binwidth*stats::dnorm(xgrid, xmean, xsd))
    p = p + ggplot2::geom_line(ggplot2::aes(df$x, df$y), col = "blue")
  }

  p = p + ggplot2::labs(title = title,x = xlab1 ,y = ylab)

  return(p)

}
#' \code{qqplot} with normal \code{qqline}
#'
#' Plot the quantile-quantile plot and quantile-quantile line using ggplot.
#'
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param add.normal Add a normal density function for comparison.
#' @param title a string with the plot's title.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @author Asael Alonzo Matamoros
#'
#' @import ggplot2
#' @importFrom stats is.ts
#' @export
#'
#' @examples
#' x = rnorm(100)
#' ggnorm(x)
#'
ggnorm = function(y,title = NULL,add.normal = TRUE){

  if (!stats::is.ts(y) && !is.numeric(y))
    stop("gghist requires a ts or numeric object")

  df = data.frame(y = as.numeric(y))

  p = ggplot2::ggplot(data = df,ggplot2::aes(sample = y)) + ggplot2::stat_qq() +
    ggplot2::labs(title = title)

  if(add.normal)
    p = p + ggplot2::stat_qq_line(color = "blue")

  return(p)
}
#' \code{acf} plot
#'
#' Plot of the auto-correlation function for a univariate time series.
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param title a string with the plot's title.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @author Asael Alonzo Matamoros
#'
#' @import forecast
#' @importFrom stats is.ts
#' @export
#'
#' @examples
#' x = rnorm(100)
#' ggacf(x)
#'
ggacf = function(y,title = NULL){
  if (!stats::is.ts(y) && !is.numeric(y))
    stop("gghist requires a ts or numeric object")

  p = forecast::ggAcf(x = y) + ggplot2::labs(title = title)

  return(p)
}
#' \code{pacf} plot.
#'
#' Plot of the partial autocorrelation function for a univariate time series.
#'
#' @param y a numeric vector or an object of the \code{ts} class containing a stationary time series.
#' @param title a string with the plot's title.
#'
#' @return None.
#'
#' @author Mitchell O'Hara-Wild and Asael Alonzo Matamoros
#'
#' @import forecast
#' @importFrom stats is.ts
#' @export
#'
#' @examples
#' x = rnorm(100)
#' ggpacf(x)
#'
ggpacf = function(y,title = NULL){
  if (!stats::is.ts(y) && !is.numeric(y))
    stop("gghist requires a ts or numeric object")

  p = forecast::ggPacf(x = y) + ggplot2::labs(title = title)

  return(p)
}
#' plot methods for varstan models.
#'
#' Preliminary plot methods for varstan models only valid for univariate time series models.
#' The function prints the fitted values time series, the trace and density plots for the
#' sampled model parameters, or the residuals' posterior mean time series.
#'
#' @param x An object of class \code{varstan}.
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \eqn{90\%} intervals (\code{prob=0.9}) rather than the traditionally used
#'   \eqn{95\%}.
#' @param ... Further arguments passed to  \code{mcmc_combo}.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  sf1 = auto.sarima(ts = birth)
#'  # fitted model
#'  plot(sf1)
#' }
#'
#' @method plot varstan
#' @import ggplot2
#' @importFrom bayesplot mcmc_combo
#' @importFrom stats quantile
#' @export
#'
plot.varstan = function(x,prob = 0.95,...){

  if( !is.varstan(x))
    stop("The current object is not a varstan class")

  g = autoplot.varstan(object = x, prob = prob,...)

  return(g)
}
#' autoplot methods for varstan models.
#'
#' Preliminary autoplot methods for varstan models only valid for univariate time series models.
#' The function prints the fitted values time series, the trace and density plots for the
#' sampled model parameters, or the residuals' posterior mean time series.
#'
#' @param object An object of class \code{varstan}.
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \code{90\%} intervals (\code{prob=0.9}) rather than the traditionally used
#'   \code{95\%}.
#' @param ... Further arguments passed to  \code{mcmc_combo}.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  sf1 = auto.sarima(ts = birth)
#'  # fitted model
#'  autoplot(sf1)
#' }
#'
#' @method autoplot varstan
#' @import ggplot2
#' @importFrom bayesplot mcmc_combo
#' @importFrom stats quantile
#' @export
#'
autoplot.varstan = function(object,prob = 0.95,...){
  if( !is.varstan(object))
    stop("The current object is not a varstan class")

  pp = posterior_predict(object)
  pI = posterior_interval(pp,prob = prob)
  pM = apply(pp, 2, FUN = mean)

  data = data.frame(time = as.numeric(time(object$ts)),y = as.numeric(object$ts))
  data$yhat = pM; data = cbind(data,pI)
  colnames(data) = c("time","y","yhat","low","hi")

  colors = c("yhat" = "#0000CC", "y" = "#000000")

  g = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~yhat),data = data)+
    ggplot2::geom_line(ggplot2::aes_(y = ~yhat,color = "yhat"))+
    ggplot2::geom_point(ggplot2::aes_(y = ~y,color = "y"))+
    ggplot2::geom_smooth(ggplot2::aes_(ymin = ~low, ymax = ~hi),fill="#333333", color="#0000CC",
                         stat = "identity",size = 0.5)+
    ggplot2::labs(x = "time",y = object$series.name,color = "Legend",title = "posterior predict") +
    ggplot2::scale_color_manual(values = colors)

  return(g)
}
#' Visual check of residuals in a \code{varstan} object.
#'
#'
#' Performs a visual check of residuals in time series models, this method is inspired in
#' the \code{check.residuals} function provided by the \code{forecast} package.
#'
#' @param object a varstan object.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return None. Function produces a ggplot2 graph.
#'
#' @author Asael Alonzo Matamoros.
#'
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' \donttest{
#'  library(astsa)
#'  sf1 = auto.sarima(ts = birth)
#'  # fitted model
#'  check_residuals(sf1)
#' }
#'
check_residuals = function(object,...){

  if( !is.varstan(object))
    stop("The current object is not a varstan class")
  res = residuals.varstan(object = object,...)

  lay = matrix(c(1,1,2,3,4,5),nrow = 3,ncol = 2,byrow = TRUE)

  p1 = autoplot(res,main = "Expected Values of the Posterior Predictive Errors",ylab = "residuals")
  p2 = gghist(y = res,add.normal = TRUE,xlab = "residuals")
  p3 = ggnorm(y = res,add.normal = TRUE,)
  p4 = ggacf(y = res)
  p5 = ggpacf(y = res)

  grob = list(p1,p2,p3,p4,p5)
  gridExtra::grid.arrange(grobs = grob,ncol=2,nrow = 3,layout_matrix = lay)

}
#' MCMC Plots Implemented in \pkg{bayesplot}
#'
#' Convenient way to call MCMC plotting functions
#' implemented in the \pkg{bayesplot} package.
#'
#' @param object An \code{varstan} object.
#' @param pars Names of parameters to be plotted,
#'   as given by a character vector or regular expressions.
#'   By default, all parameters except for group-level and
#'   smooth effects are plotted. May be ignored for some plots.
#' @param combo An array that contains the types of plot. By default
#'   combo = c("dens","trace"). Supported types are (as names) \code{hist},
#'   \code{dens}, \code{hist_by_chain}, \code{dens_overlay},
#'   \code{violin}, \code{intervals}, \code{areas}, \code{acf},
#'   \code{acf_bar},\code{trace}, \code{trace_highlight}, \code{scatter},
#'   \code{rhat}, \code{rhat_hist}, \code{neff}, \code{neff_hist}
#'   \code{nuts_acceptance}, \code{nuts_divergence},
#'   \code{nuts_stepsize}, \code{nuts_treedepth}, and \code{nuts_energy}.
#'   For an overview on the various plot types see
#'   \code{\link[bayesplot:MCMC-overview]{MCMC-overview}}.
#' @param ... Additional arguments passed to the plotting functions.
#'   See \code{\link[bayesplot:MCMC-overview]{MCMC-overview}} for
#'   more details.
#'
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object
#'   that can be further customized using the \pkg{ggplot2} package.
#'
#' @export
#' @import bayesplot
#' @examples
#' \dontrun{
#' sf1 = stan_ssm(ipc)
#'
#' # plot posterior intervals
#' mcmc_plot(sf1)
#'
#' # only show population-level effects in the plots
#' mcmc_plot(sf1, pars = "level")
#' }
#'
mcmc_plot.varstan = function(object, pars = NULL, combo = c("dens","trace"),
                              fixed = FALSE, exact_match = FALSE, ...) {

  if( !is.varstan(object))
    stop("The current object is not a varstan class")

  gp = get_parameters(object)
  code = match(x = pars,table = gp)

  if(is.element(NA,code))
    stop("pars contains an incorrect value")

  if(is.null(pars)) pars = gp

  x = as.data.frame(extract_stan(object = object,pars = pars))

  g = bayesplot::mcmc_combo(x = x,pars = pars,combo = combo)

  return(g)
}

#' @rdname mcmc_plot.varstan
#' @export
mcmc_plot = function(object, ...) {
  UseMethod("mcmc_plot")
}
