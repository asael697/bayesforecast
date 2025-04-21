#' Draw from posterior predictive h steps ahead distribution
#'
#' The posterior predictive distribution is the distribution of the outcome
#' implied by the model after using the observed data to update our beliefs
#' about the unknown parameters in the model. Simulating data from the posterior
#' predictive distribution using the observed predictors is useful for checking
#' the fit of the model. Drawing from the posterior predictive distribution at
#' interesting values of the predictors also lets us visualize how a
#' manipulation of a predictor affects (a function of) the outcome(s). With new
#' observations of predictor variables we can use the posterior predictive
#' distribution to generate predicted outcomes.
#'
#' @aliases posterior_predict
#'
#' @param object a `varstan` object
#' @param h an integer indicating the number of predictions. The default number
#' of predictions is 12.
#' @param xreg	Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param robust a bool for obtain the robust estimation.
#' @param draws an integer indicating the number of draws to return. The default
#' number of draws is 1000.
#' @param seed an optional \code{seed} to use.
#' @param ... Further arguments passed to  \code{posterior_predict}.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return
#' A \code{draws} by \code{h} data.frame of simulations from the
#' posterior predictive distribution. Each row of the data.frame
#' is a vector of predictions generated using a single draw of
#' the model parameters from the posterior distribution.
#'
#' @importFrom stats arima predict
#' @importFrom MASS mvrnorm
#'
#' @importFrom rstantools posterior_predict
#' @method posterior_predict varstan
#' @export
#' @export posterior_predict
#'
posterior_predict.varstan = function(object, h = 0, xreg = NULL, robust = FALSE,
                                     draws = 1000, seed = NULL, ...){

  if (! is.varstan(object))
    stop("The current object is not a varstan class", call. = FALSE)

  if( h == 0){
    fc = as.data.frame(extract_stan(object,"fit", permuted = TRUE) )
  }
  else{
    if(is.garch(object$model))
      fc = posterior_predict_garch(object = object, h = h, xreg = xreg,
                                   robust = robust, draws = draws, seed = seed)

    if(is.SVM(object$model) )
      fc = posterior_predict_SVM(object = object, h = h, xreg = xreg,
                                 robust = robust, draws = draws, seed = seed)

    if(is.Sarima(object$model))
      fc = posterior_predict_Sarima(object = object, h = h, xreg = xreg,
                                    robust = robust, draws = draws, seed = seed)

    if(is.naive(object$model))
      fc = posterior_predict_Sarima(object = object, h = h, xreg = xreg,
                                    robust = robust, draws = draws, seed = seed)

    if(is.ssm(object$model))
      fc = posterior_predict_ets(object = object, h = h, xreg = xreg,
                                 robust = robust, draws = draws, seed = seed)

    if(is.LocalLevel(object$model))
      fc = posterior_predict_ets(object = object, h = h, xreg = xreg,
                                 robust = robust, draws = draws, seed = seed)

    if(is.Holt(object$model))
      fc = posterior_predict_ets(object = object, h = h, xreg = xreg,
                                 robust = robust, draws = draws, seed = seed)

    if(is.Hw(object$model))
      fc = posterior_predict_ets(object = object, h = h, xreg = xreg,
                                 robust = robust, draws = draws, seed = seed)
  }
  return(fc)
}
#' Expected Values of the Posterior Predictive Distribution
#'
#' The function returns the posterior estimate of the fitted values
#' of a \code{varstan} model, similar to the fit_values functions of other
#' packages.
#'
#' @param object a `varstan` object.
#' @param robust a bool value, if its \code{TRUE} it returns the median of the
#' posterior distribution, and if its \code{FALSE} it returns the mean, by default
#' is the \code{FALSE} value.
#' @param ... Further arguments passed to  \code{posterior_predict}.
#'
#' @details
#' This function returns a time series of the predicted \emph{mean} response values.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return
#' A time series \code{(ts)} of predicted \emph{mean} response values.
#'
#' @importFrom stats median ts
#'
#' @seealso \code{posterior_predict.varstan}
#'
#' @export
#'
fitted.varstan = function(object, robust = FALSE, ...){
  if( !is.varstan(object) )
    stop("The current object is not a varstan class")

  post = as.data.frame(extract_stan(object,"fit", permuted = TRUE) )
  if(robust)
    sum1 = t(matrix(apply(post,2,stats::median), nrow = 1,byrow = TRUE))
  else
    sum1 = t(matrix(apply(post,2,mean),nrow = 1, byrow = TRUE))

  fit = stats::ts(sum1,start = stats::start(object$ts),
                  frequency = stats::frequency(object$ts))

  return(fit)
}
#' Expected Values of the Posterior Predictive Distribution
#'
#' Compute posterior samples of the expected value/mean of the posterior
#' predictive distribution. Can be performed for the data used to fit the model
#' (posterior predictive checks) or for new data. By definition, these
#' predictions have smaller variance than the posterior predictions performed by
#' the \code{\link{posterior_predict.varstan}} method. This is because only the
#' uncertainty in the mean is incorporated in the samples computed by
#' \code{posterior_epred} while any residual error is ignored. However, the
#' estimated means of both methods averaged across samples should be very
#' similar.
#'
#' @aliases posterior_epred
#'
#' @param object a `varstan` object.
#' @param h An integer indicating the number of predictions. The default number
#' of predictions is 12.
#' @param xreg	Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as ts. It should not be a data frame.
#' @param robust a bool for obtain the robust estimation.
#' @param draws a integer indicating the number of draws to return. The default
#' number of draws is 1000.
#' @param seed An optional \code{seed} to use.
#' @param ... Further arguments passed to  \code{posterior_predict}.
#'
#' @return
#' An \code{array} of predicted \emph{mean} response values. For categorical and
#' ordinal models, the output is an S x N x C array. Otherwise, the output is an
#' \code{S x N} matrix, where S is the number of posterior samples, N is the number
#' of observations, and C is the number of categories. In multivariate models, an
#' additional  dimension is added to the output which indexes along the different
#' response variables.
#'
#' @method posterior_epred varstan
#' @importFrom rstantools posterior_epred
#' @export posterior_epred
#' @export
#'
posterior_epred.varstan = function(object, h = 0, xreg = NULL, robust = FALSE,
                                   draws = 1000, seed = NULL, ...){

  if (! is.varstan(object))
    stop("The current object is not a varstan class",call. = FALSE)

  pp = posterior_predict(object = object, h = h, xreg = xreg, robust = robust,
                         draws = draws, seed = seed,...)
  if(robust)
    ppe = t(matrix(apply(pp,2,stats::median),nrow = 1,byrow = TRUE))
  else
    ppe = t(matrix(apply(pp,2,mean),nrow = 1,byrow = TRUE))

  return(ppe)
}

#############################################################################
# Internal
#############################################################################

#' Draw from posterior predictive distribution of an Seasonal ARIMA model.
#'
#' @param object a `varstan` object
#' @param h An integer indicating the number of predictions. The default number
#' of predictions is 1.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as h. It should not be a data frame.
#' @param robust a bool for obtain the robust estimation.
#' @param draws an integer indicating the number of draws to return. The default
#'    number of draws is 1000.
#' @param seed An optional \code{seed} to use.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return
#' A \code{draws} by \code{h} data.frame of simulations from the
#' posterior predictive distribution. Each row of the data.frame is a vector of
#' predictions generated using a single draw of the model parameters from the
#' posterior distribution.
#'
#' @importFrom stats arima predict rnorm
#' @noRd
#'
posterior_predict_Sarima = function(object, h = 1, xreg = NULL, robust = TRUE,
                                    draws = 1000, seed = NULL){

  if (!is.null(seed))
    set.seed(seed)

  nm = object$stan.parmaters$chains * (object$stan.parmaters$iter - object$stan.parmaters$warmup)
  draw = draws
  if( nm < draws) draw = nm

  sp = sample(1:nm,size = draw)

  # Extract the necessary lags for predict
  order = get_order(object);
  fix = NULL
  if(order$p > 0) fix =c(fix,"ar")
  if(order$q > 0) fix =c(fix,"ma")
  if(order$P > 0) fix =c(fix,"sar")
  if(order$Q > 0) fix =c(fix,"sma")
  if(order$d1 > 0)fix =c(fix,"breg")
  y = object$ts;

  # point estimate of the model parameters
  par0 = data.frame(extract_stan(object, pars = c("mu0","sigma0")))
  par0 = as.matrix(par0[sp,])

  if(!is.null(fix)){
    fix = data.frame(extract_stan(object, pars = fix))
    fix = data.frame(fix[sp,])
  }
  else{
    fix = matrix(data = 0,nrow = nm,ncol = 2)
    fix = data.frame(fix[sp,])
    order$p = 2
  }

  # The previous data
  yh = matrix(0,nrow = draw,ncol = h)
  reg = NULL
  xh = NULL

  #preliminary checks
  if( order$d1 > 0){
    reg = object$model$reg
    # Check xreg
    if(is.null(xreg)){
      warning("No xreg specified, the forecast wont be accurate \n")
      xh  = matrix(0, nrow = h, ncol = order$d1)
    }
    else if( dim(xreg)[1] != h ||  dim(xreg)[2] != order$d1){
      # Check xreg dimensions
      warning("The dimension of xreg are not correct, the forecast wont be accurate \n")
      xh = matrix(0, nrow = h ,ncol = order$d1)
    }
    else xh = xreg
  }
  for(i in 1:draw){
    modi = suppressWarnings(stats::arima(x = y,order = c(order$p, order$d, order$q),
                                         seasonal =list(order = c(order$P,
                                                                  order$D,
                                                                  order$Q),
                                                        period = stats::frequency(y)),
                                         xreg = reg,include.mean = FALSE,
                                         fixed = fix[i,]))
    yh[i,] = suppressWarnings(as.numeric(stats::predict(modi, n.ahead = h, newxreg = xh)$pred))
    for (j in 1:h) yh[i,j] = stats::rnorm(n = 1,mean = yh[i,j] + par0[i,1], sd = par0[i,2])
  }
  colnames(yh) = paste0("yh.",1:h)
  yh = as.data.frame(yh)

  return(yh);
}
#' Draw from posterior predictive distribution of an ARMA-GARCH model
#'
#' @param object a `varstan` object
#' @param h An integer indicating the number of predictions. The default number
#' of predictions is 1.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as h. It should not be a data frame.
#' @param robust a bool for obtain the robust estimation.
#' @param draws an integer indicating the number of draws to return. The default
#'    number of draws is 1000
#' @param seed An optional \code{seed} to use.
#'
#' @author Asael Alonzo Matamoros
#'
#'
#' @return
#' A \code{draws} by \code{h} data.frame of simulations from the
#' posterior predictive distribution. Each row of the data.frame is a vector of
#' predictions generated using a single draw of the model parameters from the
#' posterior distribution.
#'
#' @importFrom stats arima predict rnorm
#' @noRd
#'
posterior_predict_garch = function(object, h = 1, xreg = NULL, robust = TRUE,
                                   draws = 1000, seed = NULL){

  if (!is.null(seed))
    set.seed(seed)

  nm = object$stan.parmaters$chains * (object$stan.parmaters$iter - object$stan.parmaters$warmup)
  draw = draws
  if( nm < draws) draw = nm

  sp = sample(1:nm,size = draw)

  # Extract the necessary lags for predict
  order = get_order(object);
  n = object$model$n

  # Extract the posterior values
  sigma = extract_ts(object,"sigma", lag = object$model$n, robust = robust)
  sigma = as.numeric(sigma^2)
  y = object$ts;
  reg = NULL
  xh = NULL

  fixmean = NULL
  if(order$p > 0) fixmean =c(fixmean,"ar")
  if(order$q > 0) fixmean =c(fixmean,"ma")
  if(order$d1 > 0)fixmean =c(fixmean,"breg")


  fixsig = NULL
  if(order$s > 0) fixsig =c(fixsig,"arch")
  if(order$k > 0) fixsig =c(fixsig,"garch")

  fixmgarch = NULL

  # point estimate of the model parameters
  par0 = data.frame(extract_stan(object,pars = c("mu0","sigma0")))
  par0 = data.frame(par0[sp,])

  if(!is.null(fixsig)){
    fixsig = data.frame(extract_stan(object, pars = fixsig))
    fixsig = data.frame(fixsig[sp,])
  }

  if(!is.null(fixmean)){
    fixmean = data.frame(extract_stan(object, pars = fixmean))
    fixmean = data.frame(fixmean[sp,])
  }

  if(order$h > 0){
    fixmgarch = data.frame(extract_stan(object,pars = "mgarch"))
    fixmgarch = data.frame(fixmgarch[sp,])
  }

  #preliminary checks
  if( order$d1 > 0){
    reg = object$model$xreg
    # Check xreg
    if(is.null(xreg)){
      warning("No xreg specified, the forecast wont be accurate \n")
      xh  = matrix(0,nrow = h,ncol = order$d1)
    }
    else if( dim(xreg)[1] != h ||  dim(xreg)[2] != order$d1){
      # Check xreg dimensions
      warning("The dimension of xreg are not correct, the forecast wont be accurate \n")
      xh = matrix(0,nrow = h,ncol = order$d1)
    }
    else xh = xreg
  }

  # The previous data
  yh  =  matrix(0,nrow = draw,ncol = h)
  muh =  matrix(0,nrow = draw,ncol = h)
  sigh = matrix(0,nrow = draw,ncol = h)

  for(i in 1:draw){

    if(!is.null(fixmean)){
      modi = stats::arima(x = y,
                          order = c(order$p,0,order$q),
                          xreg = reg,include.mean = FALSE,
                          fixed = fixmean[i,])

      muh[i,] = as.numeric(stats::predict(modi,n.ahead = h,newxreg = xh)$pred)
    }
    modi = stats::arima(x = sigma,
                        order = c(order$s,0,order$k),
                        include.mean = FALSE,
                        fixed = fixsig[i,])

    sigh[i,] = as.numeric(stats::predict(modi,n.ahead = h)$pred)
    sigh[i,] = sqrt(abs(sigh[i,])+par0[i,2])

    if(order$h > 0){
      for(j in 1:h){
        for(k in 1:order$h){
          if(j-k+1 > 0)
            muh[i,j] = muh[i,j] + fixmgarch[i,k]*sigh[i,j-k+1]
          else
            muh[i,j] = muh[i,j] + fixmgarch[i,k]*sqrt(sigma[n-(j-k+1)])
        }
      }
    }
    for (j in 1:h) yh[i,j] = stats::rnorm(n = 1,mean = muh[i,j]+par0[i,1],sd = sigh[i,j])
  }
  colnames(yh) = paste0("yh.",1:h)
  yh = as.data.frame(yh)
  return(yh)
}
#' Draw from posterior predictive distribution of an ARMA-SVM model
#'
#' @param object a `varstan` object
#' @param h an integer indicating the number of predictions. The default number
#' of predictions is 1.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as h. It should not be a data frame.
#' @param robust a bool for obtain the robust estimation.
#' @param draws an integer indicating the number of draws to return. The default
#' number of draws is 1000.
#' @param seed An optional \code{seed} to use.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return
#' A \code{draws} by \code{h} data.frame of simulations from the
#' posterior predictive distribution. Each row of the data.frame is a vector of
#' predictions generated using a single draw of the model parameters from the
#' posterior distribution.
#'
#' @importFrom stats arima predict rnorm
#' @noRd
#'
posterior_predict_SVM = function(object, h = 1, xreg = NULL, robust = TRUE,
                                 draws = 1000, seed = NULL){

  if (!is.null(seed))
    set.seed(seed)

  nm = object$stan.parmaters$chains * (object$stan.parmaters$iter - object$stan.parmaters$warmup)
  draw = draws
  if( nm < draws) draw = nm

  sp = sample(1:nm,size = draw)

  # Extract the necessary lags for predict
  order = get_order(object);n = object$model$n

  # Extract the posterior values
  ht = extract_ts(object,"h",lag = object$model$n,robust = robust)

  y = object$ts;
  reg = NULL
  xh = NULL

  fixmean = NULL
  if(order$p > 0) fixmean =c(fixmean,"ar")
  if(order$q > 0) fixmean =c(fixmean,"ma")
  if(order$d1 > 0)fixmean =c(fixmean,"breg")


  fixsig =c("beta","alpha")

  fixmgarch = NULL

  # point etimate of the model parameters
  par0 = data.frame(extract_stan(object,pars = c("mu0","sigma0")))
  par0 = par0[sp,]

  fixsig = data.frame(extract_stan(object,pars = fixsig))
  fixsig = data.frame(fixsig[sp,])

  if(!is.null(fixmean)){
    fixmean = data.frame(extract_stan(object,pars = fixmean))
    fixmean = data.frame(fixmean[sp,])
  }

  #preliminary checks
  if( order$d1 > 0){
    reg = object$model$xreg
    # Check xreg
    if(is.null(xreg)){
      warning("No xreg specified, the forecast wont be accurate \n")
      xh  = matrix(0,nrow = h,ncol = order$d1)
    }
    else if( dim(xreg)[1] != h ||  dim(xreg)[2] != order$d1){
      # Check xreg dimensions
      warning("The dimension of xreg are not correct, the forecast wont be accurate \n")
      xh = matrix(0,nrow = h,ncol = order$d1)
    }
    else xh = xreg
  }

  # The previous data
  yh  =  matrix(0,nrow = draw,ncol = h)
  muh =  matrix(0,nrow = draw,ncol = h)
  sigh = matrix(0,nrow = draw,ncol = h)

  for(i in 1:draw){
    if(!is.null(fixmean)){
      modi = stats::arima(x = y,
                          order = c(order$p,0,order$q),
                          xreg = reg,include.mean = FALSE,
                          fixed = fixmean[i,])

      muh[i,] = as.numeric(stats::predict(modi,n.ahead = h,newxreg = xh)$pred)
    }

    modi = stats::arima(x = ht,
                        order = c(1,0,0),
                        include.mean = TRUE,
                        fixed = fixsig[i,])

    sigh[i,] = as.numeric(stats::predict(modi,n.ahead = h)$pred)
    for (j in 1:h) sigh[i,j] = stats::rnorm(n=1,mean = sigh[i,j],sd = par0[i,2])
    sigh[i,] = exp(sigh[i,]/2)

    for (j in 1:h) yh[i,j] = stats::rnorm(n = 1,mean = muh[i,j]+par0[i,1],sd = sigh[i,j])
  }
  colnames(yh) = paste0("yh.",1:h)
  yh = as.data.frame(yh)
  return(yh)
}
#' Draw from posterior predictive distribution of a SSM model
#'
#' @param object a `varstan` object.
#' @param h an integer indicating the number of predictions. The default number
#' of predictions is 1.
#' @param xreg Optionally, a numerical matrix of external regressors,
#' which must have the same number of rows as h. It should not be a data frame.
#' @param robust a bool for obtain the robust estimation. The default
#' @param draws an integer indicating the number of draws to return. The default
#' number of draws is 1000.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return
#' A \code{draws} by \code{h} data.frame of simulations from the
#' posterior predictive distribution. Each row of the data.frame is a vector of
#' predictions generated using a single draw of the model parameters from the
#' posterior distribution.
#'
#' @importFrom stats rnorm
#' @noRd
#'
posterior_predict_ets = function(object, h = 1, xreg = NULL, robust = TRUE,
                                 draws = 1000, seed = NULL){
  if (!is.null(seed))
    set.seed(seed)

  # correct number of draws
  nm = object$stan.parmaters$chains * (object$stan.parmaters$iter - object$stan.parmaters$warmup)
  draw = draws
  if(nm < draws) draw = nm
  #preliminary checks
  n = length(object$ts);d1 = object$model$d1;states = pars = NULL
  m = object$model$period; n1 = n-m+1;
  trendtype = object$model$is_td;damped = object$model$is_dp
  seasontype = object$model$is_ss

  if(d1 > 0){
    # Check xreg
    if(is.null(xreg)){
      warning("No xreg specified, the forecast wont be accurate \n")
      xh  = matrix(0,nrow = h,ncol = d1)
    }
    else if( dim(xreg)[1] != h ||  dim(xreg)[2] != d1){
      # Check xreg dimensions
      warning("The dimension of xreg are not correct, the forecast wont be accurate \n")
      xh = matrix(0,nrow = h,ncol = d1)
    }
    else xh = xreg

    breg = data.frame(extract_stan(object = object,pars = "breg"))
  }

  # Extract level
  l = data.frame(extract_stan(object = object,pars = "l"))
  states = cbind(states,l[,n])
  colnames(states) = "l"

  # Extract level parameters
  pars = data.frame(extract_stan(object = object,pars = c("sigma0","level")))
  colnames(pars) = c("sigma","alpha")

  # Extract trend
  if(trendtype){
    # trend states
    l = data.frame(extract_stan(object = object,pars = "b"))
    states = cbind(states,l[,n])
    colnames(states) = c("l","b")

    #trend parameters
    pars = cbind(pars,data.frame(extract_stan(object = object,pars = "trend")))
    colnames(pars) = c("sigma","alpha","beta")

    # damped
    if(damped){
      pars = cbind(pars,data.frame(extract_stan(object = object,pars = "damped")))
      colnames(pars) = c("sigma","alpha","beta","phi")
    }
  }

  # extract seasonality
  if(seasontype){
    # seasonal last states
    l = data.frame(extract_stan(object = object,pars = "s"))
    l = l[,n:n1];colnames(l) = paste0("s",1:m)
    states = cbind(states,l)

    # seasonal parameters
    pars = cbind(pars,data.frame(extract_stan(object = object,pars = "seasonal")))
    colnames(pars) = c("sigma","alpha","beta","phi","gamma")
  }

  # The previous data
  yh  =  matrix(0,nrow = draw,ncol = h)
  mu2 = rep(0,h)

  for(i in 1:draw){
    mu1 = class1(h = h,last.state = as.numeric(states[i,]),
                 trendtype = trendtype,seasontype = seasontype,
                 damped = damped,m = m,par = pars[i,])


    if(d1 > 0) mu2 = sum(xh*breg[i,])

    yh[i,] = rnorm(n = h,mean = mu1$mu + mu2,sd = mu1$var)
  }
  colnames(yh) = paste0("yh.",1:h)
  yh = as.data.frame(yh)
  return(yh)
}
#' Point forecast Kalman Filter procedure
#'
#' @param h number of predictiion
#' @param last.state a vector with the last states (l,b,s)
#' @param trendtype logical value for specify the trend
#' @param seasontype logical value for specify the  seasonality
#' @param damped logical value for specify a damped trend
#' @param m an integer for specify the periodicity
#' @param par a vreal vector with the states parameters
#'
#' @author Rob Hyndman.
#' @noRd
#'
class1 = function(h, last.state, trendtype, seasontype, damped, m, par){
  p = length(last.state);H = matrix(c(1, rep(0, p - 1)), nrow = 1)
  sigma2 = as.numeric(par["sigma"])^2

  # H matrix
  if (seasontype) H[1, p] = 1
  if (trendtype) {
    if (damped) H[1, 2] = as.numeric(par["phi"])
    else H[1, 2] = 1
  }

  # F matrix
  F =  matrix(0, p, p)
  F[1, 1] = 1

  if (trendtype) {
    if (damped)  F[1, 2] = F[2, 2] = as.numeric(par["phi"])
    else F[1, 2] = F[2, 2] = 1
  }
  if (seasontype) {
    F[p - m + 1, p] = 1
    F[(p - m + 2):p, (p - m + 1):(p - 1)] = diag(m - 1)
  }
  # G matrix
  G = matrix(0, nrow = p, ncol = 1)

  G[1, 1] = as.numeric(par["alpha"])
  if (trendtype == "A") G[2, 1] = as.numeric(par["beta"])
  if (seasontype == "A")G[3, 1] = as.numeric(par["gamma"])

  mu = numeric(h);Fj = diag(p);cj = numeric(h - 1)

  if (h > 1) {
    for (i in 1:(h - 1)){
      mu[i] = H %*% Fj %*% last.state
      cj[i] = H %*% Fj %*% G
      Fj = Fj %*% F
    }
    cj2 = cumsum(cj ^ 2)
    var = sigma2 * c(1, 1 + cj2)
  }
  else var = sigma2

  mu[h] = H %*% Fj %*% last.state

  return(list(mu = mu, var = sqrt(var), cj = cj))
}
