loglik_smooth_spline      <- function(fittedSmoothSpline) {
  ## Return the penalised loglikelihood
  #
  # FDA (Functional Data Analysis) uses the integrated squared second derivative as roughness, here we use the integrated second derivative.
  #
  # Calculate the (log) likelihood of a spline given the data used to fit the spline.
  # The likelihood consists of two main parts:
  #   1) (weighted) residuals sum of squares
  #   2) a penalty term   -   roughness is the integrated second derivative aka total curvature
  # The penalty term consists of a smoothing parameter lambda and a roughness measure of the spline
  #   J(s) = INTEGRAL s''(t) dt. 
  # The overall log likelihood is log L(s|x) = (y-s(x))'W(y-s(x)) + lambda J(s)
  
  s       <- fittedSmoothSpline
  x       <- fittedSmoothSpline$x
  y       <- fittedSmoothSpline$y
  w       <- fittedSmoothSpline$w
  yinput  <- fittedSmoothSpline$yin
  
  # weighted residuals sum of squares (weighted euclidian dist yinput to ypred) = loglikelihood without penality (but then forgets the smoothing part of the equation)
  wrss    <- sum(w * (yinput-y)^2)
  
  # smoothing parameter
  lambda  <- s$lambda
  
  # roughness score
  sDeriv  <- stats::smooth.spline(stats::predict(s, x, deriv=2))
  ab <- range(x, na.rm=TRUE)
  Js <- stats::integrate(function(x) stats::predict(sDeriv, x=x)$y,lower=ab[1], upper=ab[2], rel.tol=.Machine$double.eps^(1/8), stop.on.error=FALSE)$value
  
  # penalty term
  penalty <- -lambda * Js
  
  # penalised loglikelihood (as it's a sum, if it wasn't a log it would be multiplied)
  l <- (wrss + penalty)
  
  # The smoothing spline estimate of a function is defined to be the minimizer of l
  return(l)
}


#' Calculate the Akaike Information Criterion for a smooth.spline
#'
#' Calculate the Akaike Information Criterion (\emph{AIC}) for a fitted \code{\link[stats]{smooth.spline}}. The smaller the AIC, the better the spline fit.
#'
#' @param fittedSmoothSpline A fitted \code{\link[stats]{smooth.spline}}
#'
#' @return The AIC value.
AIC_smooth_spline         <- function(fittedSmoothSpline) {
  # AIC = -2*logLik + k*npar
  # for AIC k=2
  # npar = df
  
  df = fittedSmoothSpline$df
  AIC = -2 * -loglik_smooth_spline(fittedSmoothSpline) + 2 * df  
  return(AIC)
}


#' Calculate the Akaike Information Criterion Corrected for small observation numbers for a smooth.spline
#'
#' Calculate the Akaike Information Criterion Corrected for small observation numbers (\emph{AICc}) for a fitted \code{\link[stats]{smooth.spline}}. The smaller the AICc, the better the spline fit.
#'
#' @param fittedSmoothSpline A fitted \code{\link[stats]{smooth.spline}}
#'
#' @return The AICc value.
AICc_smooth_spline        <- function(fittedSmoothSpline) {
  df    = fittedSmoothSpline$df
  nobs  = length(fittedSmoothSpline$yin)
  AICc  = AIC_smooth_spline(fittedSmoothSpline) + 2 * df * (df + 1) / (nobs - df - 1)
  return(AICc)
}
