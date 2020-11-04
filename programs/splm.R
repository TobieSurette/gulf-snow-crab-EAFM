#' Smoothed Piecewise Linear Models
#' 
#' @description Functions for fitting 
#' 
#' @param x Numeric vector at which the function is to be evaluated.
#' @param theta Named parameter vector with an intercept parameter \code{alpha}, slope parameters \code{beta0, ..., betak}, 
#'              where \code{k} is the number of transition (i.e. break) points, \code{transition1, ..., transitionk} and 
#'              \code{window}, the transition window scale parameter (log-scale). \code{window} may also vary by transition 
#'              point.
#'              
#' @return A vector of the function evaluated at 'x'. If 'x' is not specified, then a function is returned.

#' @describeIn splm Smoothed piecewise-linear model definition and evaluation.
#' @export
splm <- function(x, theta){
   if (missing(x)){
      # Parse 'theta':
      alpha <- theta[sort(names(theta)[grep("alpha", names(theta))])]
      beta <- theta[sort(names(theta)[grep("beta", names(theta))])]
      transition <- theta[sort(names(theta)[grep("transition", names(theta))])]
      window <- exp(theta[sort(names(theta)[grep("window", names(theta))])])
      k <- length(transition)
      if (length(window) == 1) window <- rep(window, k)
      
      # Define 'splm' function:
      f <- function(x){
         y <- alpha + beta[1] * x
         for (i in 1:k) y <- y + window[i] * (beta[i+1] - beta[i]) * log(1 + exp((x - transition[i]) / window[i]))
         names(y) <- names(x)
         return(y)
      }
      return(f)   
   }else{
      return(splm(theta = theta)(x))
   }
}

