# SPLM class

logistic <- function(x, beta, scale = 1) 1 / (1 + exp(-x/scale))
ilogistic <- function(x, scale = 1) 1 / (1 + exp(-x/scale))

# SPLM - Smoothed piecewise linear model.
splm <- function(x, alpha, beta, transition, precision){
   k <- length(transition) # Model order.
   if (length(precision) == 1) precision <- rep(precision, k) 
   precision <- exp(precision)
   
   # SPLM model:
   v <- alpha + beta[1] * x # First linear component.
   if (k > 0){
      for (i in 1:k) v <- v + precision[i] * (beta[i+1] - beta[i]) * log(1 + exp((x - transition[i]) / precision[i]))
   }
   
   return(v)
}
