splm <- function(x, alpha, beta, transition, precision){
   # SPLM - Smoothed piecewise linear model.
   
   # Model order:
   k <- length(transition)
   
   # Expand precision:
   precision <- exp(precision)
   if (length(precision) == 1) precision <- rep(precision, k)

   # SPLM model:
   v <- alpha + beta[1] * x
   for (i in 1:k) v <- v + precision[i] * (beta[i+1] - beta[i]) * log(1 + exp((x - transition[i]) / precision[i]))

   return(v)
}
