library(gulf)

x0 <- 100
y0 <- 0
m <- 0.05
a <- 0.0022
b <- 0.115

clg()

t <- seq(40, 150, length = 1000)
p <- a*(t-x0)^3 + b*(t-x0)^2 + m*(t-x0) + y0

mu <- function(x, theta, scale = 0){
   # Scaling factor:
   s <- exp(scale)
   
   # Pygmy male component:
   d <- (x - theta[["x0"]])
   v <- s*(theta[["a"]]*(d^3) + theta[["b"]]*(d^2))
   v[(x >= theta[["x0"]])] <- 0
   
   # Pure-logistic component for males:
   v <- v + theta[["m"]] * (x - theta[["x0"]]) + theta[["y0"]]

   return(v)
}

loglike <- function(theta, xx, yy){
   v <- dnorm(yy, mu(xx, theta), exp(theta[["log.sigma"]]), log = TRUE) 
   return(-sum(v))
}

xx <- res$length[is.finite(res$logit)]
yy <- res$logit[is.finite(res$logit)]

theta <- c(x0 = 100, y0 = 0, m = 0.1, a = 0, b = 0, log.sigma = 0)
loglike(theta, xx, yy)
theta <- optim(theta, loglike, xx = xx, yy = yy, control = list(trace = 3, maxit = 2000))$par

clg()
#plot(t, mu(t, theta), type = "n")
v <- mu(t, theta, scale = 0)
plot(t, exp(v) / (1 + exp(v)), type = "n", xlab = "Carapace width (mm)", ylab = "Proportion", xlim = c(40, 140), ylim = c(0, 1), yaxs = "i", xaxs = "i")
index <- is.finite(res$logit)
dbarplot(res$p, res$length, add = TRUE, width = 1)
lines(t, exp(v) / (1 + exp(v)), lwd = 2, col = "red")
box()

clg()
#plot(t, mu(t, theta), type = "n")
v <- mu(t, theta, scale = +0.15)
plot(t, exp(v) / (1 + exp(v)), type = "n", xlab = "Carapace width (mm)", ylab = "Proportion", xlim = c(40, 140), ylim = c(0, 1), yaxs = "i", xaxs = "i")
index <- is.finite(res$logit)
dbarplot(res$p, res$length, add = TRUE, width = 1)
for (scale in seq(-2.5, 0.25, by = 0.01)){
   v <- mu(t, theta, scale = scale)
   lines(t, exp(v) / (1 + exp(v)), lwd = 2, col = "red")
}
box()

# TMB code:

   PARAMETER_VECTOR(pygmy_effect);
   PARAMETER(log_sigma_pygmy);
   Type sigma_pygmy = exp(log_sigma_pygmy); 
   
   for (int j = 0; j < n_station; j++){
      res -= dnorm(pygmy_effect[j], Type(0.0), sigma_pygmy, true);
   }   

   # Pygmy male component:
   d <- (x - x0)
   v <- s*(a*(d^3) +b*(d^2))
   v[(x >= x0)] <- 0
   
   # Pure-logistic component for males:
   logit_p <- v + m * (x - x0) + y0
 
   # Proportion of mature crab:
 
   logit_p_mature = eta_alpha + eta_beta[0] * cw[i] + 
   
  
                   eta_precision * (eta_beta[1] - eta_beta[0]) * log(1 + exp((cw[i] - eta_transition[0]) / eta_precision)) + 
                   eta_precision * (eta_beta[2] - eta_beta[1]) * log(1 + exp((cw[i] - eta_transition[1]) / eta_precision));
  p_mature = Type(1) / (Type(1) + exp(-logit_p_mature));

