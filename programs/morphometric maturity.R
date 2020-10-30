theta 

theta <- c(alpha_immature  
           beta_immature
           alpha_mature
           beta_mature
           log_sigma_immature
           log_sigma_mature
           log_sigma_outlier
           logit_p_xp = log(60),
           logit_p_m  # maturity.
           logit_p_outlier

            
loglike <- function(x, y, theta){
   # Unconditional proportion of mature crab:
   logit_p <- 4 * theta[["p_m"]] * (x - theta[["p_xp"]])
   p <- 1 / (1 + exp(logit_p))  
   
   # Unconditional proportion of outliers:
   p_outlier <- 1 / (1 + exp(theta[["logit_p_outlier"]]))
   
   # Allometric relations:
   mu_immature = theta[["alpha_immature"]] + theta[["beta_immature"]] * log(x)
   mu_mature = theta[["alpha_mature"]] + theta[["beta_mature"]] * log(x)
   
   # Distribution errors:
   sigma_immature <- exp(theta[["log_sigma_immature"]])
   sigma_mature <- exp(theta[["log_sigma_mature"]])
   sigma_mature <- exp(theta[["log_sigma_outlier"]])
   
   # Gaussian mixture of kurtotic densities:
   v <- log((1-p) * ((1-p_outlier) * dnorm(log(y), mu_immature, sigma_immature) +                 
                       (p_outlier) * dnorm(log(y), mu_immature, sigma_immature + sigma_outlier)) +
              (p) * ((1-p_outlier) * dnorm(log(y), mu_mature, sigma_mature) +
                       (p_outlier) * dnorm(log(y), mu_mature, sigma_mature + sigma_outlier)))
         
   return(-sum(v))
}
  


   
