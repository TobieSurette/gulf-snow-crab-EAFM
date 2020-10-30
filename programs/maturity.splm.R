maturity.splm <- function(x, y, p){
   # MATURITY.SPLM -
   
   # Define list of parameter values:
   if (missing(p)){
      # SCS 2017 fit:
      p <- list(log_alpha_immature = -2.660548,
                log_alpha_mature = -3.057006,
                log_beta_immature = c(1.211642, 1.515985),
                log_beta_mature = 1.355044,
                log_precision_immature = -1.644497,
                transition_immature = 4.596061,
                log_sigma = -2.958989,
                log_sigma_outlier = 2,
                alpha_outlier = -5.704126,
                eta_alpha = -10.40064,
                eta_beta = c(0.15893107, 0.01525848, 0.28906458),
                eta_transition = c(58.79796, 101.05293),
                log_eta_precision = 1.448468)
   }
              
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

   # Calculate predicted allometric means:
   mu_imm <- splm(log(x), p$log_alpha_immature, p$log_beta_immature, p$transition_immature, p$log_precision_immature)
   mu_mat <- p$log_alpha_mature + p$log_beta_mature * log(x)    
   
   # Determine maturity and ouliers:
   logit_pp <- splm(x, p$eta_alpha, p$eta_beta, p$eta_transition, p$log_eta_precision)
   p_mature <- 1 / (1 + exp(-logit_pp))
   
   # Error parameters:
   sigma <- exp(p$log_sigma)
   sigma_outlier <- exp(p$log_sigma_outlier)
   p_outlier <- 1 / (1 + exp(-p$alpha_outlier))  
   
   # Component densities:
   dout <- (1-p_mature) * p_outlier * dnorm(log(y), mu_imm, sigma + sigma_outlier, FALSE) +
               p_mature * p_outlier * dnorm(log(y), mu_mat, sigma + sigma_outlier, FALSE) 
   dimm <- (1-p_mature) * (1-p_outlier) * dnorm(log(y), mu_imm, sigma, FALSE)                 
   dmat <-     p_mature * (1-p_outlier) * dnorm(log(y), mu_mat, sigma, FALSE)  

   # Classify observations:
   pmat <- dmat / (dimm + dmat + dout)
   pimm <- dimm / (dimm + dmat + dout)
   pout <- dout / (dimm + dmat + dout)  
   
   # Determine which are mature:
   class <- apply(data.frame(pout, pimm, pmat), 1, which.max)
   index <- unlist(lapply(class, length)) == 0
   class[index] <- NA
   class[which(class == 1)] <- NA
   class[x <= 40] <- 2
   
   index <- which(is.na(class) & !is.na(x))
   logit_pp <- splm(x[index], p$eta_alpha, p$eta_beta, p$eta_transition, p$log_eta_precision)
   p_mature <- 1 / (1 + exp(-logit_pp))   
   class[index[p_mature >= 0.5]] <- 3
   class[index[p_mature < 0.5]] <- 2
   
   return(class == 3)
}
