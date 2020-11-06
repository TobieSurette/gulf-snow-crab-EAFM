#' Analyze Morphometric Data
#' 
#' @description Functions to analyze and visualize morphometric data.
#' 
#' @param x Predictor variable.
#' @param y Response variable.
#' @param z Binary classification variable (optional).
#' @param theta Named parameter vector.

#' @export
morphometry <- function(x, ...) UseMethod("morphometry")

#' @describeIn morphometry Perform morphometric regression analysis.
#' @export
morphometry.default <- function(x, y, species, kurtotic = FALSE, ...){
   # General case:
   m <- lm(log(y) ~ log(x))
   theta <- coef(m)
   names(theta) <- c("alpha", "beta")
   theta["alpha"] <- exp(theta["alpha"])
   
   return(theta)
}

#' @describeIn morphometry Perform morphometric regression analysis for snow crab biological data.
#' @export
morphometry.scsbio <- function(x, y, z, theta, sex, fit = TRUE, discrete = FALSE){
   if ("scsbio" %in% names(x)){
      x <- x$carapace.width
      y <- y$chela.height
      sex <- x$sex
   }
   
   # Horner's method for polynomial evaluation::
   polyval <- function(x, p){
      p <- rev(p) # Intercept is last parameter:
      y <- p[[1]] * rep(1, length(x))  
      if (length(p) > 1) for (i in 2:length(p)) y <- y * x + p[i]
      return(y)
   }   
   
   # Calculate morphometric means:
   beta_immature <- theta[sort(names(theta)[grep("beta_immature", names(theta))])] # Immature coefficients.
   beta_mature   <- theta[sort(names(theta)[grep("beta_mature", names(theta))])]   # Mature coefficients.
   mu_immature <- polyval(log(x), beta_immature) # Immature mean.
   mu_mature   <- polyval(log(x), beta_mature)   # Mature mean.
   
   # Calculate unconditional maturity proportions:
   eta <- theta[grep("^p_", names(theta))] # Maturity proportion parameters.
   logit_p_mature <- splm(x, eta)
   p_mature <- 1 / (1 + exp(-logit_p_mature))
   
   # Error parameters:
   sigma <- exp(theta[["log_sigma"]])
   sigma_kurtotic <- exp(theta[["log_sigma_kurtotic"]])
   
   # Proportion of kurtotic distributional component:
   p_kurtotic <- 1 / (1 + exp(-theta[["logit_p_kurtotic"]]))  
   
   # Compile results:
   v <- data.frame(mu_immature = mu_immature,
                   mu_mature = mu_mature,
                   p_mature = p_mature,
                   sigma = sigma,
                   sigma_kurtotic = sigma_kurtotic,
                   p_kurtotic = p_kurtotic)
      
   # Calculate log-likelihood and posterior classification probabilities:
   if (!missing(y)){
       # Mixture component densities:
       if (discrete){
          dimm <- (1-p_kurtotic) * dnorm(log(y), mu_immature, sigma) + p_kurtotic * dnorm(log(y), mu_immature, sigma + sigma_kurtotic) 
          dmat <- (1-p_kurtotic) * dnorm(log(y), mu_mature,   sigma) + p_kurtotic * dnorm(log(y), mu_mature,   sigma + sigma_kurtotic) 
       }else{
          dimm <- ((1-p_kurtotic) * pnorm(log(y+0.5), mu_immature, sigma) + p_kurtotic * pnorm(log(y+0.5), mu_immature, sigma + sigma_kurtotic)) -  
                  ((1-p_kurtotic) * pnorm(log(y-0.5), mu_immature, sigma) + p_kurtotic * pnorm(log(y-0.5), mu_immature, sigma + sigma_kurtotic))
          dmat <- ((1-p_kurtotic) * pnorm(log(y+0.5), mu_mature,   sigma) + p_kurtotic * pnorm(log(y+0.5), mu_mature,   sigma + sigma_kurtotic)) -          
                  ((1-p_kurtotic) * pnorm(log(y-0.5), mu_mature,   sigma) + p_kurtotic * pnorm(log(y-0.5), mu_mature,   sigma + sigma_kurtotic)) 
       }

       # Initialize log-likelihood:
       v$loglike <- rep(0, length(x)) 

       # Likelihood for known classification:
       if (!missing(z)){
          ix <-!is.na(z) 
          v$loglike[ix] <- (1-z[ix]) * log(1-p_mature[ix]) + z[ix] * log(p_mature[ix]) # Bernouilli likelihood.
          iy <- which(ix & (z == 0) & !is.na(x) & !is.na(y))
          v$loglike[iy]  <- v$loglike[iy] + log(dimm[iy]) # Gaussian immature.
          iy <- which(ix & (z == 1) & !is.na(x) & !is.na(y))
          v$loglike[iy]  <- v$loglike[iy] + log(dmat[iy]) # Gaussian mature.
       }
       ix <- is.na(x) | is.na(y)      # Missing morphometry.
       v$loglike[!ix]  <- v$loglike[!ix] + log((1-p_mature[!ix]) * dimm[!ix] + p_mature[!ix] * dmat[!ix])  # Add mixture contribution.
      
       # Posterior maturity probabilities:
       v$p_mature_posterior <- p_mature * dmat / ((1-p_mature) * dimm + p_mature * dmat)    
       v$p_mature_posterior[is.na(v$p_mature_posterior)] <- p_mature[is.na(v$p_mature_posterior)] 
   }
      
   return(v)
}

# Fit morphometric model:
fit.morphometry.scsbio <- function(x, y, z, sex, theta, discrete = FALSE){
   
   if (!missing(sex) & missing(theta)){
      if (sex == 1){
         # Male morphometry initial values:
         theta <- c(beta_immature = c(-2.03, 1.116, -0.06026, 0.0114), # Log-scale immature morphometric coefficients.
                    beta_mature = c(-2.858, 1.312),  # Log-scale mature morphometric coefficients.
                    log_sigma = -3.3,                # Log-scale standard error.
                    log_sigma_kurtotic = 0,          # Log-scale extra standard error for kurtotic observations.
                    logit_p_kurtotic = -2,           # Logit-scale proportion of kurtotic observations.
                    p_alpha = -11,                   # Logit-scale splm intercept parameter for mature proportions.
                    p_beta = c(0.25, 0.015, 0.25),   # Logit-scale splm slope parameters for mature proportions.
                    p_transition = c(45, 95),        # Logit-scale transition parameters for mature proportions.
                    p_window = 2.0)                  # Logit-scale splm window width parameter(s) for mature proportions.
      }
      if (sex == 2){
         # Male morphometry initial values:
         theta <- c(beta_immature = c(-2.72, 1.228), # Log-scale immature morphometric coefficients.
                    beta_mature = c(-2.80, 1.30),    # Log-scale mature morphometric coefficients.
                    log_sigma = -3,                  # Log-scale standard error.
                    log_sigma_kurtotic = 2,          # Log-scale extra standard error for kurtotic observations.
                    logit_p_kurtotic = -5.7,         # Logit-scale proportion of kurtotic observations.
                    p_alpha = -10.4,                 # Logit-scale splm intercept parameter for mature proportions.
                    p_beta = c(0.16, 0.015, 0.29),   # Logit-scale splm slope parameters for mature proportions.
                    p_transition = c(58.8, 101.1),   # Logit-scale transition parameters for mature proportions.
                    p_window = 1.45)                 # Logit-scale splm window width parameter(s) for mature proportions.
      }
   }
   
   # Negative log-likelihood function:
   loglike <- function(theta, x, y, z, fixed, discrete = FALSE){
      if (!missing(fixed)) theta <- c(theta, fixed)
      return(-sum(morphometry.scsbio(x, y, z = z, theta = theta, discrete = discrete)$loglike))
   }
   
   # Fit proportions
   cat("Fitting mature proportion parameters.\n")
   fixed <- theta[-grep("^p_", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, fixed = fixed, discrete = discrete, control = list(trace = 0, maxit = 1000))$par
   theta <- c(theta, fixed)
   
   # Fit kurtotic parameters:
   cat("Fitting kurtosis parameters.\n")
   fixed <- theta[-grep("kurtotic", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, fixed = fixed, discrete = discrete, control = list(trace = 0, maxit = 500))$par
   theta <- c(theta, fixed)   
   
   # Fit immature regression:
   cat("Fitting immature regression coefficients.\n")
   fixed <- theta[-grep("immature", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, fixed = fixed, discrete = discrete, control = list(trace = 0, maxit = 500))$par
   theta <- c(theta, fixed)      

   # Fit non-regression coefficients:
   cat("Fitting non-regression coefficients.\n")
   fixed <- theta[grep("mature", names(theta))]
   theta <- theta[setdiff(names(theta), names(fixed))]
   theta <- optim(theta, loglike, x = x, y = y, z = z, fixed = fixed, discrete = discrete, control = list(trace = 0, maxit = 2000))$par
   theta <- c(theta, fixed)   
   
   # Fit immature regression:
   cat("Fitting complete model.\n")
   theta <- optim(theta, loglike, x = x, y = y, z = z, discrete = discrete, control = list(trace = 3, maxit = 5000))$par

   return(theta)
}

#' Display morphometric plot:
#' @export 
plot.morphometry.scsbio <- function(x, y, theta, xlim = c(10, 140), log = TRUE, title, discrete = FALSE, ...){
   v <- morphometry.scsbio(x, y, theta = theta)
   x0 <- seq(xlim[1], xlim[2], len = 1000)
   v0 <- morphometry.scsbio(x0, theta = theta)
   
   # Prepare plotting area:
   layout(rbind(0, 0, cbind(0, kronecker(c(1,1,2:4), matrix(1, nrow = 5, ncol = 5)), 0), 0, 0))
   par(mar = c(0, 0, 0, 0))
   
   # Plot data and regression curves:
   plot(x, y, type = "n", xlab = "", ylab = "", xlim = xlim, ylim = c(0, 40), cex.lab = 1.25, xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   if (!discrete){
      points(x[v$p_mature_posterior >= 0.5], y[v$p_mature_posterior >= 0.5], pch = 19, col = "deepskyblue1", cex = 0.1)
      points(x[v$p_mature_posterior < 0.5], y[v$p_mature_posterior < 0.5], pch = 19, col = "chartreuse2", cex = 0.1)
   }else{
      points(jitter(x[v$p_mature_posterior >= 0.5], amount = 0.5), 
             jitter(y[v$p_mature_posterior >= 0.5], amount = 0.5), pch = 19, col = "deepskyblue1", cex = 0.1)
      points(jitter(x[v$p_mature_posterior < 0.5], amount = 0.5), 
             jitter(y[v$p_mature_posterior < 0.5], amount = 0.5), pch = 19, col = "chartreuse2", cex = 0.1)
   }
   lines(x0, exp(v0$mu_immature), lwd = 2, col = "chartreuse3")
   lines(x0, exp(v0$mu_mature), lwd = 2, col = "deepskyblue3")
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   mtext("Chela height (mm)", 2, 2.5, cex = 1.0)
   if (!missing(title)) mtext(title, 3, 1.0, cex = 1.25)
   legend("topleft",
          legend = c("Mature", "Immature"),
          col = c("deepskyblue3", "chartreuse3"),
          lwd = 2)
   box()

   # Plot maturity proportions:
   cw <- round(x)
   res <- aggregate(list(k = (v$p_mature_posterior  >= 0.5)), by = list(cw = cw), sum)
   res$n <- aggregate(list(n = (v$p_mature_posterior  >= 0.5)), by = list(cw = cw), length)$n
   res$p_mature <-  res$k / res$n
   plot(xlim, c(0, 1.1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   gbarplot(res$p_mature, res$cw, col = "grey90", width = 1, add = TRUE)
   lines(x0, v0$p_mature, col = "red", lwd = 2)
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   legend("topleft",
          legend = c("model", "classification"),
          col = c("red", "grey30"),
          pch = c(NA, 22),
          pt.bg = c(NA, "grey60"),
          pt.cex = c(1,2.5),
          lwd = c(2, 0))
   box()
   mtext("Mature Proportion", 2, 2.5, cex = 1.0)
   
   # Residual plots
   if (!discrete){
      r_imm <- (log(y[v$p_mature_posterior < 0.5]) - v$mu_immature[v$p_mature_posterior < 0.5]) / exp(theta[["log_sigma"]])
      r_mat <- (log(y[v$p_mature_posterior >= 0.5]) - v$mu_mature[v$p_mature_posterior >= 0.5]) / exp(theta[["log_sigma"]]) 
   }else{
      r_imm <- (log(jitter(y[v$p_mature_posterior < 0.5], amount = 0.5)) - v$mu_immature[v$p_mature_posterior < 0.5]) / exp(theta[["log_sigma"]])
      r_mat <- (log(jitter(y[v$p_mature_posterior >= 0.5], amount = 0.5)) - v$mu_mature[v$p_mature_posterior >= 0.5]) / exp(theta[["log_sigma"]])       
   }   
   res$p_model <- morphometry.scsbio(res$cw, theta = theta)$p_mature
   res$logit_p_model <- log(res$p_model / (1-res$p_model))
   res$logit_p <- log(res$p_mature / (1-res$p_mature)) 
   res$residuals_logit_p <- res$logit_p - res$logit_p_model
   
   # Mature residual plot:
   if (!discrete){
      plot(x[v$p_mature_posterior >= 0.5], r_mat, pch = 19, xlab = "", ylab = "", xlim = xlim, col = "deepskyblue1", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }else{
      plot(jitter(x[v$p_mature_posterior >= 0.5], amount = 0.5), r_mat, pch = 19, xlab = "", ylab = "", xlim = xlim, col = "deepskyblue1", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }
   grid()
   abline(0, 0, col = "deepskyblue1", lwd = 2)  
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   mtext("Mature residuals", 2, 2.5, cex = 1.0)
   box()
   
   # Immature residual plot:
    if (!discrete){
       plot(x[v$p_mature_posterior < 0.5], r_imm,  xlab = "",  ylab = "", xlim = xlim, pch = 19, col = "chartreuse2", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }else{
       plot(jitter(x[v$p_mature_posterior < 0.5], amount = 0.5), r_imm,  xlab = "",  ylab = "", xlim = xlim, pch = 19, col = "chartreuse2", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }
   grid()
   abline(0, 0, col = "chartreuse3", lwd = 2)
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   mtext("Immature residuals", 2, 2.5, cex = 1.0)
   box()
   
   axis(1)
   mtext("Carapace width (mm)", 1, 2.5, cex = 1.15)   
}



