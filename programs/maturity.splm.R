b <- read.scsbio(2019, survey = "regular", sex = 1)
b  <- b[which(is.new.shell(b)), ]
b <- b[!is.na(b$carapace.width) & !is.na(b$chela.height) & (b$carapace.width <= 140), ]

x <- b$carapace.width
y <- b$chela.height

# Fit morphometric model:
fit.morphometry <- function(...){
   # Negative log-likelihood function:
   loglike <- function(theta, x, y, fixed){
      if (!missing(fixed)) theta <- c(theta, fixed)
      return(-sum(morphometry(x, y, theta)$loglike))
   }
   
   # Fit proportions
   
   # Fit standard errors:
   
   # Fit immature regression coefficients (i.e. quadratic or cubic polynomial):
   
   # Fit proportions

   # Fit model:
   args <- list(...)
   theta <- optim(theta = args$theta, loglike, ..., control = list(trace = 3, maxit = 5000))$par
   if ("fixed" %in% names(args)) theta <- c(theta,  args$fixed)
   
   return(theta)
}

# Evaluate morphometric model:
# @param x Predictor variable.
# @param y Response variable.
# @param z Binary classification variable (optional).
# @param theta Named parameter vector.
#' @export
morphometry.scs <- function(x, y, z, sex, theta){
   if (!missing(sex)){
      if (sex == 1){
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
       dimm <- (1-p_kurtotic) * dnorm(log(y), mu_immature, sigma) + p_kurtotic * dnorm(log(y), mu_immature, sigma + sigma_kurtotic) 
       dmat <- (1-p_kurtotic) * dnorm(log(y), mu_mature, sigma) + p_kurtotic * dnorm(log(y), mu_mature, sigma + sigma_kurtotic) 

       
       v$loglike <- rep(0, length(x)) # Initialize log-likelihood.
       ix <- is.na(x) | is.na(y)      # Missing morphometry.
       
       # Likelihood for known classification:
       if (!missing(z)){
          ix <-!is.na(z) 
          v$loglike[ix] <- (1-z[ix]) * log(1-p_mature[ix]) + z[ix] * log(p_mature[ix]) # Bernouilli likelihood.
          iy <- which(ix & (z == 0) & !is.na(x) & !is.na(y))
          v$loglike[iy]  <- v$loglike[iy] + log(dimm[iy]) # Gaussian immature.
          iy <- which(ix & (z == 1) & !is.na(x) & !is.na(y))
          v$loglike[iy]  <- v$loglike[iy] + log(dmat[iy]) # Gaussian mature.
       }
       v$loglike[!ix]  <- v$loglike[!ix]  + log((1-p_mature[!ix]) * dimm[!ix] + p_mature[!ix] * dmat[!ix])  # Add mixture contribution.
      
       # Posterior maturity probabilities:
       v$p_mature_posterior <- dmat / (dimm + dmat)          
   }
      
   return(v)
}

# Display morphometric plot:
plot.morphometry.scs <- function(x, y, theta, log = TRUE){
   xlim <- c(40, 140)
   v <- morphometry(x, y, theta)
   x0 <- seq(xlim[1], xlim[2], len = 1000)
   v0 <- morphometry(x0, theta = theta)
   
   clg()
   gdevice(layout = 1:4)
   plot(x, y, type = "n", xlab = "", ylab = "", xlim = xlim, ylim = c(0, 40), cex.lab = 1.25, xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   col <- colorRamp(c("white", "blue"))(0.25)[1, ] / 255
   col <- rgb(col[1], col[2], col[3])
   points(x[v$p_mature_posterior >= 0.5], y[v$p_mature_posterior >= 0.5], pch = 21, bg = col, cex = 0.25)
   col <- colorRamp(c("white", "green"))(0.25)[1, ] / 255
   col <- rgb(col[1], col[2], col[3])
   points(x[v$p_mature_posterior < 0.5], y[v$p_mature_posterior < 0.5], pch = 21, bg = col, cex = 0.25)
   lines(x0, exp(v0$mu_immature), lwd = 2, col = "green")
   lines(x0, exp(v0$mu_mature), lwd = 2, col = "blue")
   mtext("Chela height (mm)", 2, 2.5, cex = 1.25)
   box()
   
   # Plot maturity proportions:
   cw <- round(x)
   res <- aggregate(list(k = (v$p_mature_posterior  >= 0.5)), by = list(cw = cw), sum)
   res$n <- aggregate(list(n = (v$p_mature_posterior  >= 0.5)), by = list(cw = cw), length)$n
   res$p_mature <-  res$k / res$n
   plot(xlim, c(0, 1.1), type = "n", xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   gbarplot(res$p_mature, res$cw, col = "grey90", width = 1, add = TRUE)
   lines(x0, v0$p_mature, col = "red", lwd = 2)
   box()
   mtext("Mature Proportion", 2, 2.5, cex = 1.25)
   
   # Residual plots
   r_imm <- (log(y[v$p_mature_posterior < 0.5]) - v$mu_immature[v$p_mature_posterior < 0.5]) / exp(theta[["log_sigma"]])
   r_mat <- (log(y[v$p_mature_posterior >= 0.5]) - v$mu_mature[v$p_mature_posterior >= 0.5]) / exp(theta[["log_sigma"]]) 
   res$p_mature_model <- morphometry(res$cw, theta = theta)$p_mature
   res$p_mature_model <- log(res$p_mature_model / (1-res$p_mature_model))
   res$logit_p <- log(res$p_mature / (1-res$p_mature)) 
   res$logit_p_mature_model <- log(res$p_mature_model / (1-res$p_mature_model)) 
   res$residuals_p <- res$logit_p - res$logit_p_model
   
   plot(x[v$p_mature_posterior < 0.5], r_imm, pch = 21, xlim = xlim, bg = "grey", cex = 0.25, xaxt = "n", xaxs = "i")
   grid()
   abline(0, 0, col = "red", lwd = 2)
   mtext("Immature residuals", 2, 2.5, cex = 1.25)
   
   plot(x[v$p_mature_posterior >= 0.5], r_mat, pch = 21, xlim = xlim, bg = "grey", cex = 0.25, xaxt = "n", xaxs = "i")
   grid()
   abline(0, 0, col = "red", lwd = 2)  
   mtext("Mature residuals", 2, 2.5, cex = 1.25)
   
   axis(1)
   mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)   
}

