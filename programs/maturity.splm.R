b <- read.scsbio(2019, survey = "regular", sex = 1)
b  <- b[which(is.new.shell(b)), ]
b <- b[!is.na(b$carapace.width) & !is.na(b$chela.height) & (b$carapace.width <= 140), ]

x <- b$carapace.width
y <- b$chela.height

# Male morphometry initial values:
theta <- c(beta_immature = c(-2.660548, 1.3),
           beta_mature = c(-2.360548, 1.35), 
           log_sigma = -2.958989,
           log_sigma_kurtotic = 2,
           logit_p_kurtotic = -5.704126,
           p_alpha = -10.40064,
           p_beta = c(0.15893107, 0.01525848, 0.28906458),
           p_transition = c(58.79796, 101.05293),
           p_window = 1.448468)
      

fit.morphometry <- function(x, y, theta, fixed){
   # Negative log-likelihood function:
   loglike <- function(theta, x, y, fixed){
      if (!missing(fixed)) theta <- c(theta, fixed)
      return(-sum(morphometry(x, y, theta)$loglike))
   }
   
   # Maximum likelihood estimation:
   if (missing(fixed)) theta <- optim(theta, loglike, x = x, y = y, control = list(maxit = 5000))$par
   if (!missing(fixed)){
      theta <- optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(maxit = 5000))$par
      theta <- c(theta, fixed)
   }
   
   return(theta)
}

# Determine morphometric maturity probability.
morphometry <- function(x, y, theta, eval = FALSE, loglike = FALSE, classify = TRUE){
   # Horner's method:
   polyval <- function(x, p){
      y <- p[1] * rep(1, length(x))  
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
   
   # Return model values:
   v <- data.frame(mu_immature = mu_immature,
                   mu_mature = mu_mature,
                   p_mature = p_mature,
                   sigma = sigma,
                   sigma_kurtotic = sigma_kurtotic,
                   p_kurtotic = p_kurtotic)
      
   if (!missing(y)){
       # Mixture component densities:
       dimm <- (1-p_kurtotic) * dnorm(log(y), mu_immature, sigma) + p_kurtotic * dnorm(log(y), mu_immature, sigma + sigma_kurtotic) 
       dmat <- (1-p_kurtotic) * dnorm(log(y), mu_mature, sigma) + p_kurtotic * dnorm(log(y), mu_mature, sigma + sigma_kurtotic) 
       dmix <- (1-p_mature) * dimm + p_mature * dmat
       v$loglike <- log(dmix)
          
       # Posterior maturity probabilities:
       v$p_mature_posterior <- dmat / (dimm + dmat)          
   }
      
   return(v)
}

# Display morphometric plot:
plot.morphometry <- function(x, y, theta, log = TRUE){
   
   
}

