year <- 1987

# - Add samplers to plots
# - SPL model for the proportions, two inflection points, but start with a single one.
# - Add discrete options for data < 1998
# - Compile model parameters for each year and compare
# - Consider formulating TMB version for a hierarchical analysis.
# - Consider treating immatures as quadratic or cubic
# - Plot proportion of matures versus cw.
# - Output summary plot for each year
# - Add more summary stats on yearly plots: samplers, # tows in survey, # crab, densities
# - Once a reliable classification is is hand, check for spatial patterns, depth patterns.

# Load dataset:
b <- read.scsbio(year, survey = "regular", sex = 1)
b  <- b[which(is.new.shell(b)), ]

# Remove empty data:
x <- log(b$carapace.width)
y <- log(b$chela.height)
index <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
x <- exp(x[index])
y <- exp(y[index])

# Define model log-likelihood function:
loglike <- function(theta, x, y, fixed, discrete = FALSE){
   if (!missing(fixed)) theta <- c(theta, fixed)
   
   logx <- log(x)
   
   # Unconditional proportion of mature crab:
   logit_p <- theta[["logit_p_beta0"]] + 
              theta[["logit_p_beta1"]] * logx + 
              theta[["logit_p_beta2"]] * (logx^2) + 
              theta[["logit_p_beta3"]] * (logx^3)
   p <- 1 / (1 + exp(-logit_p))  
   
   # Unconditional proportion of outliers:
   p_outlier <- 1 / (1 + exp(-theta[["logit_p_outlier"]]))
   
   # Allometric relations:
   mu_immature = theta[["beta0_immature"]] + theta[["beta1_immature"]] * log(x)
   mu_mature = theta[["beta0_mature"]] + theta[["beta1_mature"]] * log(x)
   
   # Distribution errors:
   sigma <- exp(theta[["log_sigma"]])
   sigma_outlier  <- exp(theta[["log_sigma_outlier"]])
   
   # Gaussian mixture of kurtotic densities:
   v <- log((1-p) * ((1-p_outlier) * dnorm(log(y), mu_immature, sigma) +                 
                       (p_outlier) * dnorm(log(y), mu_immature, sigma + sigma_outlier)) +
              (p) * ((1-p_outlier) * dnorm(log(y), mu_mature, sigma) +
                       (p_outlier) * dnorm(log(y), mu_mature, sigma + sigma_outlier)))
         
   return(-sum(v))
}

maturity.probability <- function(x, y, theta){
   logx <- log(x)
   # Unconditional proportion of mature crab:
   logit_p <- theta[["logit_p_beta0"]] + 
              theta[["logit_p_beta1"]] * logx + 
              theta[["logit_p_beta2"]] * (logx^2) + 
              theta[["logit_p_beta3"]] * (logx^3)
   p <- 1 / (1 + exp(-logit_p))  
   
   # Unconditional proportion of outliers:
   p_outlier <- 1 / (1 + exp(-theta[["logit_p_outlier"]]))
   
   # Allometric relations:
   mu_immature = theta[["beta0_immature"]] + theta[["beta1_immature"]] * log(x)
   mu_mature = theta[["beta0_mature"]] + theta[["beta1_mature"]] * log(x)
   
   # Distribution errors:
   sigma <- exp(theta[["log_sigma"]])
   sigma_outlier  <- exp(theta[["log_sigma_outlier"]])
   
   # Outlier, immature and mature posterior probabilities:
   p_imm <- (1-p) * ((1-p_outlier) * dnorm(log(y), mu_immature, sigma) + (p_outlier) * dnorm(log(y), mu_immature, sigma + sigma_outlier)) 
   p_mat <- (p) * ((1-p_outlier) * dnorm(log(y), mu_mature, sigma) + (p_outlier) * dnorm(log(y), mu_mature, sigma + sigma_outlier))
   p_mat <- p_mat / (p_imm + p_mat) 
                    
   return(p_mat)
}

# Initial parameter values:
theta <- c(beta0_immature = -2.717,   # Immature log-scale linear intercept.
           beta1_immature = 1.228,    # Immature log-scale linear slope.
           beta0_mature = -2.80,      # Mature log-scale linear intercept.
           beta1_mature = 1.30,       # Mature log-scale linear slope.
           log_sigma = -3,            # Log-scale standard error.
           log_sigma_outlier = 2,     # Log-scale outlier standard error.
           logit_p_beta0 = -21.75904, # Logit-scale proportion of mature constant coeffcient.
           logit_p_beta1 = 4.79,      # Logit-scale proportion of mature linear coefficient.
           logit_p_beta2 = 0,         # Logit-scale proportion of mature quadratic coefficient.
           logit_p_beta3 = 0,         # Logit-scale proportion of mature cubic coefficient.
           logit_p_outlier = -5)      # Logit-scale proportion of outliers.

loglike(theta, x, y, fixed = fixed)

# Fit error parameters:
fixed <- theta[setdiff(names(theta), c("log_sigma"))]
theta <- theta[setdiff(names(theta), names(fixed))]
theta <- c(optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(trace = 2, maxit = 5000))$par, fixed)

# Fit maturity proportion parameters:
fixed <- theta[setdiff(names(theta), c("logit_p_beta0", "logit_p_beta1"))]
theta <- theta[setdiff(names(theta), names(fixed))]
theta <- c(optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(trace = 2, maxit = 5000))$par, fixed)

# Fit maturity proportion parameters:
fixed <- theta[setdiff(names(theta), c("logit_p_beta0", "logit_p_beta1", "logit_p_beta2", "logit_p_beta3"))]
theta <- theta[setdiff(names(theta), names(fixed))]
theta <- c(optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(trace = 2, maxit = 5000))$par, fixed)

# Fit outlier proportion parameter:
fixed <- theta[setdiff(names(theta), c("logit_p_outlier"))]
theta <- theta[setdiff(names(theta), names(fixed))]
theta <- c(optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(trace = 2, maxit = 5000))$par, fixed)

# Fit model:
fixed <- theta["log_sigma_outlier"]
theta <- theta[setdiff(names(theta), names(fixed))]
theta <- c(optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(trace = 2, maxit = 5000))$par, fixed)

# Fit model:
theta <- optim(theta, loglike, x = x, y = y, fixed = fixed, control = list(trace = 2, maxit = 5000))$par

# Plot model:
#gdevice(width = 5, height = 5)
m <- kronecker(1:3, matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0, 0, 0, 0))

x0 <- seq(exp(0), exp(5), len = 1000)
mu_immature = theta[["beta0_immature"]] + theta[["beta1_immature"]] * log(x0)
mu_mature = theta[["beta0_mature"]] + theta[["beta1_mature"]] * log(x0)
   
plot(c(3.4, 5), c(1.5, 4), type ="n", pch = 21, cex = 0.6, xaxs = "i", yaxs = "i", xlab = "", ylab = "", xaxt = "n")
grid()
p <- maturity.probability(x, y, theta)
group = (p >= 0.5) + 1 
bg <- c("blue", "green")[group]
if (year > 1997){
   points(log(x), log(y), pch = 21, cex = 0.3, col = bg, bg  = bg, xlim = c(3.4, 5), ylim = c(1.5, 4), xaxs = "i", yaxs = "i")
}else{
   points(log(x+runif(length(x))-0.5), log(y+runif(length(x))-0.5), pch = 21, cex = 0.3, col = bg, bg  = bg, xlim = c(3.4, 5), ylim = c(1.5, 4), xaxs = "i", yaxs = "i")
}
lines(log(x0), mu_immature, col = "red")
lines(log(x0), mu_mature, col = "red")
mtext("ln(chela height)", 2, 2.5, cex = 1.25)
mtext(year, 3, 2.0, cex = 1.5)
box()

# Proportions of mature:
r <- aggregate(list(p = p), by = list(cw = round(x)), mean)
r$n <- aggregate(list(n = p), by = list(cw = round(x)), length)$n
plot(c(3.4, 5), c(0, 1), type ="n", pch = 21, cex = 0.6, xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
logit_p <- theta[["logit_p_beta0"]] + 
           theta[["logit_p_beta1"]] * log(x0) + 
           theta[["logit_p_beta2"]] * (log(x0)^2) + 
           theta[["logit_p_beta3"]] * (log(x0)^3)
p <- 1 / (1 + exp(-logit_p))  
lines(log(x0), p, lwd = 2, col = "blue")
mtext("Maturity proportion", 2, 2.5, cex = 1.25)
mtext("ln(carapace width)", 1, 2.5, cex = 1.25)
lines(log(r$cw), r$p, col = "red", lwd = 2)
box()

plot(c(3.4, 5), c(-5, 10), type ="n", pch = 21, cex = 0.6, xaxs = "i", yaxs = "i", xlab = "", ylab = "")
points(log(r$cw), log(r$p/(1-r$p)))

