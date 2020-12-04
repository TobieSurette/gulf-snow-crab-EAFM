library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

years  <- 2017
sex <- 2
step <- 0.5
if (sex == 1){
   n_instar <- 9
   xlim <- c(0, 140)
   ylim <- c(0, 40)
}else{
   n_instar <- 7
   xlim <- c(0, 100)
   ylim <- c(0, 40)   
}

source("instar.data.R")
source("instar.plot.R")

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar.cpp")
dyn.load(dynlib("instar"))

# Define initial parameters:
parameters <- list(mu0                     = 10,                       # First instar mean size.
                   log_sigma0              = log(0.8),                 # Log-scale standard error for first instar.
                   log_hiatt_slope         = log(c(0.350, 0.0920)),    # Hiatt slope parameters.
                   log_hiatt_intercept     = log(c(0.689, 8.000)),     # Hiatt intercept parameters.
                   log_growth_error        = log(c(0.01, 0.10)),       # Growth increment error inflation parameters
                   log_lambda_alpha        = 7.5,                      # Log-scale global mean density.
                   log_lambda_instar       = rep(0, n_instar),         # Log-scale instar effect.
                   log_sigma_lambda_instar = -1,                       # Log-scale instar effect error parameter.
                   logit_scale             = rep(0, n_instar))

# Fit instar abundances:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))
map$log_lambda_alpha  <- factor(1)
map$log_lambda_instar <- factor(1:length(parameters$log_lambda_instar))
map$logit_scale       <- factor(1:length(parameters$logit_scale))
map$log_sigma_lambda_instar <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_count", map = map, random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$logit_scale <- obj$report()$logit_scale

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_count", map = map, random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$logit_scale <- obj$report()$logit_scale
parameters$log_growth_error <- theta[grep("log_growth_error", names(theta))]
parameters$log_sigma0 <- theta[grep("log_sigma0", names(theta))]

obj <- MakeADFun(data, parameters, DLL = "instar_count", random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
obj$par <- theta

plot.instar(obj, data, xlim = xlim, ylim = c(0, 180))

100 * obj$report()$scale
