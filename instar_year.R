library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

# Model parameters:
n_instar <- 7
years  <- 2010:2020
sex <- 2
step <- 0.5

# Set derived quantities:
instars <- as.character(as.roman(4:(4+n_instar-1)))
n_year <- length(years)
if (sex == 1){
   n_instar <- 9
   xlim <- c(0, 140)
   ylim <- c(0, 40)
}else{
   n_instar <- 7
   xlim <- c(0, 100)
   ylim <- c(0, 40)   
}

source("instar.year.data.R")
#source("instar.plot.R")

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_year.cpp")
dyn.load(dynlib("instar_year"))

# Define initial parameters:
parameters <- list(mu0                 = 10,                        # First instar mean size.
                   log_sigma0          = log(0.8),                  # Log-scale standard error for first instar.
                   log_hiatt_slope     = log(c(0.350, 0.0920)),     # Hiatt slope parameters.
                   log_hiatt_intercept = log(c(0.689, 8.000)),      # Hiatt intercept parameters.
                   log_growth_error    = log(c(0.01, 0.10)),        # Growth increment error inflation parameters.
                   log_mu_year         = rep(0, n_instar * n_year), # Log-scale instar mean year interaction (n_instar x n_year).
                   log_sigma_mu_year   = -1,                        # Instar mean year interaction error term.
                   log_n_imm_year_0    = rep(4, n_instar),          # First year immature instar abundances (n_instar).
                   log_n_imm_instar_0  = rep(4, n_year-1),          # First instar recruitment for all years (n_year-1).
                   log_sigma_n_imm_instar_0 = -1,                   # Log-scale first instar annual recruitment error parameter.
                   log_n_skp_instar_0  = rep(0, n_instar-5),        # First year skip abundances (n_instar-5).                         
                   log_n_rec_instar_0  = rep(0, n_instar-5),        # First year mature recruit abundances (n_instar-5).
                   log_n_res_instar_0  = rep(0, n_instar-5),        # First year mature residual abundances (n_instar-5).
                   selectivity_x50     = 65,                        # Size-at-50% trawl selectivity.
                   log_selectivity_slope = -3,                      # Log-scale trawl selectivity slope.
                   logit_p_skp = c(rep(-8, 4), rep(-6, n_instar-4)), # Logit-scale skip-moulting probabilities (n_instar).
                   logit_p_mat = c(rep(-8, 4), rep(0, n_instar-6), 3, 4), # Logit-scale moult-to-maturity probabilities (n_instar).
                   logit_M_imm = 0,                                 # Logit-scale immature mortality.
                   logit_M_mat = 0)                                 # Logit-scale mature mortality.  

# Estimate initial abundance parameters:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))
logit_M_matlog_n_imm_year_0 <- factor(1:length(parameters$log_n_imm_year_0))
map$log_n_imm_instar_0 <- factor(1:length(parameters$log_n_imm_instar_0))
map$log_n_imm_year_0   <- factor(1:length(parameters$log_n_imm_year_0))
map$log_n_skp_instar_0 <- factor(1:length(parameters$log_n_skp_instar_0))
map$log_n_rec_instar_0 <- factor(1:length(parameters$log_n_rec_instar_0))
map$log_n_res_instar_0 <- factor(1:length(parameters$log_n_res_instar_0))
map$log_sigma_n_imm_instar_0 <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_mu_year", "log_n_imm_instar_0"),
                 map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")[, 1]
random <- summary(rep, "random")[, 1]
parameters$log_n_imm_instar_0 <- as.numeric(random[grep("log_n_imm_instar_0", names(random))])
parameters$log_n_imm_year_0   <- as.numeric(fixed[grep("log_n_imm_year_0", names(fixed))])
parameters$log_n_skp_instar_0 <- as.numeric(fixed[grep("log_n_skp_instar_0", names(fixed))])
parameters$log_n_rec_instar_0 <- as.numeric(fixed[grep("log_n_rec_instar_0", names(fixed))])
parameters$log_n_res_instar_0 <- as.numeric(fixed[grep("log_n_res_instar_0", names(fixed))])
parameters$log_sigma_n_imm_instar_0 <- as.numeric(theta["log_sigma_n_imm_instar_0"])

# Add mortality parameters:
map$logit_M_mat <- factor(1)
map$logit_M_imm <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_mu_year", "log_n_imm_instar_0"),
                 map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")[, 1]
random <- summary(rep, "random")[, 1]
parameters$log_n_imm_instar_0 <- as.numeric(random[grep("log_n_imm_instar_0", names(random))])
parameters$log_n_imm_year_0   <- as.numeric(fixed[grep("log_n_imm_year_0", names(fixed))])
parameters$log_n_skp_instar_0 <- as.numeric(fixed[grep("log_n_skp_instar_0", names(fixed))])
parameters$log_n_rec_instar_0 <- as.numeric(fixed[grep("log_n_rec_instar_0", names(fixed))])
parameters$log_n_res_instar_0 <- as.numeric(fixed[grep("log_n_res_instar_0", names(fixed))])
parameters$log_sigma_n_imm_instar_0 <- as.numeric(theta["log_sigma_n_imm_instar_0"])
parameters$logit_M_mat <- as.numeric(theta["logit_M_mat"])
parameters$logit_M_imm <- as.numeric(theta["logit_M_imm"])

# Add selectivity parameters:
map$selectivity_x50 <- factor(1)
map$log_selectivity_slope <-  factor(1)                
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_mu_year", "log_n_imm_instar_0"),
                 map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")[, 1]
random <- summary(rep, "random")[, 1]
parameters$log_n_imm_instar_0 <- as.numeric(random[grep("log_n_imm_instar_0", names(random))])
parameters$log_n_imm_year_0   <- as.numeric(fixed[grep("log_n_imm_year_0", names(fixed))])
parameters$log_n_skp_instar_0 <- as.numeric(fixed[grep("log_n_skp_instar_0", names(fixed))])
parameters$log_n_rec_instar_0 <- as.numeric(fixed[grep("log_n_rec_instar_0", names(fixed))])
parameters$log_n_res_instar_0 <- as.numeric(fixed[grep("log_n_res_instar_0", names(fixed))])
parameters$log_sigma_n_imm_instar_0 <- as.numeric(theta["log_sigma_n_imm_instar_0"])
parameters$logit_M_mat <- as.numeric(theta["logit_M_mat"])
parameters$logit_M_imm <- as.numeric(theta["logit_M_imm"])
parameters$selectivity_x50 <- as.numeric(theta["selectivity_x50"])
parameters$log_selectivity_slope <- as.numeric(theta["log_selectivity_slope"])
                    
# Add growth parameters:
map$log_hiatt_slope     <- factor(1:length(parameters$log_hiatt_slope))
map$log_hiatt_intercept <- factor(1:length(parameters$log_hiatt_intercept))       
map$log_growth_error    <- factor(1:length(parameters$log_growth_error))
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_mu_year", "log_n_imm_instar_0"),
                 map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")[, 1]
random <- summary(rep, "random")[, 1]
parameters$log_n_imm_instar_0 <- as.numeric(random[grep("log_n_imm_instar_0", names(random))])
parameters$log_n_imm_year_0   <- as.numeric(fixed[grep("log_n_imm_year_0", names(fixed))])
parameters$log_n_skp_instar_0 <- as.numeric(fixed[grep("log_n_skp_instar_0", names(fixed))])
parameters$log_n_rec_instar_0 <- as.numeric(fixed[grep("log_n_rec_instar_0", names(fixed))])
parameters$log_n_res_instar_0 <- as.numeric(fixed[grep("log_n_res_instar_0", names(fixed))])
parameters$log_sigma_n_imm_instar_0 <- as.numeric(theta["log_sigma_n_imm_instar_0"])
parameters$logit_M_mat <- as.numeric(theta["logit_M_mat"])
parameters$logit_M_imm <- as.numeric(theta["logit_M_imm"])
parameters$selectivity_x50 <- as.numeric(theta["selectivity_x50"])
parameters$log_selectivity_slope <- as.numeric(theta["log_selectivity_slope"])
parameters$log_hiatt_slope     <- as.numeric(fixed[grep("log_hiatt_slope", names(fixed))])
parameters$log_hiatt_intercept <- as.numeric(fixed[grep("log_hiatt_intercept", names(fixed))])
parameters$log_growth_error    <- as.numeric(fixed[grep("log_growth_error", names(fixed))])

# Add moult to maturity parameters:
map$logit_p_mat <- factor(1:length(parameters$logit_p_mat))
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_mu_year", "log_n_imm_instar_0"),
                 map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")[, 1]
random <- summary(rep, "random")[, 1]
parameters$log_n_imm_instar_0 <- as.numeric(random[grep("log_n_imm_instar_0", names(random))])
parameters$log_n_imm_year_0   <- as.numeric(fixed[grep("log_n_imm_year_0", names(fixed))])
parameters$log_n_skp_instar_0 <- as.numeric(fixed[grep("log_n_skp_instar_0", names(fixed))])
parameters$log_n_rec_instar_0 <- as.numeric(fixed[grep("log_n_rec_instar_0", names(fixed))])
parameters$log_n_res_instar_0 <- as.numeric(fixed[grep("log_n_res_instar_0", names(fixed))])
parameters$log_sigma_n_imm_instar_0 <- as.numeric(theta["log_sigma_n_imm_instar_0"])
parameters$logit_M_mat <- as.numeric(theta["logit_M_mat"])
parameters$logit_M_imm <- as.numeric(theta["logit_M_imm"])
parameters$selectivity_x50 <- as.numeric(theta["selectivity_x50"])
parameters$log_selectivity_slope <- as.numeric(theta["log_selectivity_slope"])
parameters$log_hiatt_slope     <- as.numeric(fixed[grep("log_hiatt_slope", names(fixed))])
parameters$log_hiatt_intercept <- as.numeric(fixed[grep("log_hiatt_intercept", names(fixed))])
parameters$log_growth_error    <- as.numeric(fixed[grep("log_growth_error", names(fixed))])
parameters$logit_p_mat         <- as.numeric(fixed[grep("logit_p_mat", names(fixed))])

# Add remaining growth parameters (stop)
map$mu0 <- factor(1)
map$log_sigma0 <- factor(1)
map$log_sigma_mu_year <- factor(1)
map$log_mu_year <- factor(1:length(parameters$log_mu_year))
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_mu_year", "log_n_imm_instar_0"),
                 map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")[, 1]
random <- summary(rep, "random")[, 1]
parameters$log_n_imm_instar_0 <- as.numeric(random[grep("log_n_imm_instar_0", names(random))])
parameters$log_n_imm_year_0   <- as.numeric(fixed[grep("log_n_imm_year_0", names(fixed))])
parameters$log_n_skp_instar_0 <- as.numeric(fixed[grep("log_n_skp_instar_0", names(fixed))])
parameters$log_n_rec_instar_0 <- as.numeric(fixed[grep("log_n_rec_instar_0", names(fixed))])
parameters$log_n_res_instar_0 <- as.numeric(fixed[grep("log_n_res_instar_0", names(fixed))])
parameters$log_sigma_n_imm_instar_0 <- as.numeric(theta["log_sigma_n_imm_instar_0"])
parameters$logit_M_mat <- as.numeric(theta["logit_M_mat"])
parameters$logit_M_imm <- as.numeric(theta["logit_M_imm"])
parameters$selectivity_x50 <- as.numeric(theta["selectivity_x50"])
parameters$log_selectivity_slope <- as.numeric(theta["log_selectivity_slope"])
parameters$log_hiatt_slope     <- as.numeric(fixed[grep("log_hiatt_slope", names(fixed))])
parameters$log_hiatt_intercept <- as.numeric(fixed[grep("log_hiatt_intercept", names(fixed))])
parameters$log_growth_error    <- as.numeric(fixed[grep("log_growth_error", names(fixed))])
parameters$logit_p_mat         <- as.numeric(fixed[grep("logit_p_mat", names(fixed))])

parameters$mu0  <- as.numeric(theta["mu0"])
parameters$log_sigma0  <- as.numeric(theta["log_sigma0"])
parameters$log_sigma_mu_year <- as.numeric(theta["log_sigma_mu_year"])
parameters$log_mu_year <- as.numeric(random[grep("log_mu_year", names(random))])

