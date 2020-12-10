library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

# Model parameters:
n_instar <- 7
years  <- 2006:2020
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

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_year.cpp")
dyn.load(dynlib("instar_year"))

# Define initial parameters:
parameters <- list(mu0                 = 10,                             # First instar mean size.
                   log_sigma0          = log(0.8),                       # Log-scale standard error for first instar.
                   log_hiatt_slope     = log(c(0.350, 0.0920)),          # Hiatt slope parameters.
                   log_hiatt_intercept = log(c(0.689, 8.000)),           # Hiatt intercept parameters.
                   log_growth_error    = log(c(0.01, 0.10)),             # Growth increment error inflation parameters.
                   log_mu_year         = rep(0, n_instar * n_year),      # Log-scale instar mean year interaction (n_instar x n_year).
                   log_sigma_mu_year   = -1,                             # Instar mean year interaction error term.
                   log_n_imm_year_0    = rep(4, n_instar-1),             # First year immature instar abundances (n_instar-1).
                   log_n_imm_instar_0  = rep(4, n_year),                 # First instar recruitment for all years (n_year).
                   log_sigma_n_imm_instar_0 = -1,                        # Log-scale first instar annual recruitment error parameter.
                   log_n_skp_instar_0  = rep(0, n_instar-5),             # First year skip abundances (n_instar-5).                         
                   log_n_rec_instar_0  = rep(0, n_instar-5),             # First year mature recruit abundances (n_instar-5).
                   log_n_res_instar_0  = rep(0, n_instar-5),             # First year mature residual abundances (n_instar-5).
                   selectivity_x50     = c(20, 65),                      # Size-at-50% trawl selectivity.
                   log_selectivity_slope = c(-3, -3),                    # Log-scale trawl selectivity slope.
                   logit_selectivity_proportion = -1,               
                   logit_p_skp = c(rep(-8, 4), rep(-6, n_instar-5)),     # Logit-scale skip-moulting probabilities (n_instar).
                   logit_p_mat = c(rep(-8, 4), rep(0, n_instar-6), 3),   # Logit-scale moult-to-maturity probabilities (n_instar).
                   logit_p_mat_year = rep(0, (n_instar-1) * (n_year-1)), # Logit-scale mout-to-maturity instar x year interaction (n_instar x n_year).
                   log_sigma_p_mat_year = -1,                            # Moult-to-maturity instar x year interaction error term.
                   logit_M_imm = 0,                                      # Logit-scale immature mortality.
                   logit_M_mat = 0)                                      # Logit-scale mature mortality.  


# Define random variables in model:
random <- c("log_mu_year", "log_n_imm_instar_0", "logit_p_mat_year")

# Estimate initial abundance parameters:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))
map <- update.map(map, free = c("log_n_imm_instar_0", "log_n_imm_year_0", "log_n_skp_instar_0", "log_n_rec_instar_0", "log_n_res_instar_0", "log_sigma_n_imm_instar_0"))
obj <- MakeADFun(data, parameters, DLL = "instar_year",  random = random, map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters <- update.parameters(parameters, obj)

# Add mortality parameters:
map <- update.map(map, free = c("logit_M_mat","logit_M_imm"))
obj <- MakeADFun(data, parameters, DLL = "instar_year",  random = random, map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters <- update.parameters(parameters, obj)

# Add selectivity parameters:
map <- update.map(map, free = c("selectivity_x50", "log_selectivity_slope", "logit_selectivity_proportion"))
obj <- MakeADFun(data, parameters, DLL = "instar_year",  random = random, map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters <- update.parameters(parameters, obj)

# Add moult to maturity parameters:
map <- update.map(map, free = c("logit_p_mat", "logit_p_mat_year", "log_sigma_p_mat_year"))
obj <- MakeADFun(data, parameters, DLL = "instar_year",  random = random, map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta

# Add some growth parameters:
map <- update.map(map, free = c("log_hiatt_slope", "log_hiatt_intercept"))
obj <- MakeADFun(data, parameters, DLL = "instar_year",  random = random, map = map)
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta

#map$log_growth_error    <- factor(1:length(parameters$log_growth_error))

