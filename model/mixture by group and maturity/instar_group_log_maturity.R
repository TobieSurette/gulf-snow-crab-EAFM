rm(list = ls())
library(TMB)
library(gulf.data)
library(gulf.graphics)

source("/Users/crustacean/Desktop/gulf-snow-crab-EAFM/R/TMB utilities.R")
source("/Users/crustacean/Desktop/gulf-snow-crab-EAFM/model/mixture by group/plot.instar.group.R")

compile("instar_group_log_maturity.cpp")
dyn.load(dynlib("instar_group_log_maturity"))

BSM <- FALSE
sex      <- 2           # Crab sex.
maturity <- 0           # Crab maturity.
years    <- 2017:2020   # Define survey years.

# Define instar mean sizes:
mu_instars <- c(3.22, 4.63, 6.62, 10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)
names(mu_instars) <- 1:11
mu_instars <- mu_instars[as.character(4:10)]
   
setwd("/Users/crustacean/Desktop/gulf-snow-crab-EAFM/model/mixture by group and maturity")

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

if (BSM){
   # Load BSM data:
   b <- read.csv("BSM/BSM_fem_94-98.csv")
   b <- b[which(b$maturity %in% c("0", "1")), ]
   b$maturity <- as.numeric(b$maturity)
   b$date <- as.character(date(b))
   b <- sort(b, by = c("date", "tow.number", "number"))
}else{
   # Define data:
   b <- read.scsbio(years, survey = "regular")
   b <- b[which(b$sex == sex), ]
   b <- b[which(b$carapace.width > 6), ]
   b$maturity <- is.mature(b)
   b <- b[-which(b$maturity & !is.new.shell(b)), ]
   b$year <- year(b)
   b <- sort(b, by = c("date", "tow.number", "crab.number"))
   b$carapace.width <- as.numeric(b$carapace.width)
   b <- b[!is.na(b$carapace.width), ]
}

# Create tow table and tow index:
tows  <- aggregate(list(n = b$carapace.width),  b[c("date", "tow.number")], length)
tows  <- sort(tows, by = c("date", "tow.number"))
b$tow <- match(b[c("date", "tow.number")], tows[c("date", "tow.number")]) - 1
tows  <- aggregate(list(n = b$carapace.width),  b[c("date", "tow", "tow.number")], length)
tows$n_imm <- aggregate(list(x = !b$maturity),  b[c("date", "tow.number")], sum, na.rm = TRUE)$x
tows$n_mat <- aggregate(list(x = b$maturity),  b[c("date", "tow.number")], sum, na.rm = TRUE)$x
tows  <- sort(tows, by = c("date", "tow.number"))
s <- read.scsset(year = years, valid = 1, survey = "regular")
ix <- match(tows[c("date", "tow.number")], s[c("date", "tow.number")])
tows$swept.area <- s$swept.area[ix]
   
# Set up data:
#r <- aggregate(list(f = b$year), list(x = round(log(b$carapace.width), 2), year = b$year), length)
data <- aggregate(list(f = b$year), 
                  by = list(date = b$date,
                            tow.number = b$tow.number,
                            group = b$tow, 
                            year = b$year, 
                            maturity = b$maturity,
                            x = round(log(b$carapace.width), 2)), 
                   length)
data <- as.list(data)
data$n_instar <- length(mu_instars)
data$n_group  <- max(data$group) + 1
data$precision <- rep(0.1, length(data$x))

# Define initial parameters:
parameters <- list(mu_instar_0 = as.numeric(log(mu_instars[1])),   # Mean size of the first instar.
                   log_increment = log(0.345),                         # Log-scale mean parameter associated with instar growth increments.
                   log_increment_delta = log(0.010),
                   log_sigma_mu_instar_group = -4,                 # Log-scale error for instar-group means random effect.
                   mu_instar_group = rep(0, length(mu_instars) * length(unique(data$group))),  # Instar-group means random effect.
                   mu_instar_group_mature = rep(0, length(mu_instars) * length(unique(data$group))),  # Instar-group means random effect.
                   log_sigma = -2,             # Log-scale instar standard error.
                   log_sigma_instar_group = rep(0, length(mu_instars) * length(unique(data$group))),
                   log_sigma_sigma_instar_group = -4,
                   mu_logit_p = rep(4, data$n_instar-1),           # Instar proportions log-scale mean parameter.
                   mu_logit_p_mature = rep(4, data$n_instar-1), 
                   log_sigma_logit_p_instar_group = -1,            # Instar-group proportions log-scale error parameter.
                   logit_p_instar_group = rep(0, (length(mu_instars)-1) * length(unique(data$group))), # Multi-logit-scale parameters for instar-group proportions by year.
                   logit_p_instar_group_mature = rep(0, (length(mu_instars)-1) * length(unique(data$group)))
                   ) 
  
random     <- c("mu_instar_group", "mu_instar_group_mature", "logit_p_instar_group", "logit_p_instar_group_mature", "log_sigma_instar_group")

# Fit mixture proportions for immatures:
map <- lapply(parameters, function(x) as.factor(rep(NA, length(x))))
map$mu_logit_p <- factor(1:(data$n_instar-1))
map$log_sigma_logit_p_instar_group <- factor(1)
map$logit_p_instar_group <- factor(1:length(parameters$logit_p_instar_group))
map$mu_logit_p_mature <- factor(1:(data$n_instar-1))
map$log_sigma_logit_p_instar_group <- factor(1)
map$logit_p_instar_group_mature <- factor(1:length(parameters$logit_p_instar_group))
map$log_sigma <- factor(1)

obj <- MakeADFun(data[data.cpp("instar_group_log_maturity.cpp")], 
                 parameters[parameters.cpp("instar_group_log_maturity.cpp")], 
                 random = random, map = map, DLL = "instar_group_log_maturity")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))

          
#map$log_sigma <- factor(rep(1,data$n_instar))


map$log_increment       <- factor(1)
map$log_increment_delta <- factor(1)

map$mu_instar_0 <- factor(1)
map$mu_instar_0 <- factor(NA)

map$mu_instar_group <- factor(1:length(parameters$mu_instar_group))
map$mu_instar_group_mature <- factor(1:length(parameters$mu_instar_group_mature))

map$log_sigma_mu_instar_group <- factor(1)




