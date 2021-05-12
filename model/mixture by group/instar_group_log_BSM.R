rm(list = ls())
library(TMB)
library(gulf.data)
library(gulf.graphics)

source("/Users/crustacean/Desktop/gulf-snow-crab-EAFM/R/TMB utilities.R")
source("/Users/crustacean/Desktop/gulf-snow-crab-EAFM/model/mixture by group/plot.instar.group.R")

sex    <- 2          # Crab sex.
maturity <- 0        # Crab maturity.

mu_instars <- c(10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)
names(mu_instars) <- 4:11
roman <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV")
   
setwd("model/mixture by group")

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

# Load BSM data:
b <- read.csv("BSM/BSM_fem_94-98.csv")

# Define data:
if (sex == 2){
   if (maturity == 0) mu_instars <- mu_instars[as.character(4:10)]
   if (maturity == 1) mu_instars <- mu_instars[as.character(9:11)]
} 
b <- b[which(b$maturity == "0"), ]

b$date <- as.character(date(b))
b <- sort(b, by = c("date", "tow.number", "number"))

# Create tow table and tow index:
tows <- aggregate(list(n = b$carapace.width),  b[c("date", "tow.number")], length)
tows <- sort(tows, by = c("date", "tow.number"))
b$tow <- match(b[c("date", "tow.number")], tows[c("date", "tow.number")]) - 1
tows <- aggregate(list(n = b$carapace.width),  b[c("date", "tow", "tow.number")], length)
tows <- sort(tows, by = c("date", "tow.number"))


if (maturity == 0) xlim = c(2.5, 4.25)
if (maturity == 1) xlim = c(3.5, 4.5)


mu_instars <- c(3.22, 4.63, 6.62, 10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)
names(mu_instars) <- 1:11
mu_instars <- mu_instars[as.character(1:9)]

b$carapace.width <- as.numeric(b$carapace.width)
b <- b[!is.na(b$carapace.width), ]
          
# Set up data:
#r <- aggregate(list(f = b$year), list(x = round(log(b$carapace.width), 2), year = b$year), length)
r <- aggregate(list(f = b$year), 
               by = list(date = b$date,
                         tow.number = b$tow.number,
                         group = b$tow, 
                         year = b$year, 
                         x = round(log(b$carapace.width), 2)), 
               length)

data <- list(x = r$x, 
             f = r$f)
data$group <- r$group
data$year <- r$year
data$n_instar <- length(mu_instars)
data$n_group  <- max(data$group) + 1
data$precision <- rep(0.1, length(data$x))
data$tow.number <- r$tow.number
data$date <- r$date

# Define initial parameters:
parameters <- list(mu_instar_0 = as.numeric(log(mu_instars[1])),   # Mean size of the first instar.
                   log_increment = log(0.405),                         # Log-scale mean parameter associated with instar growth increments.
                   log_increment_delta = log(0.013),
                   log_sigma_mu_instar_group = -4,                 # Log-scale error for instar-group means random effect.
                   mu_instar_group = rep(0, length(mu_instars) * length(unique(data$group))),  # Instar-group means random effect.
                   log_sigma = -2,                                 # Log-scale instar standard error.
                   log_sigma_instar_group = rep(0, length(mu_instars) * length(unique(data$group))),
                   log_sigma_sigma_instar_group = -4,
                   mu_logit_p = 4,                                 # Instar proportions log-scale mean parameter.
                   log_sigma_logit_p_instar_group = -1,            # Instar-group proportions log-scale error parameter.
                   logit_p_instar_group = rep(0, (length(mu_instars)-1) * length(unique(data$group)))) # Multi-logit-scale parameters for instar-group proportions by year.
  
parameters <- parameters[parameters.cpp("instar_group_log.cpp")]

random <- c("mu_instar_group", "logit_p_instar_group", "log_sigma_instar_group")

# Fit mixture proportions:
map <- lapply(parameters, function(x) as.factor(rep(NA, length(x))))
map$mu_logit_p <- factor(1)
map$log_sigma_logit_p_instar_group <- factor(1)
map$logit_p_instar_group <- factor(1:length(parameters$logit_p_instar_group))
obj <- MakeADFun(data[data.cpp("instar_group_log.cpp")], 
                 parameters[parameters.cpp("instar_group_log.cpp")], 
                 random = random, map = map, DLL = "instar_group_log")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))

map$log_sigma <- factor(1)

map$log_increment <- factor(1)
map$log_increment_delta <- factor(1)

map$mu_instar_0 <- factor(1)

map$mu_instar_group <- factor(1:length(parameters$mu_instar_group))

map$log_sigma_mu_instar_group <- factor(1)




