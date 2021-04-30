library(gulf)
library(TMB)

year <- 2017

setwd("U:/Snow Crab/Morphometry")
source("U:/TMB/TMB utilities.R")
source("splm.R")
source("chela.plot.splm.R")

# Compile and load 
cpp.file <- "Morphometry_splm.cpp" 
compile(cpp.file)
dll.file <- gsub(".cpp", "", cpp.file)  
dyn.load(dynlib(dll.file))

update.parameters <- function(x, fixed, random){
   # Update parameter list:
   if (!missing(fixed)){
      if (is.null(dim(fixed))){
          str <- names(fixed)
          fixed <- matrix(fixed) 
          rownames(fixed) <- str
      }
      for (i in 1:length(x)){
         if (names(x[i]) %in% rownames(fixed)){
            x[[i]] <- as.numeric(fixed[names(x[i]) == rownames(fixed),1])
         }
      }
   }
    
   if (!missing(random)){
      for (i in 1:length(x)){
         if (names(x[i]) %in% rownames(random)){
            x[[i]] <- as.numeric(random[names(x[i]) == rownames(random),1])
         }
      }
   }
   
   return(x)   
}

# Load data:
s <- read.scset(year = year)
x <- read.scbio(year = year)
x <- x[x$sex == 1, ]
x$station <- as.numeric(substr(x$tow.id, 3, 5)) # match(x$tow.id, sort(unique(x$tow.id))) #
index <- match(x$tow.id, s$tow.id)
x$longitude <- longitude(s)[index]
x$latitude <- latitude(s)[index]

# Define data:
index <- !is.na(x$chela.height) & !is.na(x$carapace.width)
data <- list(y         = x$chela.height[index],                # Chela height.
             cw        = x$carapace.width[index],             # Carapace width.
             station   = x$station[index],
             longitude = x$longitude[index],
             latitude  = x$latitude[index],
             sampler   = x$sampler[index],
             cw_imm    = x$carapace.width[which((x$carapace.width <= 40) & is.na(x$chela.height))]) # Morphometric immature.

# Morphometric plot:
plot(log(x$carapace.width), log(x$chela.height), pch = 21, bg = "grey", cex = 0.25, xlim = c(3.5, 5.0), ylim = c(1.5, 3.5))

# Initialize parameters
parameters = list(log_alpha_immature = -2.807869,                                
                  log_alpha_mature = -3.047809,
                  log_beta_immature = c(1.24, 1.3),
                  log_beta_mature = 1.353115,
                  log_precision_immature = -4.6,
                  transition_immature = 4,
                  log_sigma = -3.138193,                  
                  log_sigma_outlier = 2,         # Log-scale addiotnal variability of data outliers.         
                  alpha_outlier = -4,
                  eta_alpha = -7.406967,         # Logit-linear intercept for mixture proportions.
                  eta_beta = c(0.10268507, 0.02624539, 0.25488259), # Logit-linear slopes for mixture proportions.
                  eta_transition = c(60, 100),   # Logit-linear transition sizes for mixture proportions.
                  log_eta_precision = 1,         # Logit-linear transition precision for mixture proportions.
                  station_effect_immature = rep(0, 355),
                  log_sigma_station_immature  = -1,
                  station_effect_mature = rep(0, 355),
                  log_sigma_station_mature  = -1)
                  
names(parameters) %in% parameters.cpp(cpp.file)
parameters.cpp(cpp.file) %in% names(parameters)  


pars <- parameters

chela.plot.splm(data, pars)

random.parameters <- c("station_effect_immature", "station_effect_mature")

active.parameters <- c("eta_alpha", "eta_beta") 
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
pars <- update.parameters(pars, fixed)
chela.plot.splm(data, pars)

active.parameters <- unique(c(active.parameters, "log_alpha_immature", "log_alpha_mature", "log_beta_immature", "log_beta_mature", "log_sigma"))
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
pars <- update.parameters(pars, fixed)
chela.plot.splm(data, pars)

# Outlier parameters:
active.parameters <- unique(c(active.parameters, "alpha_outlier"))
active.parameters <- active.parameters[active.parameters %in% names(pars)]

active.parameters %in% parameters.cpp(cpp.file)
parameters.cpp(cpp.file)[!(parameters.cpp(cpp.file) %in% active.parameters)]
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)
chela.plot.splm(data, pars)

active.parameters <- unique(c(active.parameters, c("eta_transition", "transition_immature")))
active.parameters <- active.parameters[active.parameters %in% names(pars)]
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)
chela.plot.splm(data, pars)
 
active.parameters <- unique(c(active.parameters, c("log_eta_precision", "log_precision_immature")))
active.parameters <- active.parameters[active.parameters %in% names(pars)]
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)
chela.plot.splm(data, pars)

active.parameters <- c("station_effect_immature", "log_sigma_station_immature", "station_effect_mature", "log_sigma_station_mature")
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 500, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)
chela.plot.splm(data, pars)

pars$log_sigma_station_immature <- -4

active.parameters <- setdiff(parameters.cpp(cpp.file), "log_sigma_outlier")
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.parameters,
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)
chela.plot.splm(data, pars)





