library(gulf)
library(TMB)

setwd("U:/Snow Crab/Morphometry")
source("U:/TMB/TMB utilities.R")
source("C:/gulf package/gulf/R/observer.scobs.R")

# Compile and load TMB program:
cpp.file <- "Morphometry_observer.cpp" 
compile(cpp.file)
dll.file <- gsub(".cpp", "", cpp.file)  
dyn.load(dynlib(dll.file))

year <- 2017

# Load data:
x <- read.scobs(year = year)

# Format data:
x$observer <- observer(x)

x$week <- week(date(x))
index <- which(!is.na(x$fishing.grid) & (x$week >= 15) & 
               !(x$fishing.grid %in% c("GZ32", "HH50")) & !is.na(x$chela.height.right) & !is.na(x$carapace.width) &
               (x$carapace.width >= 50) & (x$carapace.width <= 150))
x <- x[index, ]

observers <- c("DOMINIC VIGNEUX", "PATRICK TREMBLAY", "GENEVIEVE MARTEL", "DEREK TRUPPNER", "CHRISTIAN THERIAULT", "ALLAN ARSENAULT", 
               "EDITH BERGERON", "YVES LAROCQUE", "STEPHANE CHIASSON", "JEAN-PHILIPPE BERTIN", "SYLVIE CHIASSON", "ANDRE PARADIS")  

#x <- x[x$observer %in% observers, ]

# Define input data:
observers <- sort(unique(x$observer))
data <- list(y = round(x$chela.height.right),
             cw = round(x$carapace.width),
             observer = match(x$observer, observers),
             grid = match(x$fishing.grid, sort(unique(x$fishing.grid))),
             week = as.integer(x$week - min(x$week) + 1), 
             trip = match(x$trip.number, sort(unique(x$trip.number))),
             observers = observers,
             grids = sort(unique(x$fishing.grid)),
             weeks = sort(unique(x$week)),
             trips = sort(unique(x$trip.number)))

# Calculate grid coordinates in kilometers:
data$grid_x = apply(grid.corners(data$grids)[, c(1, 3)], 1, mean)
data$grid_y = apply(grid.corners(data$grids)[, c(1, 3)+1], 1, mean)
tmp <- deg2km(data$grid_x, data$grid_y)
data$grid_x <- tmp[, 1]
data$grid_y <- tmp[, 2]
            
# Define initial parameters:  
pars <- list(log_alpha_immature = -2.807869,                                 
             log_beta_immature = 1.250458,
             log_alpha_mature = -3.047809,
             log_beta_mature = 1.353115,
             mu_observer = 0,
             observer_effect = rep(0, max(data$observer)),
             log_sigma_observer = -1,            
             trip_effect = rep(0, max(data$trip)),
             log_sigma_trip = -1,     
             log_sigma = -3.138193,
             observer_error_effect = rep(0, max(data$observer)),
             log_sigma_observer_error = -1,              
             log_sigma_outlier = 2,
             alpha_outlier = -4,
             observer_outlier_effect = rep(0, max(data$observer)),
             log_sigma_observer_outlier = -1,
             alpha_p = -5.02779864,
             beta_p =  0.07004936,          
             observer_effect_p = rep(0, max(data$observer)),  
             observer_effect_beta_p = rep(0, max(data$observer)),  
             log_sigma_observer_p = -1,
             log_sigma_observer_beta_p = -1,
             grid_effect_p = rep(0, max(data$grid)),             
             log_sigma_grid_p = -1,
             week_effect_p = rep(0, max(data$week)),             
             log_sigma_week_p = -1)
             
# Function to update parameters:
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
                                   
subset.data <- function(x, index){
   flag <- is.data.frame(x)
   x <- as.list(x)
   vars <- names(x)[unlist(lapply(x, length)) == max(unlist(lapply(x, length)))]
   for (i in 1:length(vars)) x[[vars[i]]] <- x[[vars[i]]][index]
   if (flag) x <- as.data.frame(x)
   return(x)
}
    
random.pars <- c("observer_effect", "observer_error_effect", "observer_outlier_effect", "observer_effect_p", "observer_effect_beta_p", "grid_effect_p", "week_effect_p", "trip_effect")
                                                                                  
index <- 1:nrow(x) # sort(sample(1:nrow(data), 10000))  #  #     
       
active.parameters <- c("mu_observer", "observer_effect", "log_sigma_observer") 
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index), 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars[parameters.cpp(cpp.file)], active.parameters),
                 random = "observer_effect",
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 200, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)
                           
# Maturity proportions:
active.parameters <- unique(c(active.parameters, "alpha_p", "beta_p"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index), 
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 300, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

# Overall outlier proportions:
active.parameters <- unique(c(active.parameters, "alpha_outlier"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index),
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 300, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

# Maturity observer effect proportions:
active.parameters <- unique(c(active.parameters, "observer_effect_p", "log_sigma_observer_p"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index),
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 300, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

# Maturity observer effect proportions slope:
active.parameters <- unique(c(active.parameters, "observer_effect_beta_p", "log_sigma_observer_beta_p"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index),
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 300, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

# Observer error terms:
active.parameters <- unique(c(active.parameters, "log_sigma", "observer_error_effect", "log_sigma_observer_error"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index),
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 300, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

# Outlier terms:
active.parameters <- unique(c(active.parameters, "log_sigma_outlier", "observer_outlier_effect", "log_sigma_observer_outlier"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index),
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars, active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 500, trace = 3))
#theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

# Grid and week effects:
active.parameters <- unique(c(active.parameters, "grid_effect_p", "log_sigma_grid_p", "week_effect_p", "log_sigma_week_p"))
obj <- MakeADFun(data = subset.data(data[data.cpp(cpp.file)], index),
                 parameters = pars[parameters.cpp(cpp.file)],
                 map = control.parameters(pars[parameters.cpp(cpp.file)], active.parameters),
                 random = random.pars,
                 DLL = dll.file)
                 
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
pars <- update.parameters(pars, fixed, random)

clg()
for (i in 1:24) chela.plot(data, pars, observer = i)


setdiff(names(pars), active.parameters)


p <- pars
    # Plot observer effects:
    vars <- names(p)[intersect(grep("observer", names(p)), grep("_effect", names(p)))]
    m <- kronecker(c(1:length(vars), 0), matrix(rep(1, 25), nrow = 5))
    m <- rbind(0, cbind(0, m, 0))
    windows(height = 10, width = 7)
    layout(m)
    par(mar =c(0, 0, 0, 0))
    for (i in 1:length(vars)){
       dbarplot(p[[vars[i]]], xaxt = "n", xlab = "", ylab = "") 
       grid()
       mtext(vars[i], 2, 2.5)
    }
    axis(1, at = 1:length(data$observers), labels = data$observers, las = 2)

# Chela bias plot:
dbarplot(pars$mu_observer + pars$observer_effect)
text(1:length(data$observers), 0.5 * (pars$mu_observer + pars$observer_effect), data$observers, srt = 90)

# Chela error versus bias plot
sigma <- exp(pars$log_sigma + pars$observer_error_effect)
plot(pars$mu_observer + pars$observer_effect, 3.2 * sqrt(exp(sigma^2)-1), pch = 21, cex = 1.40, bg = "red",
     xlab = "Chela mesurement bias (mm)", ylab = "Variation (mm)")
     #,   xlim = c(-5, 2), ylim = c(0.15, 0.2))
text(pars$mu_observer + pars$observer_effect, 3.2 * sqrt(exp(sigma^2)-1), data$observers, pos = 2, cex = 0.65)    
lines(c(0, 0), par("usr")[3:4], lwd = 2)    

# Observer effect:
r <- random[(rownames(random) == "observer_effect"), ]
dbarplot(r[, 1], ylim = c(min(r[,1] - 2*r[,2]), max(r[,1] + 2*r[,2])), ylab = "Effect value")
for (i in 1:nrow(r)){
   lines(c(i, i), c(r[i,1] - 1.96 * r[i,2], r[i,1] + 1.96 * r[i,2]), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] - 1.96 * r[i,2], 2), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] + 1.96 * r[i,2], 2), lwd = 2, col = "black")
}

r <- random[(rownames(random) == "observer_effect_p"), ]
dbarplot(r[, 1], ylim = c(min(r[,1] - 2*r[,2]), max(r[,1] + 2*r[,2])), ylab = "Effect value")
for (i in 1:nrow(r)){
   lines(c(i, i), c(r[i,1] - 1.96 * r[i,2], r[i,1] + 1.96 * r[i,2]), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] - 1.96 * r[i,2], 2), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] + 1.96 * r[i,2], 2), lwd = 2, col = "black")
}

r <- random[(rownames(random) == "observer_effect_beta_p"), ]
dbarplot(r[, 1], ylim = c(min(r[,1] - 2*r[,2]), max(r[,1] + 2*r[,2])), ylab = "Effect value")
for (i in 1:nrow(r)){
   lines(c(i, i), c(r[i,1] - 1.96 * r[i,2], r[i,1] + 1.96 * r[i,2]), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] - 1.96 * r[i,2], 2), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] + 1.96 * r[i,2], 2), lwd = 2, col = "black")
}


# Grid effect:
r <- random[(rownames(random) == "grid_effect_p"), ]
dbarplot(r[, 1], ylim = c(min(r[,1] - 2*r[,2]), max(r[,1] + 2*r[,2])), ylab = "Effect value")
for (i in 1:nrow(r)){
   lines(c(i, i), c(r[i,1] - 1.96 * r[i,2], r[i,1] + 1.96 * r[i,2]), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] - 1.96 * r[i,2], 2), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] + 1.96 * r[i,2], 2), lwd = 2, col = "black")
}

# Week effect:
r <- random[(rownames(random) == "week_effect_p"), ]
dbarplot(r[, 1], ylim = c(min(r[,1] - 2*r[,2]), max(r[,1] + 2*r[,2])), ylab = "Effect value")
for (i in 1:nrow(r)){
   lines(c(i, i), c(r[i,1] - 1.96 * r[i,2], r[i,1] + 1.96 * r[i,2]), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] - 1.96 * r[i,2], 2), lwd = 2, col = "black")
   lines(c(i-0.15, i+0.15), rep(r[i,1] + 1.96 * r[i,2], 2), lwd = 2, col = "black")
}


for (i in 1:length(data$observers)) chela.plot(data, pars, observer = i)

