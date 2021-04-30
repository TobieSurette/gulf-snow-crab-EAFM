library(gulf)
library(TMB)

setwd("U:/Snow Crab/Morphometry")
source("U:/TMB/TMB utilities.R")

# Function to update parameters:
update.parameters <- function(x, fixed, random){
   if (!missing(fixed)){
      for (i in 1:length(x)){
         if (names(x[i]) %in% rownames(fixed)){
            x[[i]] <- as.numeric(fixed[names(x[i]) == rownames(fixed),1])
         }
      }
   }
    
   if (!missing(random)){
      for (i in 1:length(x)){
         if (names(x[i]) %in% rownames(random)){
            x[[i]] <- as.numeric(random[grep(names(x[i]), rownames(random)), 1])
         }
      }
   }
   
   return(x)   
}

x <- read.scbio(year = 2017)
x <- x[x$sex == 1, ]
x$maturity <- is.mature(x, probability = TRUE) >= 0.5
x$station <- as.numeric(substr(x$tow.id, 3, 5)) # match(x$tow.id, sort(unique(x$tow.id))) #

# Proportion of matures:
res <- aggregate(list(n = x$maturity), by = list(cw = round(x$carapace.width)), length)
res$k <- aggregate(list(k = x$maturity), by = list(cw = round(x$carapace.width)), sum)$k
res$p <- res$k / res$n
windows();
dbarplot(res$p, res$cw, xlim = c(40, 145))

# Define data:
index <- !is.na(x$chela.height) & !is.na(x$carapace.width)
data <- data.frame(y = x$chela.height[index],                # Chela height.
                   cw = x$carapace.width[index],             # Carapace width.
                   maturity = x$maturity[index],             # Morphometric maturity.
                   station = x$station[index])             

imm <- lm(log(y) ~ log(cw), data = data[!data$maturity, ])
mat <- lm(log(y) ~ log(cw), data = data[data$maturity, ])

# Morphometric plot:
plot(log(x$carapace.width), log(x$chela.height), pch = 21, bg = "grey", cex = 0.25, xlim = c(3.5, 5.0), ylim = c(1.5, 3.5))
abline(imm, col = "red", lwd = 2)
abline(mat, col = "blue", lwd = 2)

# Initialize parameters
parameters = list(log_alpha_immature = -2.879,                                 
                  log_alpha_mature = -3.134,
                  log_beta_immature = 1.265,
                  log_beta_mature = 1.371,
                  log_sigma_immature = -2.3,
                  log_sigma_mature = -2.3,
                  eta = c(0, 0),
                  logit_p_kurtotic = -2,
                  log_sigma_kurtotic = -2.3, 
                  log_sigma_station_immature = -2,
                  log_sigma_station_mature = -2,
                  station_immature_effect = rep(0, 355),
                  station_mature_effect = rep(0, 355))
                  
 
# Compile and load 
cpp.file <- "Morphometry.cpp" 
compile(cpp.file)
dll.file <- gsub(".cpp", "", cpp.file)  
dyn.load(dynlib(dll.file))
                  
active.parameters <- c("log_sigma_immature", "log_sigma_mature", "eta") 
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(parameters, active.parameters),
                 random = c("station_immature_effect", "station_mature_effect"),
                 DLL = dll.file)

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

parameters$log_sigma_immature <- as.numeric(theta$par["log_sigma_immature"] )
parameters$log_sigma_mature <- as.numeric(theta$par["log_sigma_mature"])
parameters$eta <- as.numeric(theta$par[names(theta$par) == "eta"])

active.parameters <- setdiff(names(parameters), c("log_sigma_station_immature", "station_immature_effect", "log_sigma_station_mature", "station_mature_effect"))
obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 map = control.parameters(parameters, active.parameters),
                 random = c("station_immature_effect", "station_mature_effect"),
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

plot(log(x$carapace.width), log(x$chela.height), pch = 21, bg = "grey", cex = 0.25, xlim = c(3.5, 5.0), ylim = c(1.5, 3.5))
abline(theta$par["log_alpha_immature"], theta$par["log_beta_immature"], col = "red", lwd = 2)
abline(theta$par["log_alpha_mature"], theta$par["log_beta_mature"], col = "blue", lwd = 2)

parameters$log_sigma_immature <- as.numeric(theta$par["log_sigma_immature"] )
parameters$log_sigma_mature <- as.numeric(theta$par["log_sigma_mature"])
parameters$eta <- as.numeric(theta$par[names(theta$par) == "eta"])

cw <- seq(0, 150, len = 1000)
p <- 1 / (1 + exp(- parameters$eta[1] - parameters$eta[2] * cw))
plot(cw, p)

plot(log(x$carapace.width), log(x$chela.height), pch = 21, bg = "grey", cex = 0.25, xlim = c(3.5, 5.0), ylim = c(1.5, 3.5))
abline(theta$par["log_alpha_immature"], theta$par["log_beta_immature"], col = "red", lwd = 2)
abline(theta$par["log_alpha_mature"], theta$par["log_beta_mature"], col = "blue", lwd = 2)

fixed <- matrix(theta$par) 
rownames(fixed) <- names(theta$par)
parameters <- update.parameters(parameters, fixed = fixed)

maturity <- function(x, y, theta){
   p_kurtotic = 1 / (1 + exp(-theta$logit_p_kurtotic)) 
   p = 1 / (1 + exp(-theta$eta[1] -theta$eta[2] * x))
   mu_immature <- theta$log_alpha_immature +  theta$log_beta_immature * log(x)
   mu_mature <- theta$log_alpha_mature +  theta$log_beta_mature * log(x)
   sigma_immature <- exp(theta$log_sigma_immature)
   sigma_mature <- exp(theta$log_sigma_mature)
   sigma_kurtotic <- exp(theta$log_sigma_kurtotic)
   
   
   p_imm <- (1-p) * ((1-p_kurtotic) * dnorm(log(y), mu_immature, sigma_immature, FALSE) +                 
                       (p_kurtotic) * dnorm(log(y), mu_immature, sigma_immature + sigma_kurtotic, FALSE))
   p_mat <- (p) * ((1-p_kurtotic) * dnorm(log(y), mu_mature, sigma_mature, FALSE) +
                     (p_kurtotic) * dnorm(log(y), mu_mature, sigma_mature + sigma_kurtotic, FALSE))

   return(p_mat / (p_imm + p_mat))
}

# Perform maturity classification:
mat <- maturity(data$cw, data$y, parameters)

plot(log(data$cw), log(data$y),  xlim = c(3.5, 5.0), ylim = c(1.5, 3.5), type = "n")
index <- mat >= 0.5
points(log(data$cw[index]), log(data$y[index]), pch = 21, bg = "grey", cex = 0.25)
points(log(data$cw[!index]), log(data$y[!index]), pch = 21, bg = "red", cex = 0.25)
abline(theta$par["log_alpha_immature"], theta$par["log_beta_immature"], col = "red", lwd = 2)
abline(theta$par["log_alpha_mature"], theta$par["log_beta_mature"], col = "blue", lwd = 2)

# Proportion of matures:
res <- aggregate(list(n = mat), by = list(cw = round(data$cw)), length)
res$k <- aggregate(list(k = mat), by = list(cw = round(data$cw)), sum)$k
res$p <- res$k / res$n
windows()
dbarplot(res$p, res$cw, xlim = c(40, 145))


ri <- log(data$y[mat < 0.5]) - theta$par["log_alpha_immature"] - theta$par["log_beta_immature"] * log(data$cw[mat < 0.5])
rm <- log(data$y[mat >= 0.5]) - theta$par["log_alpha_mature"] - theta$par["log_beta_mature"] * log(data$cw[mat >= 0.5])

resi <- aggregate(ri, by = list(data$station[mat < 0.5]), mean)
resm <- aggregate(rm, by = list(data$station[mat >= 0.5]), mean)

resi[, 3] <- NA
for (i in 1:nrow(resi)){
   index <- which(resi[i, 1] == resm[,1])
   if (length(index) > 0)  resi[i,3] <- resm[index, 2]
}
plot(resi[, 2], resi[,3])
 abline(0, 0)
 lines(c(0,0), c(-1, 1))


plot(log(data$cw[mat < 0.5]), ri, ylim = c(-0.2, 0.2))
z <- log(data$cw[mat < 0.5])
z <- round(z * 10) / 10
boxplot(ri ~ z, ylim = c(-0.1, 0.1))

obj <- MakeADFun(data = data[data.cpp(cpp.file)], 
                 parameters = parameters[parameters.cpp(cpp.file)],
                 random = c("station_immature_effect", "station_mature_effect"),
                 DLL = dll.file)
theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

rep <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
parameters <- update.parameters(parameters,  fixed, random)

# Maturity proportions diagnostic plot:
mat <- maturity(data$cw, data$y, parameters)
res <- aggregate(list(n = mat), by = list(cw = round(data$cw)), length)
res$k <- aggregate(list(k = mat), by = list(cw = round(data$cw)), sum)$k
res$p <- res$k / res$n
eta <- fixed[rownames(fixed) == "eta", 1]
p <- eta[1] + eta[2] * res$cw
p <- exp(p) / (1 + exp(p))
dbarplot(res$p, res$cw, xlim = c(40, 145), ylim = c(0, 1), xlab = "Carapace width (mm)", ylab = "Proportion mature")
lines(res$cw, p, col = "red", lwd = 2)

# Logit linearity plot:
windows()
r <- log(res$p / (1-res$p))
plot(res$cw, r)

# Spatial plot
ran <- random[rownames(random) == "station_immature_effect", 1]
gulf.map(sea = TRUE)
s <- read.scset(year = unique(x$year), valid = 1)
index <- match(1:length(ran), as.numeric(substr(s$tow.id, 3, 5)))
ii <- random[, 1] >= 0
points(longitude(s)[index][ii], latitude(s)[index][ii], cex = 25*sqrt(ran[ii]), pch = 21, bg = "red")
ii <- random[, 1] < 0
points(longitude(s)[index][ii], latitude(s)[index][ii], cex = 25*sqrt(abs(ran[ii])), pch = 21, bg = "grey")

ran <- random[rownames(random) == "station_mature_effect", 1]
windows()
gulf.map(sea = TRUE)
s <- read.scset(year = unique(x$year), valid = 1)
index <- match(1:length(ran), as.numeric(substr(s$tow.id, 3, 5)))
ii <- random[, 1] >= 0
points(longitude(s)[index][ii], latitude(s)[index][ii], cex = 25*sqrt(ran[ii]), pch = 21, bg = "red")
ii <- random[, 1] < 0
points(longitude(s)[index][ii], latitude(s)[index][ii], cex = 25*sqrt(abs(ran[ii])), pch = 21, bg = "grey")


