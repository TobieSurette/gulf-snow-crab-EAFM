library(TMB)
library(gulf.data)
library(gulf.graphics)

years  <- 1998:2020
instars <- 4:11
mu_instars <- c(10.0, 13.9, 19.4, 27.1, 37.8, 48.8, 66.8, 87.0)
names(mu_instars) <- 4:11
mu_instars <- mu_instars[as.character(instars)]

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_year.cpp")
dyn.load(dynlib("instar_year"))

# Define data:
b <- read.scsbio(years, survey = "regular")
b <- b[which(is.category(b, "MI")), ]
b <- b[which(b$carapace.width > 10), ]
b$year <- year(b)

r <- aggregate(list(f = b$year), list(x = round(log(b$carapace.width), 2), year = b$year), length)
data <- list(x = r$x, f = r$f, year = r$year - min(r$year))
   
# Define initial parameters:
parameters <- list(mu_instar_0 = as.numeric(log(mu_instars[1])),   # Mean size of the first instar.
                   log_increment = as.numeric(log(diff(log(mu_instars[1:length(instars)])))), # Vector of log-scale growth increments.
                   log_mu_increment = -1.165,          # Log-scale mean parameter associated with instar growth increments.
                   log_sigma_increment = -3.42,        # Log-scale error parameter associated with instar growth increments.
                   mu_year = rep(0, length(years)),    # Annual instar mean deviations.
                   log_sigma_mu_year = -4,             # Log-scale error for annual instar mean deviations.
                   log_sigma_mu_instar_year = -4, 
                   mu_instar_year = rep(0, length(instars) * length(years)), 
                   mu_log_sigma_instar = -2.24,        # Instar errors log-scale mean. 
                   log_sigma_log_sigma_instar = -2.44, # Instar errors log-scale error.
                   log_sigma_instar = rep(-2.41, length(instars)),   # Log-scale instar error parameters.
                   mu_logit_p = 4,                     # Instar proportions log-scale mean parameter.
                   log_sigma_logit_p_instar_year = -1, # Instar proportions log-scale error parameter.
                   logit_p_instar_year = rep(0, (length(instars)-1) * length(years))) # Multi-logit-scale parameters for instar proportions by year.
   
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("mu_year", "mu_instar_year", "log_sigma_instar", "log_increment", "logit_p_instar_year")) 
   
# Estimate parameters:
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   
# Parse output:
obj$par <- theta
rep <- sdreport(obj)
phi <- summary(rep, "random")
theta <- summary(rep, "fixed")
   

clg()
m <- kronecker(matrix(1:24, ncol = 3), matrix(1, ncol = 4, nrow = 4))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
m <- s <- p <- NULL
for (i in 1:length(years)){   
   # Plot output:
   gbarplot(data$f[data$year == i-1], data$x[data$year == i-1], border = "grey50", width = 0.01, xlim = c(2, 5), xaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   vline(obj$report()$mu[,i], col = "red", lwd = 2)
   x0 <- seq(0, 5, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:length(obj$report()$mu_instar)){
      lines(x0, 0.01 * obj$report()$p[j,i] * sum(data$f[data$year == i-1]) * dnorm(x0, obj$report()$mu[j,i], obj$report()$sigma_instar[j]), lwd = 2, col = "blue")
      d <- d + .01 * obj$report()$p[j,i] * sum(data$f[data$year == i-1]) * dnorm(x0, obj$report()$mu[j,i], obj$report()$sigma_instar[j]) 
   }
   lines(x0, d, col = "darkgreen", lwd = 2)
   
   vline(obj$report()$mu_instar, col = "green", lwd = 2)
   
   if ((i / length(years)) <= 0.5) axis(2)
   if (i %in% c(length(years)/2, length(years))) axis(1)
   if (i == length(years)/2) mtext("log(cw)", 1, 2.5, at = par("usr")[2])
   if (i == 3) mtext("Frequency", 2, 2, cex = 1.25, at = 0)
   text(par("usr")[1] + 0.1 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.85 * diff(par("usr")[3:4]), years[i], cex = 1.2)
}


mu_instar_year <- phi[grep("mu_instar_year", rownames(phi)), 1]
dim(mu_instar_year) <- c(length(years), length(instars))
dimnames(mu_instar_year) <- list(year = years, instar = instars)

image(years, instars, mu_instar_year, breaks = seq(-0.25, 0.25, len = 101), col = colorRampPalette(c("red", "white", "blue"))(100))

clg()
v <- obj$report()$mu
plot(range(years), c(10, 100), type = "n")
for (i in 1:length(years)){
   points(rep(years[i], nrow(v)), exp(v[, i]))
}


