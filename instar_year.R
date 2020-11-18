library(TMB)
library(gulf.data)
library(gulf.graphics)

category <- "MI"
years  <- 1998:2020
mu_instars <- c(10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)
names(mu_instars) <- 4:11
roman <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV")
   
# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_year.cpp")
dyn.load(dynlib("instar_year"))

# Define data:
b <- read.scsbio(years, survey = "regular")
b <- b[which(is.category(b, category)), ]
b <- b[which(b$carapace.width > 6), ]
if (b$sex[1] == "2"){
   b <- b[which(b$carapace.width <= 90), ] 
   mu_instars <- mu_instars[as.character(4:9)]
}else{
   b <- b[which(b$carapace.width <= 140), ]
} 
b$year <- year(b)

# Set up data:
r <- aggregate(list(f = b$year), list(x = round(log(b$carapace.width), 2), year = b$year), length)
data <- list(x = r$x, f = r$f, year = r$year - min(r$year))
   
# Define initial parameters:
parameters <- list(mu_instar_0 = as.numeric(log(mu_instars[1])),   # Mean size of the first instar.
                   log_increment = as.numeric(log(diff(log(mu_instars)))), # Vector of log-scale growth increments.
                   log_mu_increment = -1.165,          # Log-scale mean parameter associated with instar growth increments.
                   log_sigma_increment = -3.42,        # Log-scale error parameter associated with instar growth increments.
                   log_sigma_mu_instar_year = -4, 
                   mu_instar_year = rep(0, length(mu_instars) * length(years)), 
                   beta_log_sigma_instar = c(-2.67067, 0.09134),   # Log-scale instar error parameters.
                   log_sigma_log_sigma_instar_year = -4,
                   log_sigma_instar_year = rep(0, length(mu_instars) * length(years)), 
                   mu_logit_p = 4,                     # Instar proportions log-scale mean parameter.
                   log_sigma_logit_p_instar_year = -1, # Instar proportions log-scale error parameter.
                   logit_p_instar_year = rep(0, (length(mu_instars)-1) * length(years))) # Multi-logit-scale parameters for instar proportions by year.
   
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_increment", "mu_instar_year", "log_sigma_instar_year", "logit_p_instar_year"),
                 map = list(mu_instar_0 = factor(NA), 
                            log_increment = factor(rep(NA, length(parameters$log_increment))),
                            beta_log_sigma_instar = factor(rep(NA, length(parameters$beta_log_sigma_instar))),
                            log_sigma_mu_instar_year = factor(NA),
                            mu_instar_year = factor(rep(NA, length(parameters$mu_instar_year))),
                            log_sigma_log_sigma_instar_year = factor(NA),
                            log_sigma_instar_year = factor(rep(NA, length(parameters$log_sigma_instar_year))))
                 ) 
   


# Estimate parameters:
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))
plot.instar.year(obj, data)

# Estimate parameters: 
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_increment", "mu_instar_year", "log_sigma_instar_year", "logit_p_instar_year"),
                 map = list(mu_instar_0 = factor(NA), 
                            log_increment = factor(rep(NA, length(parameters$log_increment))),
                            log_sigma_mu_instar_year = factor(NA),
                            mu_instar_year = factor(rep(NA, length(parameters$mu_instar_year))),
                            log_sigma_log_sigma_instar_year = factor(NA),
                            log_sigma_instar_year = factor(rep(NA, length(parameters$log_sigma_instar_year))))) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
obj$par <- theta
rep <- sdreport(obj)
parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))

# Estimate parameters: 
obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_increment", "mu_instar_year", "log_sigma_instar_year", "logit_p_instar_year"),
                 map = list(log_sigma_log_sigma_instar_year = factor(NA),
                            log_sigma_instar_year = factor(rep(NA, length(parameters$log_sigma_instar_year))))) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
obj$par <- theta
rep <- sdreport(obj)
parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))

plot.instar.year(obj, data)


obj <- MakeADFun(data, parameters, DLL = "instar_year", 
                 random = c("log_increment", "mu_instar_year", "log_sigma_instar_year", "logit_p_instar_year"))
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
rep <- sdreport(obj)
parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))

plot.instar.year(obj, data)

mu_instar_year <- obj$report()$mu_instar_year
dim(mu_instar_year) <- c(length(years), length(mu_instars))
dimnames(mu_instar_year) <- list(year = years, instar = names(mu_instars))

# Instar proportions:
p <- obj$report()$p
dimnames(p) <- list(instar = names(mu_instars), year = years)
image(years, as.numeric(names(mu_instars)), t(p), ylab = "Instar")

# Instar standard errors:
sigma <- obj$report()$sigma
dimnames(sigma) <- list(instar = names(mu_instars), year = years)
image(years, as.numeric(names(mu_instars)), t(sigma), ylab = "Instar")


plot(apply(mu, 1, mean), apply(sigma, 1, mean))

sigma_instar_year <- obj$report()$sigma_instar_year
dim(sigma_instar_year) <- c(length(years), length(mu_instars))
dimnames(sigma_instar_year) <- list(year = years, instar = names(mu_instars))
image(years, as.numeric(names(mu_instars)), 
      log(sigma_instar_year), breaks = seq(-0.25, 0.25, len = 101), 
      col = colorRampPalette(c("red", "white", "blue"))(100), 
      ylab = "Year")


image(years, as.numeric(names(mu_instars)), 
      mu_instar_year, breaks = seq(-0.25, 0.25, len = 101), 
      col = colorRampPalette(c("red", "white", "blue"))(100), 
      ylab = "Year")

clg()
mu <- obj$report()$mu
dimnames(mu) <- list(names(mu_instars), years)

plot(range(years), c(10, 55), type = "n", xlab = "Years", ylab = "Carapace width (mm)")
for (i in 1:length(years)){
   points(rep(years[i], nrow(mu)), exp(mu[, i]))
}


dmu <- exp(mu) - repvec(apply(exp(mu), 1, mean), ncol = ncol(mu))
clg()
colorbar(levels = seq(-3, 3, by = 1), color = colorRampPalette(c("red", "white", "blue"))(6),
         caption = c("difference(mm)"), smooth = TRUE)

image(years, as.numeric(names(mu_instars)), t(dmu), 
      breaks = seq(-3, 3, len = 101), 
      col = colorRampPalette(c("red", "white", "blue"))(100),
      xlab = "Years", ylab = "Instar")
box()





