library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

category <- "MI"
years  <- 2019
n_instar <- 9
step <- 0.1

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_count.cpp")
dyn.load(dynlib("instar_count"))

plot.instar <- function(obj, data, xlim = c(0, 120)){
   mu <- obj$report()$mu
   sigma <- obj$report()$sigma
   n_instar <- length(mu)
   lambda <- exp(obj$report()$log_lambda)
   gbarplot(data$f_imm, data$x_imm, width = 0.1, border = "grey70", grid = TRUE)
   
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (k in 1:n_instar){
      p <- pnorm(x0+0.05, mu[k], sigma[k]) - pnorm(x0-0.05, mu[k], sigma[k])
      lines(x0, lambda[k] * p, lwd = 1, lty = "dashed", col = "blue")
      d  <- d + lambda[k] * p 
   }
   lines(x0, d, col = "blue", lwd = 2)
   vline(mu, col = "red", lwd = 1, lty = "dashed")
   mtext("Frequency", 2, 2.25, cex = 1.25)
   mtext("Carapace width(mm)", 1, 2.5, cex = 1.25)
   axis(1)
   axis(2)
   axis(3, at = mu, label = as.roman(4:(4+n_instar-1)))
   box()
}

# Define data:
s <- read.scsset(years, survey = "regular", valid = 1)
b <- read.scsbio(years, survey = "regular", sex = 1)
b$maturity <- morphometric.maturity(b)
b <- b[which(b$carapace.width >= 2), ]
import(s, fill = 0) <- freq(b[which(!b$maturity), ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
fi <- apply(s[fvars], 2, mean)
fi <- length(which(!b$maturity)) * (fi / sum(fi))
fi[setdiff(as.character(seq(5, 120, by = step)), names(fi))] <- 0
fi <- fi[order(as.numeric(names(fi)))]
    
s <- read.scsset(years+1, survey = "regular", valid = 1)
b <- read.scsbio(years+1, survey = "regular", sex = 1)
b$maturity <- morphometric.maturity(b)
b <- b[which(b$carapace.width >= 2), ]
ix <- which(b$maturity & is.new.shell(b))
import(s, fill = 0) <- freq(b[ix, ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
fm <- apply(s[fvars], 2, mean)
fm <- length(ix) * (fm / sum(fm))
fm[setdiff(as.character(seq(40, 140, by = step)), names(fm))] <- 0
fm <- fm[order(as.numeric(names(fm)))]


# Define data:
data <- list(x_imm = as.numeric(names(fi)),
             f_imm = as.numeric(fi))
            
# Define initial parameters:
parameters <- list(mu0                     = 10,                       # First instar mean size.
                   log_sigma0              = log(0.8),                 # Log-scale standard error for first instar.
                   log_hiatt_slope         = log(c(0.350, 0.0920)),    # Hiatt slope parameters.
                   log_hiatt_intercept     = log(c(0.689, 8.000)),     # Hiatt intercept parameters.
                   log_growth_error        = log(c(0.01, 0.10)),       # Growth increment error inflation parameters
                   log_lambda_alpha        = 7.5,                        # Log-scale global mean density.
                   log_lambda_instar       = rep(0, n_instar),         # Log-scale instar effect.
                   log_sigma_lambda_instar = -1)                       # Log-scale instar effect error parameter.

# Fit instar abundances:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))
map$log_lambda_alpha        <- factor(1)
map$log_lambda_instar       <- factor(1:n_instar)
map$log_sigma_lambda_instar <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_count", map = map, random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_count", map = map, random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 200))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$log_growth_error <- theta[grep("log_growth_error", names(theta))]
parameters$log_sigma0 <- theta[grep("log_sigma0", names(theta))]

obj <- MakeADFun(data, parameters, DLL = "instar_count", random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta

plot.instar(obj, data)





