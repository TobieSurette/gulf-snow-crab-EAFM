library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

years  <- 2015
n_instar <- 9
step <- 0.1

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_count.cpp")
dyn.load(dynlib("instar_count"))

plot.instar <- function(obj, data, xlim = c(0, 140), ylim = c(0, 50)){
   mu <- obj$report()$mu
   sigma <- obj$report()$sigma
   n_instar <- length(mu)
   lambda_imm <- obj$report()$lambda_imm
   lambda_mat <- obj$report()$lambda_mat 
    
   m <- kronecker(matrix(1:2), matrix(1, nrow = 5, ncol = 5))
   m <- rbind(0, cbind(0, m, 0), 0)
   layout(m)
   par(mar = c(0,0,0,0))
   
   # Immature:
   gbarplot(data$f_imm, data$x_imm, width = data$dx, border = "grey70", 
            grid = TRUE, xlim = xlim, ylim = ylim, xaxs = "i", xaxt = "n")
   grid()
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (k in 1:(n_instar-1)){
      p <- pnorm(x0+data$dx/2, mu[k], sigma[k]) - pnorm(x0-data$dx/2, mu[k], sigma[k])
      p[is.na(p)] <- 0
      lines(x0, lambda_imm[k] * p, lwd = 1, lty = "dashed", col = "blue")
      d <- d + lambda_imm[k] * p 
   }
   lines(x0, d, col = "blue", lwd = 2)
   vline(mu, col = "red", lwd = 1, lty = "dashed")
   mtext("Frequency", 2, 2.25, cex = 1.25)
   axis(2)
   axis(3, at = mu, label = as.roman(4:(4+n_instar-1)))
   text(120, ylim[1] + 0.8 * diff(ylim), years, cex = 1.5)
   box()
   mtext("Immature", 4, 1.5, cex = 1.25, srt = 180)
   
   # Mature:
   gbarplot(data$f_mat, data$x_mat, width = data$dx, border = "grey70", 
            grid = TRUE, xlim = xlim, ylim = ylim, xaxs = "i")
   grid()
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (k in 1:(n_instar-1)){
      p <- pnorm(x0+data$dx/2, mu[k+1], sigma[k+1]) - pnorm(x0-data$dx/2, mu[k+1], sigma[k+1])
      lines(x0, lambda_mat[k+1] * p, lwd = 1, lty = "dashed", col = "blue")
      p[is.na(p)] <- 0
      d <- d + lambda_mat[k+1] * p 
   }
   lines(x0, d, col = "blue", lwd = 2)
   vline(mu, col = "red", lwd = 1, lty = "dashed")
   mtext("Frequency", 2, 2.25, cex = 1.25)
   axis(2)
   text(120, ylim[1] + 0.8 * diff(ylim), years + 1, cex = 1.5)
   box()
   mtext("Carapace width(mm)", 1, 2.5, cex = 1.25)
   axis(1)
   mtext("New-shelled mature", 4, 1.5, cex = 1.25, srt = 180)
}

# Define data:
s <- read.scsset(years, survey = "regular", valid = 1)
b <- read.scsbio(years, survey = "regular", sex = 1)
b$tow.id <- tow.id(b)
b$maturity <- morphometric.maturity(b)
b <- b[which(b$carapace.width >= 2), ]
import(s, fill = 0) <- freq(b[which(!b$maturity), ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
fi <- apply(s[fvars], 2, mean)
#fi <- length(which(!b$maturity)) * (fi / sum(fi))
fi[setdiff(as.character(seq(5, 120, by = step)), names(fi))] <- 0
fi <- fi[order(as.numeric(names(fi)))]
    
s <- read.scsset(years+1, survey = "regular", valid = 1)
b <- read.scsbio(years+1, survey = "regular", sex = 1)
b$tow.id <- tow.id(b)
b$maturity <- morphometric.maturity(b)
b <- b[which(b$carapace.width >= 35), ]
ix <- which(b$maturity & is.new.shell(b))
import(s, fill = 0) <- freq(b[ix, ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
fm <- apply(s[fvars], 2, mean)
#fm <- length(ix) * (fm / sum(fm))
fm[setdiff(as.character(seq(10, 140, by = step)), names(fm))] <- 0
fm <- fm[order(as.numeric(names(fm)))]

# Define data:
data <- list(x_imm = as.numeric(names(fi)),
             f_imm = as.numeric(fi),
             x_mat = as.numeric(names(fm)),
             f_mat = as.numeric(fm))
            
# Define initial parameters:
parameters <- list(mu0                     = 10,                       # First instar mean size.
                   log_sigma0              = log(0.8),                 # Log-scale standard error for first instar.
                   log_hiatt_slope         = log(c(0.350, 0.0920)),    # Hiatt slope parameters.
                   log_hiatt_intercept     = log(c(0.689, 8.000)),     # Hiatt intercept parameters.
                   log_growth_error        = log(c(0.01, 0.10)),       # Growth increment error inflation parameters
                   log_lambda_alpha        = 7.5,                      # Log-scale global mean density.
                   log_lambda_instar       = rep(0, n_instar),         # Log-scale instar effect.
                   log_sigma_lambda_instar = -1,                       # Log-scale instar effect error parameter.
                   logit_scale             = rep(0, n_instar))

# Fit instar abundances:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))
map$log_lambda_alpha  <- factor(1)
map$log_lambda_instar <- factor(1:length(parameters$log_lambda_instar))
map$logit_scale       <- factor(1:length(parameters$logit_scale))
map$log_sigma_lambda_instar <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_count", map = map, random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$logit_scale <- obj$report()$logit_scale

plot.instar(obj, data)

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_count", map = map, random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$logit_scale <- obj$report()$logit_scale
parameters$log_growth_error <- theta[grep("log_growth_error", names(theta))]
parameters$log_sigma0 <- theta[grep("log_sigma0", names(theta))]

obj <- MakeADFun(data, parameters, DLL = "instar_count", random = "log_lambda_instar") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
obj$par <- theta

plot.instar(obj, data, ylim = c(0, 40))


100 * obj$report()$scale



