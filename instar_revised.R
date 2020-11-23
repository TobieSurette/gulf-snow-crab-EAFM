library(TMB)
library(gulf.data)
library(gulf.graphics)

years  <- 2009:2020

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_revised.cpp")
dyn.load(dynlib("instar_revised"))

   # Define data:
   b <- read.scsbio(years, survey = "regular")
   b <- b[which(is.category(b, "MI")), ]
   b <- b[which(b$carapace.width >= 2), ]
   x <- table(round(b$carapace.width, 1))
   data <- list(x = as.numeric(names(x)), f = as.numeric(x), n_instar = 9)
   
   # Define initial parameters:
   parameters <- list(mu0 = 10,              # First instar mean size.
                      log_sigma0 = log(0.8), # Log-scale standard error for first instar.
                      log_hiatt_slope     = log(c(0.350, 0.080)), # Hiatt slope parameters.
                      log_hiatt_intercept = log(c(0.689, 9.000)), # Hiatt intercept parameters.
                      log_growth_error    = log(c(0.05, 0.22)),   # Growth increment error inflation parameters
                      logit_p_instar      = rep(4.5, data$n_instar-1))
   
   # Fit proportions:
   obj <- MakeADFun(data, parameters, 
                    DLL = "instar_revised", 
                    map = list(mu0 = factor(NA), 
                               log_sigma0 = factor(NA),
                               log_hiatt_slope = factor(c(NA, NA)),
                               log_hiatt_intercept = factor(c(NA, NA)),
                               log_growth_error = factor(c(NA, NA)))) 
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   obj$par <- theta
   parameters$logit_p_instar <- as.numeric(theta)
   
   # Add growth parameters:
   obj <- MakeADFun(data, parameters, 
                    DLL = "instar_revised", 
                    map = list(mu0 = factor(NA), 
                               log_sigma0 = factor(NA),
                               log_hiatt_slope = factor(c(1, 2)),
                               log_hiatt_intercept = factor(c(1, 2)),
                               log_growth_error = factor(c(NA, NA)))) 
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   obj$par <- theta
   parameters$logit_p_instar <- theta[grep("logit_p_instar", names(theta))]
   parameters$log_hiatt_slope <- theta[grep("log_hiatt_slope", names(theta))]
   parameters$log_hiatt_intercept <- theta[grep("log_hiatt_intercept", names(theta))]
   
   # Add initial instar mean:
   obj <- MakeADFun(data, parameters, 
                    DLL = "instar_revised", 
                    map = list(mu0 = factor(1), 
                               log_sigma0 = factor(NA),
                               log_hiatt_slope = factor(c(1, 2)),
                               log_hiatt_intercept = factor(c(1, 2)),
                               log_growth_error = factor(c(NA, NA)))) 
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   obj$par <- theta
   parameters$mu0                 <- as.numeric(theta[grep("mu0", names(theta))])
   parameters$logit_p_instar      <- as.numeric(theta[grep("logit_p_instar", names(theta))])
   parameters$log_hiatt_slope     <- as.numeric(theta[grep("log_hiatt_slope", names(theta))])
   parameters$log_hiatt_intercept <- as.numeric(theta[grep("log_hiatt_intercept", names(theta))])
   
   # Add error parameters:
   obj <- MakeADFun(data, parameters, 
                    DLL = "instar_revised", 
                    map = list(mu0 = factor(1), 
                               log_sigma0 = factor(1),
                               log_hiatt_slope = factor(c(1, 2)),
                               log_hiatt_intercept = factor(c(1, 2)),
                               log_growth_error = factor(c(1, 2)))) 
   for (i in 1:10){
      theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
      obj$par <- theta
   }
   parameters$mu0                 <- as.numeric(theta[grep("mu0", names(theta))])
   parameters$log_sigma0          <- as.numeric(theta[grep("log_sigma0", names(theta))])
   parameters$logit_p_instar      <- as.numeric(theta[grep("logit_p_instar", names(theta))])
   parameters$log_hiatt_slope     <- as.numeric(theta[grep("log_hiatt_slope", names(theta))])
   parameters$log_hiatt_intercept <- as.numeric(theta[grep("log_hiatt_intercept", names(theta))])
   parameters$log_growth_error    <- as.numeric(theta[grep("log_growth_error", names(theta))])
   
   clg()
   dev.new(width = 8.5, height = 11)
   m <- s <- p <- NULL
   
   gbarplot(data$f, data$x, border = "grey60", col = "grey85", 
            width = 0.01, xlim = c(0, 120), 
            xaxs = "i", xaxt = "n", yaxt = "n", lwd = 0.5)
   grid()
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:data$n_instar){
      lines(x0, 0.1 * obj$report()$p[j] * sum(data$f) * dnorm(x0, obj$report()$mu[j], obj$report()$sigma[j]), lwd = 1, lty = "dashed", col = "blue")
      d  <- d + 0.1 * obj$report()$p[j] * sum(data$f) * dnorm(x0, obj$report()$mu[j], obj$report()$sigma[j]) 
   }
   lines(x0, d, col = "blue", lwd = 2)
   
   vline(obj$report()$mu, col = "red", lwd = 1, lty = "dashed")
   mtext("Frequency", 2, 2.5, cex = 1.25)
   mtext("Carapace width(mm)", 1, 2.5, cex = 1.25)
   axis(1)
   axis(2)
   box()

   p <- obj$report()$p
   mu <- obj$report()$mu
   sigma <- obj$report()$sigma
   
   axis(3, at = mu, as.roman(4:12))
   