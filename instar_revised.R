library(TMB)
library(gulf.data)
library(gulf.graphics)

years  <- 2009:2020

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_revised.cpp")
dyn.load(dynlib("instar_revised"))

clg()
m <- kronecker(matrix(1:length(years), ncol = 2), matrix(1, ncol = 4, nrow = 4))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
m <- s <- p <- NULL
for (i in 11:11){
   # Define data:
   b <- read.scsbio(years[i], survey = "regular")
   b <- b[which(is.category(b, "MI")), ]
   b <- b[which(b$carapace.width > 10), ]
   x <- table(round(b$carapace.width, 1))
   data <- list(x = as.numeric(names(x)), f = as.numeric(x), n_instar = 9)
   
   # Define initial parameters:
   parameters <- list(mu0 = 10,              # First instar mean size.
                      log_sigma0 = log(0.8), # Log-scale standard error for first instar.
                      log_hiatt_slope     = log(c(0.350, 0.103)), # Hiatt slope parameters.
                      log_hiatt_intercept = log(c(0.689, 9.000)), # Hiatt intercept parameters.
                      log_growth_error    = log(c(0.05, 0.26)),   # Growth increment error inflation parameters
                      logit_p_instar      = rep(4.5, data$n_instar-1))
   
   obj <- MakeADFun(data, parameters, 
                    DLL = "instar_revised", 
                    random = "logit_p_instar",
                    map = list(mu0 = factor(NA), 
                               log_sigma0 = factor(NA),
                               log_hiatt_slope = factor(c(1, 1)),
                               log_hiatt_intercept = factor(c(NA, NA)),
                               log_growth_error = factor(c(NA, NA)))) 
   
   # Estimate parameters:
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   rep <- sdreport(obj)
   
   # Parse output:
   obj$par <- theta
   theta <- summary(rep, "fixed")
   
   clg()
   dev.new(width = 8.5, height = 11)
   m <- s <- p <- NULL
   
   gbarplot(data$f, data$x, border = "grey50", width = 0.01, xlim = c(0, 120), xaxs = "i", xaxt = "n", yaxt = "n", lwd = 0.5)
   grid()
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:data$n_instar){
      lines(x0, 0.1 * obj$report()$p[j] * sum(data$f) * dnorm(x0, obj$report()$mu[j], obj$report()$sigma[j]), lwd = 0.5, lty = "dashed", col = "blue")
      d  <- d + 0.1 * obj$report()$p[j] * sum(data$f) * dnorm(x0, obj$report()$mu[j], obj$report()$sigma[j]) 
   }
   lines(x0, d, col = "blue", lwd = 0.5)
   
   # 
   vline(obj$report()$mu, col = "red", lwd = 0.5)
   mtext("Frequency", 2, 2.5, cex = 1.25, at = 0)
   mtext("Carapace width(mm)", 1, 2.5, cex = 1.25)
   box()
}

