library(TMB)
library(gulf.data)
library(gulf.graphics)

years  <- 2009:2020

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar.cpp")
dyn.load(dynlib("instar"))

clg()
m <- kronecker(matrix(1:length(years), ncol = 2), matrix(1, ncol = 4, nrow = 4))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
m <- s <- p <- NULL
for (i in 1:length(years)){
   # Define data:
   b <- read.scsbio(years[i], survey = "regular")
   b <- b[which(is.category(b, "MI")), ]
   b <- b[which(b$carapace.width > 10), ]
   x <- table(round(log(b$carapace.width), 2))
   data <- list(x = as.numeric(names(x)), f = as.numeric(x))

   # Define initial parameters:
   parameters <- list(mu_instar_0 = 2.3, 
                      log_increment = c(-1.1, -1.1, -1.1, -1.10, -1.36, -1.16, -1.33),
                      log_mu_increment = -1.165,
                      log_sigma_increment = -3.42,
                      mu_log_sigma_instar = -2.24,   
                      log_sigma_log_sigma_instar = -2.44, 
                      log_sigma_instar = rep(-2.41, 8),
                      logit_p_instar = c( 3.08, 4.44, 4.92, 5.20, 6.45, 6.08, 5.52 ))

   obj <- MakeADFun(data, parameters, DLL = "instar", random = c("log_sigma_instar", "log_increment")) 

   # Estimate parameters:
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   
   # Parse output:
   obj$par <- theta
   rep <- sdreport(obj)
   phi <- summary(rep, "random")
   theta <- summary(rep, "fixed")

   # Plot output:
   gbarplot(data$f, data$x, border = "grey50", width = 0.01, xlim = c(2, 5), xaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   vline(obj$report()$mu_instar, col = "red", lwd = 2)
   x0 <- seq(0, 5, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:length(obj$report()$mu_instar)){
      lines(x0, 0.01 * obj$report()$p_instar[j] * sum(data$f) * dnorm(x0, obj$report()$mu_instar[j], obj$report()$sigma_instar[j]), lwd = 2, col = "blue")
      d <- d + .01 * obj$report()$p_instar[j] * sum(data$f) * dnorm(x0, obj$report()$mu_instar[j], obj$report()$sigma_instar[j]) 
   }
   lines(x0, d, col = "darkgreen", lwd = 2)
   
   if ((i / length(years)) <= 0.5) axis(2)
   if (i %in% c(length(years)/2, length(years))) axis(1)
   if (i == length(years)/2) mtext("log(cw)", 1, 2.5, at = par("usr")[2])
   if (i == 3) mtext("Frequency", 2, 2, cex = 1.25, at = 0)
   text(par("usr")[1] + 0.1 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.85 * diff(par("usr")[3:4]), years[i], cex = 1.2)

   # Store output:
   m <- rbind(m, obj$report()$mu_instar)
   s <- rbind(s, obj$report()$sigma_instar)
   p <- rbind(p, obj$report()$p_instar)
}
rownames(m) <- years
rownames(s) <- years
rownames(p) <- years

#clg()
#v <- obj$report()$sigma
#plot(range(years), c(10, 120), type = "n")
#for (i in 1:length(years)){
#   points(rep(years[i], ncol(mu)), exp(mu[i, ]))
#}


           