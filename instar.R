library(TMB)
library(gulf.data)
library(gulf.graphics)

years  <- 1998:2009
#Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")
compile("instar.cpp")
dyn.load(dynlib("instar"))

clg()
m <- kronecker(matrix(1:length(years), ncol = 2), matrix(1, ncol = 4, nrow = 4))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
m <- s <- p<- NULL
for (i in 1:length(years)){
   # Define data:
   b <- read.scsbio(years[i], survey = "regular")
   b <- b[which(is.category(b, "MI")), ]
   b <- b[which(b$carapace.width > 10), ]
   x <- table(round(log(b$carapace.width), 2))
   data <- list(x = as.numeric(names(x)), n = as.numeric(x))

   # Define initial parameters:
   parameters <- list(log_mu_0 = 0.88,
                      log_sigma_mu = -2.41,
                      log_sigma_sigma = -2,
                      log_sigma = rep(-2.41, 8),
                      log_mu_inc = -1.25,
                      log_sigma_inc = -2.3,
                      log_inc = c(-1.21, -1.22, -1.27, -1.20, -1.36, -1.16, -1.33),
                      log_p = c(3.86, 4.96, 5.09, 5.44, 5.65, 5.92, 6.00))

   obj <- MakeADFun(data, parameters, DLL = "instar", random = c("log_inc", "log_sigma")) 
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   obj$par <- theta

   rep <- sdreport(obj)
   phi <- summary(rep, "random")
   theta <- summary(rep, "fixed")

   gbarplot(data$n, data$x, border = "grey50", width = 0.01, xlim = c(2, 5), xaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   vline(obj$report()$mu, col = "red", lwd = 2)
   x0 <- seq(0, 5, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:length(obj$report()$mu)){
      lines(x0, 0.01 * obj$report()$p[j] * sum(data$n) * dnorm(x0, obj$report()$mu[j], obj$report()$sigma[j]), lwd = 2, col = "blue")
      d <- d + .01 * obj$report()$p[j] * sum(data$n) * dnorm(x0, obj$report()$mu[j], obj$report()$sigma[j]) 
   }
   lines(x0, d, col = "darkgreen", lwd = 2)
   
   if ((i / length(years)) <= 0.5) axis(2)
   if (i %in% c(length(years)/2, length(years))) axis(1)
   if (i == length(years)/2) mtext("log(cw)", 1, 2.5, at = par("usr")[2])
   if (i == 3) mtext("Frequency", 2, 2, cex = 1.25, at = 0)
   text(par("usr")[1] + 0.1 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.85 * diff(par("usr")[3:4]), years[i], cex = 1.2)
   axis(2)
   
   # Store output:
   m <- rbind(m, obj$report()$m)
   s <- rbind(s, obj$report()$s)
   p <- rbind(p, obj$report()$p)
}
rownames(mu) <- years

clg()
v <- obj$report()$sigma
plot(range(years), c(10, 120), type = "n")
for (i in 1:length(years)){
   points(rep(years[i], ncol(mu)), exp(mu[i, ]))
}


           