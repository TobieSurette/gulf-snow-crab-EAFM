library(TMB)
library(gulf.data)
library(gulf.graphics)

category <- "MI"
years  <- 2019
n_instar <- 9

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_revised_tow.cpp")
dyn.load(dynlib("instar_revised_tow"))

# Define data:
b <- read.scsbio(years, survey = "regular")
b <- b[which(is.category(b, category)), ]
b <- b[which(b$carapace.width >= 2), ]
b$x <- round(b$carapace.width, 1)
b$tow <- match(b[c("date", "tow.id")], unique(b[c("date", "tow.id")])) - 1
data <- as.list(aggregate(list(f = b$tow), b[c("tow", "x")], length))

# Define initial parameters:
parameters <- list(mu0 = 10,              # First instar mean size.
                   log_sigma0 = log(0.8), # Log-scale standard error for first instar.
                   log_hiatt_slope     = log(c(0.350, 0.080)), # Hiatt slope parameters.
                   log_hiatt_intercept = log(c(0.689, 9.000)), # Hiatt intercept parameters.
                   log_growth_error    = log(c(0.05, 0.22)),   # Growth increment error inflation parameters
                   logit_p_instar      = rep(4.5, n_instar-1),
                   log_sigma_logit_p_instar_tow = -2,
                   logit_p_instar_tow = rep(0, (n_instar-1)*(max(data$tow)+1)))

# Fit proportions:
map = list(mu0 = factor(NA), 
           log_sigma0 = factor(NA),
           log_hiatt_slope = factor(c(NA, NA)),
           log_hiatt_intercept = factor(c(NA, NA)),
           log_growth_error = factor(c(NA, NA)),
           log_sigma_logit_p_instar_tow = factor(NA),
           logit_p_instar_tow = factor(rep(NA, length(parameters$logit_p_instar_tow))))
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$logit_p_instar <- as.numeric(theta)

# Add growth parameters:
map$log_hiatt_slope = factor(c(1, 2))
map$log_hiatt_intercept = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$logit_p_instar <- theta[grep("logit_p_instar", names(theta))]
parameters$log_hiatt_slope <- theta[grep("log_hiatt_slope", names(theta))]
parameters$log_hiatt_intercept <- theta[grep("log_hiatt_intercept", names(theta))]

# Add initial instar mean:
map$mu0 = factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$mu0                 <- as.numeric(theta[grep("mu0", names(theta))])
parameters$logit_p_instar      <- as.numeric(theta[grep("logit_p_instar", names(theta))])
parameters$log_hiatt_slope     <- as.numeric(theta[grep("log_hiatt_slope", names(theta))])
parameters$log_hiatt_intercept <- as.numeric(theta[grep("log_hiatt_intercept", names(theta))])

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
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

# Add tow-level proportions:
map$log_sigma_logit_p_instar_tow <- factor(1)
map$logit_p_instar_tow <- factor(1:length(parameters$logit_p_instar_tow))
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
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
parameters$log_sigma_logit_p_instar_tow <- as.numeric(theta[grep("log_sigma_logit_p_instar_tow", names(theta))])

plot.instar <- function(x, f, p, mu, sigma, xlim = c(0, 120), n_instar = 9){
   clg()
   dev.new(width = 8.5, height = 11)
   w <- min(diff(sort(unique(x))))
   gbarplot(f, x, border = "grey60", col = "grey85", xlim = xlim, width = w, xaxs = "i", xaxt = "n", yaxt = "n", lwd = 0.5)
   grid()
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:n_instar){
      lines(x0, w * p[j] * sum(f) * dnorm(x0, mu[j], sigma[j]), lwd = 1, lty = "dashed", col = "blue")
      d  <- d + w * p[j] * sum(f) * dnorm(x0, mu[j], sigma[j]) 
   }
   lines(x0, d, col = "blue", lwd = 2)
   vline(obj$report()$mu, col = "red", lwd = 1, lty = "dashed")
   mtext("Frequency", 2, 2.5, cex = 1.25)
   mtext("Carapace width(mm)", 1, 2.5, cex = 1.25)
   axis(1)
   axis(2)
   axis(3, at = obj$report()$mu, as.roman(4:12))
   box()
}

# Proportion matrix:
image(4:12,1:321, obj$report()$p, xaxt = "n",  xlab = "", ylab = "")
mtext("Instar", 1, 2.5, cex = 1.25)
mtext("Tow", 2, 2.5, cex = 1.25)
axis(1, at = 4:12, as.roman(4:12))
box()

tab <- unique(b[c("date", "tow.id")])

library(gulf.spatial)

s <- read.scsset(2019, valid = 1, survey = "regular")
index <- match(tab, s)
p <- obj$report()$p
rownames(p) <- as.character(as.roman(4:12))
colnames(p) <- substr(s$tow.id, 3, 5)[index]
rownames(s) <- substr(s$tow.id, 3, 5)

map.new()
map("coast")
points(longitude(s[colnames(p), ]), latitude(s[colnames(p), ]), cex = 8 * p[9, ] / max(p[9, ]))

r <- aggregate(data["f"], by = data["tow"], sum)
which(r[, 2] > 300)

# Individual tow plot:
i <- 167
p <- obj$report()$p[, i]
mu <- obj$report()$mu
sigma <- obj$report()$sigma
x <- rep(data$x[data$tow == i-1], data$f[data$tow == i-1])
t <- table(round(x))
plot.instar(as.numeric(names(t)), as.numeric(t), p = p, mu, sigma)


# Catch-weighted sum of mixtures:


# Map abundances:
s <- read.scsset(2019, valid = 1, survey = "regular")
rownames(s) <- substr(s$tow.id, 3, 5)
tab <- unique(b[c("date", "tow.id")])
index <- match(tab, s)
p <- obj$report()$p
rownames(p) <- as.character(as.roman(4:12))
colnames(p) <- substr(s$tow.id, 3, 5)[index]
a <- p * NA
for (i in 1:ncol(p)) a[, i] <- sum(data$f[data$tow == (i-1)]) * p[, i]
map.new()
map("coast")
points(longitude(s[colnames(a), ]), latitude(s[colnames(a), ]), cex = 8 * a["XI", ] / max(a["XI", ]))

# Aggregate mixture:
P <- apply(a, 1, sum)
t <- table(rep(data$x, data$f))
plot.instar(as.numeric(names(t)), as.numeric(t), P / sum(t), mu, sigma)





