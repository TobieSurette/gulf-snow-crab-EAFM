library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

category <- "MI"
years  <- 2010:2019
n_instar <- 9

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_revised_tow.cpp")
dyn.load(dynlib("instar_revised_tow"))

plot.instar <- function(x, f, p, mu, sigma, xlim = c(0, 120), n_instar = 9){
   n_instar <- length(mu)
   clg()
   dev.new(width = 8.5, height = 11)
   w <- min(diff(sort(unique(x))))
   gbarplot(f, x, border = "grey60", col = "grey85", xlim = xlim, width = w, xaxs = "i", xaxt = "n", yaxt = "n", lwd = 0.5)
   grid()
   x0 <- seq(0, 140, len = 1000)
   d <- rep(0, length(x0))
   for (j in 1:n_instar){
      print(c(j, p[j]))
      lines(x0, w * p[j] * sum(f) * dnorm(x0, mu[j], sigma[j]), lwd = 1, lty = "dashed", col = "blue")
      d  <- d + w * p[j] * sum(f) * dnorm(x0, mu[j], sigma[j]) 
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
b <- read.scsbio(2019, survey = "regular", sex = 1)
b$maturity <- morphometric.maturity(b)

b <- b[which(is.category(b, category)), ]
b <- b[which(b$carapace.width >= 2), ]

b$x <- round(b$carapace.width, 1)
b$tow <- match(b[c("date", "tow.id")], s[c("date", "tow.id")]) - 1
data <- as.list(aggregate(list(f = b$tow), b[c("tow", "x")], length))
data$swept_area <- s$swept.area

# Define initial parameters:
parameters <- list(mu0 = 10,              # First instar mean size.
                   log_sigma0 = log(0.8), # Log-scale standard error for first instar.
                   log_hiatt_slope     = log(c(0.350, 0.0920)), # Hiatt slope parameters.
                   log_hiatt_intercept = log(c(0.689, 8.000)), # Hiatt intercept parameters.
                   log_growth_error    = log(c(0.01, 0.10)),   # Growth increment error inflation parameters
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

t <- table(rep(data$x, data$f))
plot.instar(as.numeric(names(t)), as.numeric(t), 
            p = apply(obj$report()$p, 1, mean), 
            mu = obj$report()$mu, 
            sigma = obj$report()$sigma)

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
obj$par <- theta
parameters$log_sigma0          <- as.numeric(theta[grep("log_sigma0", names(theta))])
parameters$logit_p_instar      <- as.numeric(theta[grep("logit_p_instar", names(theta))])
parameters$log_growth_error    <- as.numeric(theta[grep("log_growth_error", names(theta))])

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
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 100))$par
obj$par <- theta
parameters$mu0                 <- as.numeric(theta[grep("mu0", names(theta))])
parameters$logit_p_instar      <- as.numeric(theta[grep("logit_p_instar", names(theta))])
parameters$log_hiatt_slope     <- as.numeric(theta[grep("log_hiatt_slope", names(theta))])
parameters$log_hiatt_intercept <- as.numeric(theta[grep("log_hiatt_intercept", names(theta))])

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_revised_tow", map = map, random = "logit_p_instar_tow") 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
obj$par <- theta
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
for (i in 1:1){
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
   obj$par <- theta
}
parameters$mu0                 <- as.numeric(theta[grep("mu0", names(theta))])
parameters$log_sigma0          <- as.numeric(theta[grep("log_sigma0", names(theta))])
parameters$logit_p_instar      <- as.numeric(theta[grep("logit_p_instar", names(theta))])
parameters$log_hiatt_slope     <- as.numeric(theta[grep("log_hiatt_slope", names(theta))])
parameters$log_hiatt_intercept <- as.numeric(theta[grep("log_hiatt_intercept", names(theta))])
parameters$log_growth_error    <- as.numeric(theta[grep("log_growth_error", names(theta))])
parameters$log_sigma_logit_p_instar_tow <- as.numeric(theta[grep("log_sigma_logit_p_instar_tow", names(theta))])



# Proportion matrix:
image(4:13,1:nrow(s), obj$report()$p, xaxt = "n",  xlab = "", ylab = "")
mtext("Instar", 1, 2.5, cex = 1.25)
mtext("Tow", 2, 2.5, cex = 1.25)
axis(1, at = 4:13, as.roman(4:12))
box()


p <- obj$report()$p
rownames(p) <- as.character(as.roman(4:13))
colnames(p) <- substr(s$tow.id, 3, 5)


map.new()
map("coast")
points(longitude(s), latitude(s), cex = 5 * p["VI", ] / max(p["VI", ]))

r <- aggregate(data["f"], by = data["tow"], sum)
r[which(r[, 2] > 300), 1]

# Individual tow plot:
i <- 68
p <- obj$report()$p[, i]
mu <- obj$report()$mu
sigma <- obj$report()$sigma
x <- rep(data$x[data$tow == i-1], data$f[data$tow == i-1])
t <- table(round(x))
plot.instar(as.numeric(names(t)), as.numeric(t), p = p, mu, sigma)

# Catch-weighted sum of mixtures:

# Map instar abundances:
clg()
dev.new(height = 11, width = 8.5)
m <- kronecker(matrix(1:10, ncol = 2), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0,0,cbind(0,0,m,0),0,0)
layout(m)
par(mar = c(0,0,0,0))
p <- obj$report()$p
rownames(p) <- as.character(as.roman(4:12))
colnames(p) <- substr(s$tow.id, 3, 5)
a <- p * NA
for (i in 1:ncol(p)) a[, i] <- sum(data$f[data$tow == (i-1)]) * p[, i]
for (i in 2:n_instar){
   map.new()
   grid()
   map("coast")
   points(longitude(s), latitude(s), cex = 2.5 * sqrt(a[i, ] / max(a)), lwd = 0.5)
   box(lwd = 0.5)
   if ((i-1) %in% 1:4) map.axis(2)
   if ((i-1) %in% c(4,8)) map.axis(1)
   if (i == 2) mtext(years, 3, 1, at = par("usr")[2], cex = 1.5)
   text(-61, 48.75, rownames(p)[i], cex = 1.5)
}

# Aggregate mixture:
P <- apply(a, 1, sum)
t <- table(rep(data$x, data$f))
plot.instar(as.numeric(names(t)), as.numeric(t), P / sum(t), mu, sigma)

# Make predictions!
mu <- obj$report()$mu
sigma <- obj$report()$sigma
mu[10] = exp(parameters$log_hiatt_intercept[2]) + mu[9] + exp(parameters$log_hiatt_slope[2]) * mu[9]
sigma[10] = exp(log(1 + exp(parameters$log_hiatt_slope[2]) + exp(parameters$log_growth_error[2])) + log(sigma[9]))
names(mu) <- as.roman(4:13)
names(sigma) <- as.roman(4:13)

a <- obj$report()$abundance_instar

b <- read.scsbio(years + 1, survey = "regular")
b <- b[which(is.category(b, "MM")), ]
b <- b[which(b$carapace.width >= 2), ]

# Determine scale length-frequency for mature males in the following year:
s <- read.scsset(years + 1, valid = 1, survey = "regular")
import(s, fill = 0, by = "tow.id") <- freq(b, by = "tow.id")
fvars <- names(s)[gsub("[0-9]", "", names(s)) == ""]  
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
gbarplot(apply(s[fvars], 2, mean), width = 1, border = "grey50")

# Do weighted average of for the predicted mixture:
x0 <- seq(0, 140, by = 1)
pp <- matrix(NA, nrow = length(mu)-1, ncol = length(x0))
for (i in 1:length(a)){
   pp[i, ] <- a[i] * dnorm(x0, mu[i+1], sigma[i+1])
}

M <- 0.5
mat <- c(0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.3, 0.5, 1)
gbarplot(apply(s[fvars], 2, mean), width = 1, border = "grey50")
for (i in 1:nrow(pp)){
   lines(x0, mat[i] * M * pp[i, ])
}

