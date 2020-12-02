library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

category <- "MI"
years  <- 2019
n_instar <- 9

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar_tow.cpp")
dyn.load(dynlib("instar_tow"))

# Define data:
s <- read.scsset(years, survey = "regular", valid = 1)
b <- read.scsbio(years, survey = "regular", sex = 1)
b$maturity <- morphometric.maturity(b)

b <- b[which(!b$maturity), ]
b <- b[which(b$carapace.width >= 2), ]

b$x <- round(b$carapace.width, 1)
b$tow <- match(b[c("date", "tow.id")], s[c("date", "tow.id")]) - 1
data <- as.list(aggregate(list(f = b$tow), b[c("tow", "x")], length))
data$swept_area <- s$swept.area
n_tow <- nrow(s)

# Define initial parameters:
parameters <- list(mu0                   = 10,                       # First instar mean size.
                   log_sigma0            = log(0.8),                 # Log-scale standard error for first instar.
                   log_hiatt_slope       = log(c(0.350, 0.0920)),    # Hiatt slope parameters.
                   log_hiatt_intercept   = log(c(0.689, 8.000)),     # Hiatt intercept parameters.
                   log_growth_error      = log(c(0.01, 0.10)),       # Growth increment error inflation parameters
                   log_lambda_alpha      = 1,                        # Log-scale global mean density.
                   log_lambda_instar     = rep(0, n_instar),         # Log-scale instar effect.
                   log_lambda_tow        = rep(0, n_tow),            # Log-scale tow effect.
                   log_lambda_instar_tow = rep(0, n_instar * n_tow), # Log-scale instar x tow interaction effect.
                   log_sigma_lambda_instar = -1,                     # Log-scale instar effect error parameter.
                   log_sigma_lambda_tow = -1,                        # Log-scale tow effect error parameter.
                   log_sigma_lambda_instar_tow = -1)                 # Log-scale instar x tow interaction effect error parameter.


# Define random effects:
random <- c("log_lambda_instar", "log_lambda_tow", "log_lambda_instar_tow")

# Fit instar abundances:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))
map$log_lambda_alpha        <- factor(1)
map$log_lambda_instar       <- factor(1:n_instar)
map$log_sigma_lambda_instar <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_tow", map = map, random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar

# Fit tow abundances:
map$log_lambda_tow       <- factor(1:n_tow)
map$log_sigma_lambda_tow <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_tow", map = map, random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$log_sigma_lambda_tow <- theta[["log_sigma_lambda_tow"]]
parameters$log_lambda_tow <- obj$report()$log_lambda_tow

# Fit instar x tow interaction:
map$log_lambda_instar_tow <- factor(1:length(map$log_lambda_instar_tow))
map$log_sigma_lambda_instar_tow <- factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_tow", map = map, random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$log_sigma_lambda_tow <- theta[["log_sigma_lambda_tow"]]
parameters$log_lambda_tow <- obj$report()$log_lambda_tow
parameters$log_sigma_lambda_instar_tow <- theta[["log_sigma_lambda_instar_tow"]]
parameters$log_lambda_instar_tow <- obj$report()$log_lambda_instar_tow

# Add error parameters:
map$log_sigma0 = factor(1)
map$log_growth_error = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_tow", map = map, random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$log_sigma_lambda_tow <- theta[["log_sigma_lambda_tow"]]
parameters$log_lambda_tow <- obj$report()$log_lambda_tow
parameters$log_sigma_lambda_instar_tow <- theta[["log_sigma_lambda_instar_tow"]]
parameters$log_lambda_instar_tow <- obj$report()$log_lambda_instar_tow
parameters$log_growth_error <- theta[grep("log_growth_error", names(theta))]
parameters$log_sigma0 <- theta[grep("log_sigma0", names(theta))]

# Add growth parameters:
map$log_hiatt_slope = factor(c(1, 2))
map$log_hiatt_intercept = factor(c(1, 2))
obj <- MakeADFun(data, parameters, DLL = "instar_tow", map = map, random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$log_sigma_lambda_tow <- theta[["log_sigma_lambda_tow"]]
parameters$log_lambda_tow <- obj$report()$log_lambda_tow
parameters$log_sigma_lambda_instar_tow <- theta[["log_sigma_lambda_instar_tow"]]
parameters$log_lambda_instar_tow <- obj$report()$log_lambda_instar_tow
parameters$log_growth_error <- theta[grep("log_growth_error", names(theta))]
parameters$log_sigma0 <- theta[grep("log_sigma0", names(theta))]
parameters$log_hiatt_slope <- theta[grep("log_hiatt_slope", names(theta))]
parameters$log_hiatt_intercept <- theta[grep("log_hiatt_intercept", names(theta))]

# Add initial instar mean:
map$mu0 = factor(1)
obj <- MakeADFun(data, parameters, DLL = "instar_tow", map = map, random = random) 
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
obj$par <- theta
parameters$log_lambda_alpha <- theta[["log_lambda_alpha"]]
parameters$log_sigma_lambda_instar <- theta[["log_sigma_lambda_instar"]]
parameters$log_lambda_instar <- obj$report()$log_lambda_instar
parameters$log_sigma_lambda_tow <- theta[["log_sigma_lambda_tow"]]
parameters$log_lambda_tow <- obj$report()$log_lambda_tow
parameters$log_sigma_lambda_instar_tow <- theta[["log_sigma_lambda_instar_tow"]]
parameters$log_lambda_instar_tow <- obj$report()$log_lambda_instar_tow
parameters$log_growth_error <- theta[grep("log_growth_error", names(theta))]
parameters$log_sigma0 <- theta[grep("log_sigma0", names(theta))]
parameters$log_hiatt_slope <- theta[grep("log_hiatt_slope", names(theta))]
parameters$log_hiatt_intercept <- theta[grep("log_hiatt_intercept", names(theta))]
parameters$mu0 <- theta[["mu0"]]



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

