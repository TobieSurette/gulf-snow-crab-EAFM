library(gulf.data)
library(gulf.graphics)

# Generate Instar I
b <- read.scsbio(2020, sex = 1)
b <- b[which(!is.mature(b)), ]

theta <- c(log_intercept = log(0.25),
           transition = 38.200,
           log_slope = log(c(0.349, 0.126)),
           log_window = log(1.6),
           log_sigma = -2.0)

x <- seq(0, 140, by = 0.5)
f <- dnorm(x, 3.15, 0.15)
names(f) <- x
G <- growth.matrix(as.numeric(names(f)), theta = theta)
ff <- freq(b, step = 0.5)
gbarplot(ff, width = 0.5, border = "grey50")
for (i in 1:11){
   f <- (f %*% G[names(f), ])[1, ]
   f <- f[names(f) %in% rownames(G)]

   f <-  ff[which.max(f)] * f / max(f)

   lines(as.numeric(names(f)), f, lwd = 2, col = "red")
   vline(as.numeric(names(f)[which.max(f)]), lty = "dashed", col = "red")
}

x <- freq(b, step = 0.5)

loglike <- function(theta, x, fixed){
   # Convert to frequency table:
   if (is.null(x)) x <- table(x)

   # Parameter transform:
   if (!missing(fixed)) theta <- c(theta, fixed)
   names(theta) <- tolower(names(theta))
   ix <- grep("^log_", names(theta))
   theta[ix] <- exp(theta[ix])
   names(theta) <- gsub("^log_", "", names(theta))

   # Parse mixture proportions:
   p <- theta[grep("^p[0-9]", names(theta))]
   p <- exp(p) / (1 + sum(exp(p)))
   p <- c(1-sum(p), p)

   # Number of instars:
   k <- length(p)

   # Define growth matrix:
   G <- growth.matrix(c(1, as.numeric(names(x))), theta = theta)

   # Define first mixture component:
   d <- matrix(0, nrow = k, ncol = ncol(G))
   dimnames(d) <- list(component = 1:k, colnames(G))
   d[1,] <- dnorm(as.numeric(colnames(G)), 3.15, 0.15)
   d[1,] <- d[1,] / sum(d[1,]) # Normalize.

   # Iterate other components through G:
   for (i in 1:(k-1)){
      vars <- colnames(d)[colnames(d) %in% rownames(G)]
      d[i+1, vars] <- (d[i,vars] %*% G[vars, ])[1, vars]
      d[i+1, ] <- d[i+1, ] / sum(d[i+1, ])
   }

   # Mixture log-likelihood:
   v <- as.numeric(x) * log(p %*% d[, names(x)])

   return(-sum(v))
}

plot.instar <- function(x, theta){
   # Parse mixture proportions:
   p <- theta[grep("^p[0-9]", names(theta))]
   p <- exp(p) / (1 + sum(exp(p)))
   p <- c(1-sum(p), p)

   # Number of instars:
   k <- length(p)

   # Define growth matrix:
   G <- growth.matrix(c(1, as.numeric(names(x))), theta = theta)

   # Define first mixture component:
   d <- matrix(0, nrow = k, ncol = ncol(G))
   dimnames(d) <- list(component = 1:k, colnames(G))
   d[1,] <- dnorm(as.numeric(colnames(G)), 3.15, 0.15)
   d[1,] <- d[1,] / sum(d[1,]) # Normalize.

   # Iterate other components through G:
   for (i in 1:(k-1)){
      vars <- colnames(d)[colnames(d) %in% rownames(G)]
      d[i+1, vars] <- (d[i,vars] %*% G[vars, ])[1, vars]
      d[i+1, ] <- d[i+1, ] / sum(d[i+1, ])
   }

   gbarplot(x, width = 0.5, border = "grey50")
   v <- (p %*% d[, names(x)])[1, ]
   lines(as.numeric(names(v)), sum(x)*v, col = "red")
   for (i in 1:k){
      lines(as.numeric(names(v)), p[i] * sum(x) * d[i, names(x)], col = "red")
   }

}

theta <- c(p = c(-2.28,-2.84,-1.88,-1.65,-0.04,1.55,2.68,1.82,2.34,2.47,2.97),   # Component proportions.
           log_intercept = log(0.25),
           transition = 38.200,
           log_slope = log(c(0.349, 0.126)),
           log_window = log(1.6),
           log_sigma = c(-2.5, -1.5))

loglike(theta, x)

v <- growth(1:120, theta = theta, error = TRUE)

plot(1:120, v$mu, type = "l", xlim = c(0, 120), xaxs = "i", yaxs = "i")
lines(1:120, v$mu - v$sigma, lty = "dashed")
lines(1:120, v$mu + v$sigma, lty = "dashed")


growth(1:120, theta = theta)

G <- growth.matrix(1:120, theta = theta)
image(G)

theta["log_sigma"] <- -1.8
plot.instar(x, theta)


# Estimate proportions:
fixed <- theta[-grep("^p", names(theta))]
theta <- theta[-which(names(theta) %in% names(fixed))]
loglike(theta, x, fixed = fixed)
theta <- optim(theta, loglike, x = x, fixed = fixed, control = list(trace = 3))$par
theta <- c(theta, fixed)

# Estimate error parameters:
fixed <- theta[-grep("sigma", names(theta))]
theta <- theta[-which(names(theta) %in% names(fixed))]
loglike(theta, x, fixed = fixed)
theta <- optim(theta, loglike, x = x, fixed = fixed, control = list(trace = 3))$par
theta <- c(theta, fixed)

# Fit growth parameters:
fixed <- theta[c(grep("^p", names(theta)), grep("intercept", names(theta)),
                 grep("sigma", names(theta)), grep("window", names(theta)))]
theta <- theta[-which(names(theta) %in% names(fixed))]
loglike(theta, x, fixed = fixed)
theta <- optim(theta, loglike, x = x, fixed = fixed, control = list(trace = 3, maxit = 5000))$par
theta <- c(theta, fixed)


theta <- optim(theta, loglike, x = x, control = list(trace = 3, maxit = 5000))$par

