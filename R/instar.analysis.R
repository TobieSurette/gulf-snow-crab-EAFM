library(gulf.data)
library(gulf.graphics)

# Generate Instar I
b <- read.scsbio(1998, sex = 1)
b <- b[which(!is.mature(b)), ]

gbarplot(table(round(log(b$carapace.width), 2)), width = 0.01, xlim = c(2, 5), xaxs = "i")
vline(2.7 -  0.35, lty = "dashed", col = "blue")
vline(2.7 + (0:3) * 0.33, lty = "dashed", col = "red")
vline(2.7 + 3 * 0.33 + 0.25 * 1:4, lty = "dashed", col = "green")

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
   m <- growth(theta = theta)
   s <- growth(theta = theta, error = TRUE)

   # Calculate gamma to match moments:
   phi <- sigma^2 / mu # Scale parameter.
   k <- mu^2 / sigma^2 # Shape parameter.


   # Define first mixture component:
   d <- matrix(0, nrow = k, ncol = ncol(G))
   dimnames(d) <- list(component = 1:k, colnames(G))
   d[1,] <- dnorm(as.numeric(colnames(G)), 3.15, 0.15)
   d[1,] <- d[1,] / sum(d[1,]) # Normalize.

   # Iterate other components through G:
   for (i in 1:(k-1)){
      vars <- colnames(d)[colnames(d) %in% rownames(G)]


      # Calculate corresponding gamma parameters:
      phi <- s^2 / m # Scale parameter.
      k <- m^2 / s^2 # Shape parameter.
      d[i+1, vars] <- (d[i,vars] %*% G[vars, ])[1, vars] # Growth happens here.


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

