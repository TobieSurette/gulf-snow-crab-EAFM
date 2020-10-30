
alpha <- 0
beta <- c(0, 1, -1)
transition = c(1, 2)
window <- 0.08


windows(width = 10)
m <- rbind(matrix(1, nrow = 3, ncol = 5), matrix(2, nrow = 3, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0, 0, 0, 0))

step <- function(x) as.numeric(x >= 0)
logistic <- function(x, scale = 1) 1 / (1 + exp(-x/scale))
xx <- seq(0, 3, len = 1000)
vv <- rep(beta[1], length(xx))
ll <- vv
k <- length(transition)
for (i in 1:k) vv <- vv + (beta[i+1] - beta[i]) * step(xx - transition[i])
for (i in 1:k) ll <- ll + (beta[i+1] - beta[i]) * logistic(xx - transition[i], window)


plot(xx, vv, type = "l", lwd = 1, col = "black", ylab = "Derivative value", xaxt = "n", cex.lab = 1.25)
lines(xx, ll, lwd = 2, col = "red")
mtext("Derivative value", 2, 2.5, cex = 1.25)
for (i in 1:k) lines(rep(transition[i], 2), par("usr")[3:4], lwd = 1, lty = "dashed", col = "black")

fun <- function(x){
   v <- log(1 + exp(x))
   v[v == Inf] <- x[v == Inf]
   return(v)
}

ss <- alpha + beta[1] * xx
for (i in 1:k) ss <- ss + 0.001 * (beta[i+1] - beta[i]) * fun((xx - transition[i]) / 0.001)

ii <- alpha + beta[1] * xx
for (i in 1:k) ii <- ii + window * (beta[i+1] - beta[i]) * log(1 + exp((xx - transition[i]) / window))

plot(xx, ss, type = "l", lwd = 1, col = "black", ylab = "", xaxt = "n", cex.lab = 1.25)
lines(xx, ii, lwd = 2, col = "red")
axis(1)
mtext("x", 1, 2.5, cex = 1.25)
mtext("Function value", 2, 2.5, cex = 1.25)
for (i in 1:k) lines(rep(transition[i], 2), par("usr")[3:4], lwd = 1, lty = "dashed", col = "black")




