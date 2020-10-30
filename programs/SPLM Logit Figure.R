
alpha <- -10.06351
beta <- c(0.14367326, 0.02815601, 0.23687980)
transition = c(60.98740, 98.36155)
window <- exp(1.0)
xx <- seq(40, 120, len = 1000)

windows()
m <- rbind(matrix(1, nrow = 3, ncol = 5), matrix(2, nrow = 3, ncol = 5), matrix(3, nrow = 3, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0, 0, 0, 0))

step <- function(x) as.numeric(x >= 0)
logistic <- function(x, scale = 1) 1 / (1 + exp(-x/scale))

vv <- rep(beta[1], length(xx))
ll <- vv
k <- length(transition)
for (i in 1:k) vv <- vv + (beta[i+1] - beta[i]) * step(xx - transition[i])
for (i in 1:k) ll <- ll + (beta[i+1] - beta[i]) * logistic(xx - transition[i], window)

plot(xx, vv, type = "l", lwd = 2, col = "blue", ylab = "Derivative value", xaxt = "n", cex.lab = 1.25)
lines(xx, ll, lwd = 2, col = "red")
mtext("Derivative value", 2, 2.5, cex = 1.25)
for (i in 1:k) lines(rep(transition[i], 2), par("usr")[3:4], lwd = 2, lty = "dashed", col = "black")
text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Piecewise Constant", cex = 1.5)
     
fun <- function(x){
   v <- log(1 + exp(x))
   v[v == Inf] <- x[v == Inf]
   return(v)
}

ss <- alpha + beta[1] * xx
for (i in 1:k) ss <- ss + 0.001 * (beta[i+1] - beta[i]) * fun((xx - transition[i]) / 0.001)

ii <- alpha + beta[1] * xx
for (i in 1:k) ii <- ii + window * (beta[i+1] - beta[i]) * log(1 + exp((xx - transition[i]) / window))

plot(xx, ss, type = "l", lwd = 2, col = "blue", ylab = "", xaxt = "n", cex.lab = 1.25)
lines(xx, ii, lwd = 2, col = "red")
mtext("Function value", 2, 2.5, cex = 1.25)
for (i in 1:k) lines(rep(transition[i], 2), par("usr")[3:4], lwd = 2, lty = "dashed", col = "black")
text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Piecewise Linear", cex = 1.5)

plot(xx, logistic(ss), type = "l", lwd = 2, col = "blue", ylab = "", xaxt = "n", cex.lab = 1.25, ylim = c(0, 1), yaxs = "i")
lines(xx, logistic(ii), lwd = 2, col = "red")
axis(1, cex.axis = 1.25)
mtext("x", 1, 3, cex = 1.5)
mtext("Logistic transform", 2, 2.5, cex = 1.25)
for (i in 1:k) lines(rep(transition[i], 2), par("usr")[3:4], lwd = 2, lty = "dashed", col = "black")
text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Logit-Piecewise Linear", cex = 1.5)

legend("bottomright", legend = c("Regular", "Logistic-Smoothed"), lwd = 2, col = c("blue", "red"), bg = "white", cex = 1.4)

