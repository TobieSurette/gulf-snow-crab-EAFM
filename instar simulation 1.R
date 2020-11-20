theta <- c(intercept = 0.689, transition = 38.2, slope = c(0.34, 0.083), window = 1.6, sigma = c(0.08, 0.25))

mu <- growth(theta = theta)
sigma <- growth(theta = theta, error = TRUE)

x0 <- seq(10, 140, length = 10000)
#plot(x0, mu(x0), type = "l", lwd = 2, ylim = c(0, 40), yaxs = "i")
#lines(x0, mu(x0) + sigma(x0), lwd = 1, lty = "dashed", col = "red")
#lines(x0, mu(x0) - sigma(x0), lwd = 1, lty = "dashed", col = "red")
#lines(x0, mu(x0) + 1.96 * sigma(x0), lwd = 1, lty = "dotted", col = "red")
#lines(x0, mu(x0) - 1.96 * sigma(x0), lwd = 1, lty = "dotted", col = "red")

x <- rnorm(150000, 10, 0.8)
instars <- list(x)
for (i in 1:8){
   instars[[i+1]] <- instars[[i]] + mu(instars[[i]]) + sigma(instars[[i]]) * rnorm(length(instars[[i]]))
}

for (i in 1:length(instars)) instars[[i]] <- sample(instars[[i]], round(min((i-1)/7, 1) * length(instars[[i]])))
#gbarplot(table(round(unlist(instars), 1)), width = 0.1, border = "lightblue")

clg()
gbarplot(table(round(log(unlist(instars)), 2)), width = 0.01, xlim = c(2.3, 4.8), xaxs = "i")

z <- round(unlist(lapply(instars, mean)), 1)
names(z) <- c("IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")
print(z)
names(instars) = c("IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")

logz <- unlist(lapply(instars, function(x) mean(log(x))))
names(logz) <- c("IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")
print(logz)
vline(logz, col = "red")

print(round(unlist(lapply(instars, sd)), 1))



#9.5 15.5 20.1 27.3 36.3 49.5 66.1 90.1
