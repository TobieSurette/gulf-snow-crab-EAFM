library(gulf.data)
library(gulf.graphics)

theta <- c(intercept = 0.689, transition = 38.2, slope = c(0.34, 0.083), window = 1.6, sigma = c(0.08, 0.25))

mu <- 10
sigma <- 0.8

for (i in 5:8){
   mu[i-3] <- 0.689 + 1.34 * mu[i-4]
   sigma[i-3] <- (1.34 + 0.05) * sigma[i-4]
}
for (i in 9:12){
   mu[i-3] <- 9 + 1.103 * mu[i-4]
   sigma[i-3] <- (1.083 + 0.26) * sigma[i-4]
}

instars <- list()
for (i in 1:9){
   instars[[i]] <- rnorm(100000, mu[i], sigma[i])
}

for (i in 1:length(instars)) instars[[i]] <- sample(instars[[i]], round(min((i-1)/7, 1) * length(instars[[i]])))
#gbarplot(table(round(unlist(instars), 1)), width = 0.1, border = "lightblue")

clg()
gbarplot(table(round(log(unlist(instars)), 2)), width = 0.01, xlim = c(2.3, 4.8), 
         xaxs = "i", border = "grey75", col = "grey90")
box()

z <- round(unlist(lapply(instars, mean)), 1)
names(z) <- as.roman(4:(length(z) + 3))
print(z)
names(instars) = as.roman(4:(length(z) + 3))

logz <- unlist(lapply(instars, function(x) mean(log(x[x > 0]))))
names(logz) <- as.roman(4:(length(z) + 3))
print(logz)
vline(logz, col = "red")

#clg()
#gbarplot(table(round(as.numeric(unlist(instars)), 1)), width = 0.01, xlim = c(0, 140), xaxs = "i", border = "grey75", col = "grey90")
#vline(z, col = "red")
#box()


#NaN 14.1 19.6 26.9 36.6 48.9 63.3 79.1 96.1
#print(round(unlist(lapply(instars, sd)), 1))


#  IV    V   VI  VII VIII   IX    X   XI  XII
#  NA  1.1  1.6  2.2  2.9  4.1  5.8  7.4  9.1
#      1.2  1.7  2.4  3.4  4.8  6.6  9.1 12.6


