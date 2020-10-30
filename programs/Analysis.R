library(gulf.data)
library(gulf.graphics)

clg()

years <- 2011:2020

for (i in 1:length(years)){
   b <- read.scsbio(years[i], survey = "regular", sex = 1)

   windows()
   x <- log(b$carapace.width)
   y <- log(b$chela.height)
   index <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
   x <- x[index]
   y <- y[index]
   #m <- lm(y~x)
   #r <- residuals(m)
   r <-  y - 1.366435 * x + 3.234782 
   plot(x, r, pch = 21, cex = 0.6, xlim = c(3.4, 5), ylim = c(-0.5, 0.5), xaxs = "i", yaxs = "i")
   mtext(years[i], 3, 2)
   v <- seq(10, 100, by = 10)
   vline(log(v), col = "red")
}
