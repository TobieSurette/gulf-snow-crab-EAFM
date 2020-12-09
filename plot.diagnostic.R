

p <- parameters

# Selectivity plot:
t <- seq(0, 140, len = 1000)
p0 <- 1 / (1 + exp(-exp(p$log_selectivity_slope[1])  * (t - p$selectivity_x50[1])));
p1 <- 1 / (1 + exp(-exp(p$log_selectivity_slope[2])  * (t - p$selectivity_x50[2])));
w <- 1 / (1 + exp(-p$logit_selectivity_proportion))
y <- w * p0 + (1-w) * p1
plot(t, y, type = "l", ylim = c(0, 1), yaxs = "i", 
     lwd = 2, col = "blue", xaxs = "i", xlab = "", ylab = "")
grid()
mtext("Carapace width(mm)", 1, 2.5, cex = 1.25) 
mtext("Selectivity", 2, 2.5, cex = 1.25) 

box()

clg()
dev.new(height = 8.5, width = 11)
m <- kronecker(matrix(1:12, ncol = 3), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, 0, m, 0),0, 0)
layout(m)
par(mar = c(0,0,0,0))

for (i in 1:length(years)){
   ix <- data$year_imm == years[i] - min(years)
   eta <- obj$report()$eta_imm[ix]
   x <- data$x_imm[ix]
   f <- data$f_imm[ix]
   plot(c(0, xlim[2]-10), c(0, 300), type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   lines(x, f, lwd = 2, col = "grey60")
   lines(x, eta, lwd = 2, col = "red")
   
   # Matures
   ix <- data$year_mat == years[i] - min(years)
   eta <- obj$report()$eta_mat[ix]
   x <- data$x_mat[ix]
   f <- data$f_mat[ix]
   lines(x, f, lwd = 2, col = "lightblue")
   lines(x, eta, lwd = 2, col = "blue")
   
   text(par("usr")[1] + 0.85 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.85 * diff(par("usr")[3:4]), years[i], cex = 1.25)
   box()
   if (i %in% c(4, 8, length(years))) axis(1)
   if (i %in% 1:4) axis(2)
   if (i == 2) mtext("Density", 2, 2.5, at = 0, cex = 1.4)
   if (i == 8) mtext("Carapace width(mm)", 1, 3, cex = 1.4)
}

