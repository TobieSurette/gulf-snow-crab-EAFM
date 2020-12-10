
# Selectivity plot:
# female trawl selectivity.png
p <- parameters
t <- seq(0, 140, len = 1000)
p0 <- 1 / (1 + exp(-exp(p$log_selectivity_slope[1])  * (t - p$selectivity_x50[1])));
p1 <- 1 / (1 + exp(-exp(p$log_selectivity_slope[2])  * (t - p$selectivity_x50[2])));
w <- 1 / (1 + exp(-p$logit_selectivity_proportion))
y <- w * p0 + (1-w) * p1
plot(t, y, type = "l", ylim = c(0, 1), yaxs = "i", 
     lwd = 2, col = "blue", xaxs = "i", xlab = "", ylab = "", xlim = xlim)
grid()
mtext("Carapace width(mm)", 1, 2.5, cex = 1.25) 
mtext("Trawl selectivity", 2, 2.5, cex = 1.25) 
box()

# female model and length-frequencies 2010-2020.pdf
clg()
p <- parameters
dev.new(height = 11, width = 8.5)
m <- kronecker(matrix(1:n_year, ncol = 3), matrix(1, ncol = 5, nrow = 5))
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
   if (i %in% ((n_year/3)*(1:3))) axis(1)
   if (i %in% 1:(n_year/3)) axis(2)
   if (i == 3) mtext("Density", 2, 2.5, cex = 1.4)
   if (i == round(2*n_year/3)) mtext("Carapace width(mm)", 1, 3, cex = 1.4)
}


# female moulting probability VIII to IX.png
p_mat <- obj$report()$p_mat
p_mat <- p_mat[-nrow(p_mat), -ncol(p_mat)]
gbarplot(p_mat[5, ], years[-length(years)], xlab = "", ylab =  "")
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Moult-to-maturity probability", 2, 2.5, cex = 1.25)
mtext("Instar VIII to IX", 3, 1.0, cex = 1.5)


# Length-frequency data
clg()
dev.new(height = 8.5, width = 11)
m <- kronecker(matrix(1:2, ncol = 1), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, m, 0),0,0)
layout(m)
par(mar = c(0,0,0,0))
image(years, seq(40, 80, by = 0.5), fm[, as.character(seq(40, 80, by = 0.5))], 
      col = grey(seq(1, 0, len = 100)))
text(par("usr")[1] + 0.88 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Mature", cex = 1.5, font = 2)
axis(4, at = apply(obj$report()$mu, 1, mean), labels = as.character(as.roman(4:(n_instar+3))))
box()
mtext("Carapace width (mm)", 2, 3, at = par("usr")[3], cex = 1.25)
image(years, seq(5, 75, by = 0.5), fi[, as.character(seq(5, 75, by = 0.5))], 
      col = grey(seq(1, 0, len = 100)), yaxt = "n")
mtext("Year", 1, 3, cex = 1.25)
axis(2, at = seq(10, 70, by = 10))
text(par("usr")[1] + 0.88 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Immature", cex = 1.5, font = 2)
box()
axis(4, at = apply(obj$report()$mu, 1, mean), labels = as.character(as.roman(4:(n_instar+3))))



