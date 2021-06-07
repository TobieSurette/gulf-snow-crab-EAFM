
# Immatures
clg()
plot(c(0, ncol(mu)), c(0, 75), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
for (i in 1:nrow(mu)){
   xx <- 1:ncol(mu)
   yy <- exp(mu[i, ])
   col <- rainbow(nrow(mu))[i]
   ix <- which(abs(yy - exp(r$mu_instar[i])) > 0.05)
   points(xx[ix], yy[ix], pch = 21, bg = col, cex = 0.5)
   hline(exp(r$mu_instar[i]), col = col, lwd = 2)
}
mtext("Instar mean size (mm)", 2, 2.5, cex = 1.25)
mtext("Tow", 1, 2.5, cex = 1.25)
axis(4, at = exp(r$mu_instar), labels = as.roman(as.numeric(names(mu_instars))))
delta <- c(0, which(diff(year(tows))>0), nrow(tows))
vline(delta, lwd = 2)
axis(3, at = delta[1:(length(delta)-1)] + diff(delta)/2, years)
box()

# Matures
clg()
plot(c(0, ncol(mu)), c(0, 75), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
for (i in 1:nrow(mu)){
   xx <- 1:ncol(mu)
   yy <- exp(mu_mature[i, ])
   col <- rainbow(nrow(mu))[i]
   ix <- which(abs(yy - exp(r$mu_instar_mature[i])) > 0.05)
   points(xx[ix], yy[ix], pch = 21, bg = col, cex = 0.5)
   hline(exp(r$mu_instar_mature[i]), col = col, lwd = 2)
}
mtext("Instar mean size (mm)", 2, 2.5, cex = 1.25)
mtext("Tow", 1, 2.5, cex = 1.25)
axis(4, at = exp(r$mu_instar_mature), labels = as.roman(as.numeric(names(mu_instars))))
delta <- c(0, which(diff(year(tows))>0), nrow(tows))
vline(delta, lwd = 2)
axis(3, at = delta[1:(length(delta)-1)] + diff(delta)/2, years)
box()

# Annual log-scale size differences:
clg()
dev.new(width = 8.5, height = 11)
m <- kronecker(matrix(1:7, ncol = 1), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, 0, cbind( 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:nrow(mu)){
   if (i == 4) mtext("Log-scale size anomaly", 2, 2.5, cex = 1.25, at = par("usr")[3])
   z <- mu[nrow(mu)-i+1, ] - r$mu_instar[nrow(mu)-i+1]
   boxplot(z ~ year(tows$date), cex = 0.2, xlab = "", ylab = "", xaxt = "n", ylim = c(-.08, 0.08))
   mtext(as.roman(nrow(mu)-i+1+3), 4, 1, cex = 1.25)
}
axis(1, at = seq(3, 33, by = 5), labels = years[seq(3, 33, by = 5)])
mtext("Year", 1, 2.75, cex = 1.5)

# Pairwise plots:
clg()
dev.new(width = 8.5, height = 11)
m <- kronecker(matrix(1:36, ncol = 6), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:6){
   for (j in 1:6){
      plot(mu[i, ], mu[j, ], cex = 0.1)
   }
}

