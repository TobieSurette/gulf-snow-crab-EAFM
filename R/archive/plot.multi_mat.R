
# Selectivity plot:
clg()
file <- paste0(sex(sex), "_trawl_selectivity.pdf")
pdf(file = file, width = 5, height = 5)
p <- parameters
t <- seq(0, 140, len = 1000)
p <- 1 / (1 + exp(-exp(p$log_selectivity_slope)  * (t - p$selectivity_x50)));
plot(t, p, type = "l", ylim = c(0, 1), yaxs = "i", 
     lwd = 2, col = "blue", xaxs = "i", xlab = "", ylab = "", xlim = xlim)
grid()
mtext("Carapace width(mm)", 1, 2.5, cex = 1.25) 
mtext("Trawl selectivity", 2, 2.5, cex = 1.25) 
box()
dev.off()

# female model and length-frequencies 2010-2020.pdf
clg()
file <- paste0(sex(sex), "_length-frequencies_", min(years), "-", max(years), ".pdf")
pdf(file = file, width = 8.5, height = 11)
p <- parameters
m <- kronecker(matrix(1:n_year, ncol = 3), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, 0, m, 0),0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:length(years)){
   ix <- data$year_imm == years[i] - min(years)
   eta <- obj$report()$eta_imm[ix]
   x <- data$x_imm[ix]
   f <- data$f_imm[ix]
   plot(c(0, xlim[2]-10), c(0, 400), type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   lines(x, f, lwd = 2, col = "grey60")
   lines(x, eta, lwd = 2, col = "red")
   
   # Total matures:
   ix <- data$year_mat == years[i] - min(years)
   lines(data$x_mat[ix], data$f_mat[ix], lwd = 2, col = "lightblue")
   lines(data$x_mat[ix], obj$report()$eta_mat[ix], lwd = 2, col = "blue")

   vline(obj$report()$mu_imm[1:6,i], lty = "dashed", col = "red")
   
   # Mature recruitment:
   ix <- data$year_rec == years[i] - min(years)
   lines(data$x_rec[ix], data$f_rec[ix], lwd = 2, col = "lightblue")
   lines(data$x_rec[ix], obj$report()$eta_rec[ix], lwd = 2, col = "blue")
   
   vline(obj$report()$mu_mat[6:7,i,1], lty = "dashed", col = "blue")
   
   # Year label:
   text(par("usr")[1] + 0.15 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.85 * diff(par("usr")[3:4]), years[i], cex = 1.25)
   box()
   
   # Axes with labels:
   if (i %in% ((n_year/3)*(1:3))) axis(1)
   if (i %in% 1:(n_year/3)) axis(2)
   if (i == 3) mtext("Density", 2, 2.5, cex = 1.4)
   if (i == round(2*n_year/3)) mtext("Carapace width(mm)", 1, 3, cex = 1.4)
}
dev.off()

# Residual plot:
clg()
file <- paste0(sex(sex), "_length-frequencies_", min(years), "-", max(years), ".pdf")
pdf(file = file, width = 8.5, height = 11)
p <- parameters
m <- kronecker(matrix(1:2, ncol = 1), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, m, 0),0,0)
layout(m)
par(mar = c(0,0,0,0))
plot(range(years), c(0, 100), type = "n", xlab = "", xaxt = "n", ylab = "")
for (i in 1:length(years)){
   # Total immatures:
   ix <- data$year_imm == years[i] - min(years)
   delta <- log(data$f_imm[ix]) - log(obj$report()$eta_imm[ix])
   
   iy <- which(delta >= 0 & is.finite(delta))
   points(rep(years[i],length(iy)), data$x_imm[ix][iy], cex = 2 * sqrt(delta[iy] / max(abs(delta[is.finite(delta)]))), pch = 21, bg = "white")
   y <- data$x_imm[ix][-iy]
   points(rep(years[i],length(y)), data$x_imm[ix][-iy], cex = 2 * sqrt(delta[-iy] / max(abs(delta[is.finite(delta)]))), pch = 21, bg = "white")
}
plot(range(years), c(0, 140), type = "n", xlab = "", ylab = "")
for (i in 1:length(years)){
   # Total matures:
   ix <- data$year_mat == years[i] - min(years)
   delta <- log(data$f_mat[ix]) - log(obj$report()$eta_mat[ix])
   iy <- which(delta >= 0 & is.finite(delta))
   points(data$x_mat[ix][iy], delta[iy], cex = 1.2 * sqrt(delta[iy] / max(abs(delta))), col = "white")
   points(data$x_mat[ix][-iy], delta[-iy], cex = 1.2 * sqrt(delta[-iy] / max(abs(delta))), col = "white")
}
dev.off()

# Length-frequency data
clg()
file <- paste0(sex(sex), "_length-frequencies_grey-scale_", min(years), "-", max(years), ".pdf")
pdf(file = file, width = 8.5, height = 11)
m <- kronecker(matrix(1:2, ncol = 1), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, m, 0),0,0)
layout(m)
par(mar = c(0,0,0,0))
image(years, seq(40, 80, by = 0.5), f_mat[, as.character(seq(40, 80, by = 0.5))], col = grey(seq(1, 0, len = 100)))
for (i in 1:length(instars)) points(years, obj$report()$mu_imm[i, ], pch = 21, bg = "red", col = "red", cex = .8)
text(par("usr")[1] + 0.88 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Mature", cex = 1.5, font = 2)
axis(4, at = apply(obj$report()$mu_mat[,,1], 1, mean), labels = as.character(as.roman(4:(n_instar+3))))
box()
hline(apply(obj$report()$mu_mat[,,1], 1, mean), col = "red", lty = "dashed")
mtext("Carapace width (mm)", 2, 3, at = par("usr")[3], cex = 1.25)
image(years, seq(5, 75, by = 0.5), f_imm[, as.character(seq(5, 75, by = 0.5))], 
      col = grey(seq(1, 0, len = 100)), yaxt = "n")
mtext("Year", 1, 3, cex = 1.25)
axis(2, at = seq(10, 70, by = 10))
text(par("usr")[1] + 0.88 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Immature", cex = 1.5, font = 2)
for (i in 1:length(instars)) points(years, obj$report()$mu_imm[i, ], pch = 21, bg = "red", col = "red", cex = .8)
hline(apply(obj$report()$mu_imm, 1, median), col = "red", lty = "dashed")
box()
axis(4, at = apply(obj$report()$mu_imm, 1, median), labels = as.character(as.roman(4:(n_instar+3))))
dev.off()

# Recruitment plot:
clg()
file <- paste0(sex(sex), "_recruitment_", min(years), "-", max(years), ".pdf")
pdf(file = file, width = 8.5, height = 11)
gbarplot(exp(p$log_n_imm_instar_0), years, grid = TRUE)
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Abundance of instar IV", 2, 2.5, cex = 1.25)
dev.off()

# Moulting probability:
clg()
dev.new(height = 8.5, width = 11)
m <- kronecker(matrix(1:2, ncol = 1), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, m, 0),0,0)
layout(m)
par(mar = c(0,0,0,0))
p_mat <- obj$report()$p_mat
gbarplot(p_mat[5, ], years[-length(years)], xlab = "", ylab =  "",  ylim = c(0, 1), xaxt = "n", grid = TRUE)
mtext("Moult-to-maturity probability", 2, 2.5, at = 0, cex = 1.25)
mtext("Instar VIII to IX", 4, 1.25, cex = 1.25)
box()
gbarplot(p_mat[6, ], years[-length(years)], xlab = "", ylab =  "", ylim = c(0, 1), grid = TRUE, yaxt = "n")
mtext("Year", 1, 3.0, cex = 1.25)
mtext("Instar IX to X", 4, 1.25, cex = 1.25)
axis(2, at = seq(0, 0.8, by = 0.2))
box()

# Population abundance plot:
clg()
dev.new(height = 8.5, width = 11)
m <- kronecker(matrix(1:2, ncol = 1), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0,cbind(0, m, 0),0,0)
layout(m)
par(mar = c(0,0,0,0))
image(as.numeric(years), 4:(n_instar+3), t(obj$report()$n_imm), col = grey(seq(1, 0, len = 100)))
text(par("usr")[1] + 0.88 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Immature", cex = 1.5, font = 2)
box()
mtext("Instar", 2, 3, at = par("usr")[3], cex = 1.25)
image(as.numeric(years), 4:(n_instar+3), t(obj$report()$n_mat), col = grey(seq(1, 0, len = 100)))
text(par("usr")[1] + 0.88 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), "Mature", cex = 1.5, font = 2)
box()
mtext("Year", 1, 3, cex = 1.25)
box()

# Growth models:
alpha <- exp(parameters$log_hiatt_intercept)
beta  <- exp(parameters$log_hiatt_slope)
plot(c(0, xlim[2]), c(0, 30), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
abline(alpha[1], beta[1], col = "red", lwd = 2)
abline(alpha[2], beta[2], col = "green", lwd = 2)
mtext("Growth-per-moult", 2, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
legend("topleft", c("Immature", "Adolescent"), lwd = 2, col = c("red", "green"), cex = 1.4)

abline(-0.227, 0.429)
#abline(-2.864, 0.246)
abline(5.099, 0.071)
# Alunno-Bruscia & Sainte-Marie (1998) 
# Immatures (from settlement to onset of physiological maturation, PM): CW i+1 = 1.429 CW i – 0.227
#Prepubescent stage ɛ (from onset of PM to next prepubescent stage): CW i+1 = 1.246 CW i – 2.864
#Prepubescent stage ɸ (from prepubescent stage ɛ to pubescent): CW i+1 = 0.828 CW i – 20.585
#Pubescent to maturity: CW i+1 = 1.071 CW i – 5.099



# Year effect
clg()
year_effect <- obj$report()$year_effect
names(year_effect) <- years
gbarplot(year_effect, xaxt = "n", grid = TRUE)
hline(1, col = "red", lwd = 2)
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Relative catchability", 2, 2.5, cex = 1.25)
axis(1)

   
# Mortality:
gbarplot(obj$report()$M_mat)
obj$report()$M_imm

# Mortality:
gbarplot(obj$report()$n_mat)


plot(years, obj$report()$mu_imm[1, ])
plot(years, obj$report()$mu_imm[2, ])
plot(years, obj$report()$mu_imm[3, ])
plot(years, obj$report()$mu_imm[4, ])
plot(years, obj$report()$mu_imm[5, ])
plot(years, obj$report()$mu_imm[6, ])
plot(years, obj$report()$mu_imm[7, ])

# Immature growth anomalies:
delta <- obj$report()$mu_imm - repvec(obj$report()$mu, ncol = n_year)
image(years, 4:(n_instar+3), t(delta), col = colorRampPalette(c("red", "white", "blue"))(100), 
      zlim = c(-10, 10), xlab = "", ylab = "")
mtext("Years", 1, 2.5, cex = 1.25)
mtext("Instars", 2, 2.5, cex = 1.25)
   
# Mature growth anomalies:
delta <- obj$report()$mu_mat - repvec(obj$report()$mu, ncol = n_year)
image(years, 4:(n_instar+3), t(delta), col = colorRampPalette(c("red", "white", "blue"))(100), 
      zlim = c(-5, 5), xlab = "", ylab = "")
mtext("Years", 1, 2.5, cex = 1.25)
mtext("Instars", 2, 2.5, cex = 1.25)

# Immature instar means:
image(years, 4:(n_instar+3), t(obj$report()$mu_imm), col = colorRampPalette(c("red", "white", "blue"))(100), 
      zlim = c(0, 70), xlab = "", ylab = "")
mtext("Years", 1, 2.5, cex = 1.25)
mtext("Instars", 2, 2.5, cex = 1.25)

# Mature instar means:
image(years, 4:(n_instar+3), t(obj$report()$mu_mat), col = colorRampPalette(c("red", "white", "blue"))(100), 
      zlim = c(0, 70), xlab = "", ylab = "")
mtext("Years", 1, 2.5, cex = 1.25)
mtext("Instars", 2, 2.5, cex = 1.25)



