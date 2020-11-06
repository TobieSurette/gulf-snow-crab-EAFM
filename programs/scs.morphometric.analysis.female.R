library(gulf.data)
library(gulf.graphics)
library(mgcv)

years <- 1990:2020

# Load data:
b <- read.scsbio(years, survey = "regular", sex = 2)
b  <- b[which(is.new.shell(b)), ]
b$maturity <- is.mature(b)
b <- b[!is.na(b$maturity), ]
b <- b[!is.na(b$carapace.width) & (b$carapace.width <= 90), ]
b$year <- as.factor(year(b))

m <- gam(maturity ~ s(carapace.width, by = year), 
         family = binomial, 
         data = b[sample(1:nrow(b), 10000), ])
x0 <- seq(20, 90, len = 1000)
plot(x0, predict(m, newdata = list(carapace.width = x0)))

# Compile maturity statistics:
logit.p <- matrix(0, nrow = length(20:90), ncol = length(years))
dimnames(logit.p) <- list(cw = 20:90, year = years)
p <- logit.p
logit.p.sd <- p
s <- list(mu = NULL, sigma = NULL, n = NULL)
for (i in 1:length(years)){
   print(years[i])
   m <- gam(maturity ~ carapace.width + s(carapace.width), 
         family = binomial, 
         data = b[b$year == years[i], ])
   tmp <- predict(m, newdata = list(carapace.width = as.numeric(row.names(p))), se.fit = TRUE)
   logit.p[,i] <- tmp$fit
   logit.p.sd[,i] <- tmp$se.fit
   
   # Recalculate mean size:
   cw <- b$carapace.width[which(b$maturity & b$year == years[i])]
   s$mu[i] <- mean(cw, na.rm = TRUE)
   s$sigma[i] <- sd(cw, na.rm = TRUE)
   s$n[i] <- length(cw)
}  

# Mean size of new matures:
gdevice("pdf", file = "results/figures/sGSL SC female - new-shelled mature size")
plot(range(years), c(50, 70), type = "n", yaxs = "i", xlab = "", ylab = "")
grid()
gbarplot(s$mu, years, add = TRUE, width = 1, col = "grey90", border = "grey50")  
error.bar(years, lower = s$mu - 1.96 * s$sigma / sqrt(s$n), upper = s$mu + 1.96 * s$sigma / sqrt(s$n))
mtext("Survey year", 1, 2, cex = 1.5)
mtext("Carapace width (mm)", 2, 2, cex = 1.5)
mtext("sGSL female snow crab - new-shelled mature size", 3, 1, cex = 1.5)
box()
dev.off()

# Calculate average global proportions:
p <- 1 / (1 + exp(-logit.p))
r <- logit.p - repvec(apply(logit.p, 1, mean), ncol = length(years))
r[logit.p.sd > 3] <- NA
p[logit.p.sd > 3] <- NA

# Maturity proportion matrix figure:
clg()
gdevice("pdf", file = "results/figures/sGSL SC female maturity proportions")
colorbar(round(seq(0, 1, by = 0.2), 1), col = c("blue", "white", "red"), caption = c("logit-scale", "deviation"), smooth = TRUE)
image(as.numeric(colnames(p)), as.numeric(rownames(p)), 
      t(p), xlab = "", ylab = "", 
      breaks = seq(0, 1, by = 0.1), 
      zlim = c(0, 1),
      col = colorRampPalette(c("blue", "white", "red"))(10),
      ylim = c(30, 75))
mtext("Survey year", 1, 2, cex = 1.5)
mtext("Carapace width (mm)", 2, 2, cex = 1.5)
mtext("Proportion of matures female snow crab", 3, 1, cex = 1.5)
box()
dev.off()

# Maturity anomalies figure:
clg()
gdevice("pdf", file = "results/figures/sGSL SC female maturity anomalies")
colorbar(round(seq(-4, 4, by = 0.5), 1), col = c("blue", "white", "red"), caption = c("logit-scale", "deviation"), smooth = TRUE)
image(as.numeric(colnames(r)), as.numeric(rownames(r)), 
      t(r), xlab = "", ylab = "", zlim = c(-0.5, 0.5),
      breaks = seq(-4, 4, by = 0.1), 
      col = colorRampPalette(c("blue", "white", "red"))(80),
      ylim = c(30, 75))
mtext("Survey year", 1, 2, cex = 1.5)
mtext("Carapace width (mm)", 2, 2, cex = 1.5)
mtext("Maturation anomalies for female snow crab", 3, 1, cex = 1.5)
box()
dev.off()

# Maturity curves overlapped:
clg()
gdevice("pdf", file = "results/figures/sGSL SC female maturity curves")
p[logit.p.sd > 3] <- NA
plot(c(20, 90), c(0, 1), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
cols <- colorRampPalette(c("grey90", "grey60"))(length(years))
for (i in 1:ncol(p)){
   lines(as.numeric(row.names(p)), p[, i], lwd = 1, col = cols[i])
}   
lines(as.numeric(row.names(p)), 1 / (1 + exp(-apply(logit.p, 1, mean))),  lwd = 2, col = "red")
mtext("Maturity proportion", 2, 2, cex = 1.5)
mtext("Carapace width (mm)", 1, 2, cex = 1.5)
mtext("Maturity curve variability", 3, 1, cex = 1.5)
legend("topleft", 
       legend = c("Average", "Annual"),
       lwd = 2, 
       col = c("red", "grey60"),
       bg = "white")
box()
dev.off()
