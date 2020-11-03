library(gulf.data)
library(gulf.graphics)
library(mgcv)

years <- 1990:2020

#gdevice(layout = 1:5)


# Compile maturity statistics:
k <- matrix(0, nrow = length(20:90), ncol = length(years))
dimnames(k) <- list(cw = 20:90, year = years)
n <- k
s <- list(mu = NULL, sigma = NULL, n = NULL)
for (i in 1:length(years)){
   print(years[i])
   # Load dataset:
   b <- read.scsbio(years[i], survey = "regular", sex = 2)
   b  <- b[which(is.new.shell(b)), ]
   b$maturity <- is.mature(b)

   # Calculate maturity proportions:
   res <- aggregate(list(k = b$maturity), by = list(cw = round(b$carapace.width)), function(x) return(length(which(x == 1))))
   res$n <- aggregate(list(n = b$maturity), by = list(cw = round(b$carapace.width)), function(x) return(sum(!is.na(x))))$n
   res$p <- res$k / res$n
   res$logit.p <- log(res$p/(1-res$p))
   
   res <- res[res$cw %in% row.names(k), ]
   
   k[as.character(res$cw), i] <- res$k
   n[as.character(res$cw), i] <- res$n
   
   s$mu[i] <- mean(b$carapace.width[which(b$maturity)], na.rm = TRUE)
   s$sigma[i] <- sd(b$carapace.width[which(b$maturity)], na.rm = TRUE)
   s$n[i] <- length(which(b$maturity))
}  

gdevice("pdf", file = "results/figures/sGSL female snow crab - new-shelled mature size")

# Mean size of new matures:
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
logit <- function(x) return(log(x/(1-x)))
pg <- apply(p, 1, function(x) return(mean(x[is.finite(x)], na.rm = TRUE)))
pg[as.numeric(names(pg)) > 80] <- 1 # Ad hoc fix.
p <- k/n
r <- logit(p) - repvec(logit(pg), ncol = length(years))

# Maturity anomalies figure:
#gdevice("pdf", file = "results/figures/sGSL female snow crab maturity anomalies")
clg()
colorbar(round(seq(-4, 4, by = 0.5), 1), col = c("blue", "white", "red"),
         caption = c("logit-scale", "deviation"), smooth = TRUE)

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
gdevice("pdf", file = "results/figures/sGSL female snow crab maturity curve variability")
plot(c(20, 90), c(0, 1), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
cols <- colorRampPalette(c("grey90", "grey60"))(length(years))
for (i in 1:ncol(p)){
   pp <- p[, i]
   pp[n[, i] < 20] <- NA
   lines(as.numeric(names(pp)), pp, lwd = 1, col = cols[i])
}   
lines(as.numeric(names(pg)), pg, lwd = 2, col = "red")
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
