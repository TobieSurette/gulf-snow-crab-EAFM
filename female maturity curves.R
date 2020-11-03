library(gulf.data)
library(gulf.graphics)

years <- 1990:2020

#gdevice(layout = 1:5)


# Compile maturity statistics:
k <- matrix(0, nrow = length(20:90), ncol = length(years))
dimnames(k) <- list(cw = 20:90, year = years)
n <- k
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
}  

# Calculate average global proportions:
pg <- apply(k/n, 1, function(x) return(mean(x[is.finite(x)], na.rm = TRUE)))
pg[as.numeric(names(pg)) > 80] <- 1 # Ad hoc fix.
p <- k/n
r <- k/n - repvec(pg, ncol = length(years))

# Maturity anomalies figure:
clg()
gdevice("pdf", file = "results/figures/sGSL female snow crab maturity anomalies")
colorbar(seq(-0.5, 0.5, by = 0.1), col = c("blue", "white", "red"),
         caption = c("logit-scale", "deviation"), smooth = TRUE)

image(as.numeric(colnames(r)), as.numeric(rownames(r)), 
      t(r), xlab = "", ylab = "", zlim = c(-0.5, 0.5),
      breaks = seq(-0.5, 0.5, by = 0.01), 
      col = colorRampPalette(c("blue", "white", "red"))(100),
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
