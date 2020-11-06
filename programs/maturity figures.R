# Maturity proportion plots:

library(gulf.data)
library(gulf.graphics)

sex <- 2
type <- "empirical"
path <- "results/tables/"
file <- paste0("sGSL ", sex(sex),  " maturity proportions - ", type, ".csv")
output <- "pdf"

# Read probabilities:
p <- read.csv(paste0(path, file))
rownames(p) <- p$carapace.width
p <- p[, -1]
colnames(p) <- gsub("[a-z,A-Z]", "", colnames(p))

if (sex == 1){
   ylim <- c(10, 140) 
   range <- c(-2, 2)
}else{
   ylim <- c(30, 75)
   range <- c(-4, 4)
}

# Maturity proportion matrix figure:
clg()
file <- paste0("results/figures/sGSL ", sex(sex), " maturity proportions - ", type)
gdevice(output, file = file)
colorbar(round(seq(0, 1, by = 0.2), 1), col = c("blue", "white", "red"), caption = c("logit-scale", "deviation"), smooth = TRUE)
image(as.numeric(colnames(p)), as.numeric(rownames(p)), 
      t(p), xlab = "", ylab = "", 
      breaks = seq(0, 1, by = 0.1), 
      zlim = c(0, 1),
      col = colorRampPalette(c("blue", "white", "red"))(10),
      ylim = ylim)
mtext("Survey year", 1, 2, cex = 1.5)
mtext("Carapace width (mm)", 2, 2, cex = 1.5)
mtext(paste0("Proportion of mature ", sex(sex), " snow crab"), 3, 1, cex = 1.5)
box()
dev.off()

# Maturity anomalies figure:
clg()
file <- paste0("results/figures/sGSL ", sex(sex), " maturity anomalies - ", type)
gdevice(output, file = file)
p <- as.matrix(p)
logit.p <- log(p/(1-p))
logit.p[!is.finite(logit.p)] <- NA
r <- logit.p - repvec(apply(logit.p, 1, mean, na.rm = TRUE), ncol = ncol(logit.p))
breaks <- seq(range[1], range[2], by = 0.1)
colorbar(round(seq(range[1], range[2], by = 0.5), 1), col = c("blue", "white", "red"), caption = c("logit-scale", "deviation"), smooth = TRUE)
image(as.numeric(colnames(r)), as.numeric(rownames(r)), 
      t(r), xlab = "", ylab = "", zlim = c(-0.5, 0.5),
      breaks = seq(range[1], range[2], by = 0.1), 
      col = colorRampPalette(c("blue", "white", "red"))(length(breaks)-1),
      ylim = ylim)
mtext("Survey year", 1, 2, cex = 1.5)
mtext("Carapace width (mm)", 2, 2, cex = 1.5)
mtext(paste0("Maturation anomalies for ", sex(sex), " snow crab"), 3, 1, cex = 1.5)
box()
dev.off()

# Maturity curves overlapped:
clg()
file <- paste0("results/figures/sGSL ", sex(sex), " maturity curves - ", type)
gdevice(output, file = file)
plot(ylim, c(0, 1), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
cols <- colorRampPalette(c("grey90", "grey60"))(ncol(p))
for (i in 1:ncol(p)){
   lines(as.numeric(row.names(p)), p[, i], lwd = 1, col = cols[i])
}   
lines(as.numeric(row.names(p)), 1 / (1 + exp(-apply(logit.p, 1, mean, na.rm = TRUE))),  lwd = 2, col = "red")
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
