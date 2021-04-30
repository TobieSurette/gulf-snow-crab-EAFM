# Calculate mean size annual trends:

library(gulf.data)
library(gulf.graphics)
library(mgcv)

years <- 1990:2020
sex <- 2
path <- "results/figures/"
output <- "pdf"
if (sex == 1){
   xlim <- c(10, 140) 
   range <- c(-2, 2)
}else{
   xlim <- c(30, 75)
   range <- c(-4, 4)
}

# Load data:
b <- read.scsbio(years, survey = "regular", sex = sex)
b  <- b[which(is.new.shell(b)), ]
b$maturity <- is.mature(b)
b <- b[!is.na(b$maturity), ]
b <- b[!is.na(b$carapace.width) & (b$carapace.width <= xlim[2]), ]
b$year <- as.factor(year(b))

# Compile maturity statistics:
s <- list(mu = NULL, sigma = NULL, n = NULL)
for (i in 1:length(years)){
   # Recalculate mean size:
   cw <- b$carapace.width[which(b$maturity & b$year == years[i])]
   s$mu[i] <- mean(cw, na.rm = TRUE)
   s$sigma[i] <- sd(cw, na.rm = TRUE)
   s$n[i] <- length(cw)
}  

# Mean size of new matures:
gdevice(output, file = paste0("results/figures/sGSL ", sex(sex), " - new-shelled mature size"))
if (sex == 1) ylim <- c(70, 110) else ylim <- c(50, 70)
plot(range(years), ylim, type = "n", yaxs = "i", xlab = "", ylab = "")
grid()
gbarplot(s$mu, years, add = TRUE, width = 1, col = "grey90", border = "grey50")  
error.bar(years, lower = s$mu - 1.96 * s$sigma / sqrt(s$n), upper = s$mu + 1.96 * s$sigma / sqrt(s$n))
mtext("Survey year", 1, 2, cex = 1.5)
mtext("Carapace width (mm)", 2, 2, cex = 1.5)
mtext(paste0("sGSL ", sex(sex), " snow crab - new-shelled mature size"), 3, 1, cex = 1.5)
box()
dev.off()
