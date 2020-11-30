library(gulf.graphics)
library(gulf.stats)
# MAke sure this works!

#x <- read.scsbio(2019, sex = 1, survey = "regular")

y <- morphometric.maturity(x)

year <- 2019
xlim = c(0, 135)
ylim = c(0, 400)
m <- kronecker(matrix(1:3), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0,cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))

x <- read.scsbio(year+1, category = "MISC12", survey = "regular")
gbarplot(freq(x), border = "grey60", width = 1, ylim = ylim, xlim = xlim, xaxs = "i", grid = TRUE, xaxt = "n")
text(120, ylim[1] + 0.8*diff(ylim), year + 1, cex = 1.5)
text(mean(xlim), ylim[1] + 0.9*diff(ylim), "Immature", , cex = 1.5)
box()

x <- read.scsbio(year, category = "MISC12", survey = "regular")
gbarplot(freq(x), border = "grey60", width = 1, ylim = ylim, xlim = xlim, xaxs = "i", grid = TRUE, xaxt = "n")
text(120, ylim[1] + 0.8*diff(ylim), year, cex = 1.5)
text(mean(xlim), ylim[1] + 0.9*diff(ylim), "Immature", , cex = 1.5)
box()
mtext("Frequency", 2, 3, cex = 1.25)

x <- read.scsbio(year+1, category = "MM", survey = "regular")
gbarplot(freq(x), border = "grey60", width = 1, ylim = ylim, xlim = xlim, xaxs = "i", grid = TRUE)
text(120, ylim[1] + 0.8*diff(ylim), year + 1, cex = 1.5)
text(mean(xlim), ylim[1] + 0.9*diff(ylim), "Mature", , cex = 1.5)
box()

mtext("Carapace width (mm)", 1, 3, cex = 1.25)

