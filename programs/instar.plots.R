library(gulf.data)
library(gulf.graphics)

clg()

years <- 1998:2009
years <- 2009:2020
category <- "FI"

m <- kronecker(matrix(1:length(years), ncol = 2), matrix(1, ncol = 4, nrow = 4))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:length(years)){
   # Generate Instar I
   b <- read.scsbio(years[i], survey = "regular")
   b <- b[which(is.category(b, category = category)), ]

   gbarplot(table(round(log(b$carapace.width), 2)), width = 0.01, xlim = c(2, 5), lwd = 0.5, border = "grey50",
            xaxs = "i", xaxt = "n", yaxt = "n", ylim = c(0, 200))
   vline(2.7 -  0.35, lty = "dashed", col = "blue")
   text(2.7 -  0.35, par("usr")[3] + 0.5 * diff(par("usr")[3:4]), "IV")
   vline(2.7 + (0:3) * 0.33, lty = "dashed", col = "red")
   text(2.7 + (0:3) * 0.33, par("usr")[3] + 0.8 * diff(par("usr")[3:4]), c("V", "VI", "VII", "VIII"))
   vline(2.7 + 3 * 0.33 + 0.25 * 1:4, lty = "dashed", col = "green") 
   text(2.7 + 3 * 0.33 + 0.25 * 1:4, par("usr")[3] + 0.8 * diff(par("usr")[3:4]), c("IX", "X", "XI", "XII"))
   
   if ((i / length(years)) <= 0.5) axis(2)
   if (i %in% c(length(years)/2, length(years))) axis(1)
   if (i == length(years)/2) mtext("log(cw)", 1, 2.5, at = par("usr")[2])
   if (i == 3) mtext("Frequency", 2, 2, cex = 1.25, at = 0)
   text(par("usr")[1] + 0.1 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.85 * diff(par("usr")[3:4]), years[i], cex = 1.2)
   
}
