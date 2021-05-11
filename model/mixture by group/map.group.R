library(gulf.spatial)

r <- obj$report()

vars <- c("date", "tow.number")
tows <- unique(b[c(vars, "tow")])

mu <- r$mu
dimnames(mu) <- list(instar = names(mu_instars), tow = 0:(ncol(mu)-1))
p <- r$p
dimnames(p) <- list(instar = names(mu_instars), tow = 0:(ncol(mu)-1))
sigma <- r$sigma
dimnames(sigma) <- list(instar = names(mu_instars), tow = 0:(ncol(mu)-1))

# Calculate instar counts:
n <- p * repvec(tows$n, nrow = data$n_instar)
d <- 1000000 * n / repvec(tows$swept.area, nrow = data$n_instar)
#c(10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)

# Plot instar size anomalies:
years <- sort(unique(data$year))
instars <- 4:9
for (j in 1:length(years)){
   for (i in 1:length(instars)){
      scale <- 1.5 * sqrt((50 / mu_instars[as.character(instars[i])]))
      limit <- 2.5 #round(3 * (50 / mu_instars[as.character(instars[i])]), 1)
      s <- read.scsset(years[j], valid = 1, survey = "regular")   
      
      # Instars:
      s$tow <- tows$tow[match(s[vars], tows[vars])] 
      ix <- !is.na(as.character(s$tow))
      s$value <- NA
      s$value[ix] <- exp(mu[as.character(instars[i]), as.character(s$tow[ix])]) - exp(r$mu_instar)[instars[i]-3]

      file <- paste0("maps/instar size/Female immature ", years[j], " - instar ", as.roman(instars[i]), " size anomaly.pdf")
      pdf(file = file, width = 8.5, height = 8.5)
      map.new()
      #map("bathymetry")
      map("coast")
      ix <- s$value > 0
      points(lon(s)[ix], lat(s)[ix], pch = 21, cex = scale * sqrt(s$value[ix]), bg = "blue")
      points(lon(s)[!ix], lat(s)[!ix], pch = 21, cex = scale * sqrt(-s$value[!ix]), bg = "red")
      ix <- which(is.na(s$value)) 
      points(lon(s)[ix], lat(s)[ix], pch = "x", lwd = 2)
      box()
      map.axis(1:2)
      v <- round(seq(-limit, limit, len = 7),1)
      cex = scale * sqrt(abs(v))
      cex[cex == 0] <- 1
      legend("bottomleft", 
             legend = paste0(rev(v), " mm"),
             pch = c(21,21,21, 4, 21, 21, 21),
             pt.bg = c("blue", "blue", "blue", "black", "red", "red", "red"),
             pt.cex = cex, bg = "white")
     # mtext(paste0(unique(data$year), " survey, instar ", as.roman(i), ", mean size = ", round(exp(r$mu_instar)[i-3], 1), " mm"),
   #         3, 0.5, cex = 1.5)
      dev.off()
   }
}
   
# Plot instar densities:
years <- sort(unique(data$year))
instars <- 4:9
for (j in 1:length(years)){
   for (i in 1:length(instars)){
      scale <- 0.02
      limit <- 20000
      s <- read.scsset(years[j], valid = 1, survey = "regular")   
      
      # Instars:
      s$tow <- tows$tow[match(s[vars], tows[vars])] 
      ix <- !is.na(as.character(s$tow))
      s$value <- NA
      s$value[ix] <- d[as.character(instars[i]), as.character(s$tow[ix])]

      file <- paste0("maps/instar density/Female immature ", years[j], " - instar ", as.roman(instars[i]), " density.pdf")
      pdf(file = file, width = 8.5, height = 8.5)
      map.new()
      #map("bathymetry")
      map("coast")
      ix <- s$value > 0
      points(lon(s)[ix], lat(s)[ix], pch = 21, cex = scale * sqrt(s$value[ix]), bg = "red3")
      ix <- which(is.na(s$value)) 
      points(lon(s)[ix], lat(s)[ix], pch = "x", lwd = 2)
      box()
      map.axis(1:2)
      v <- round(seq(0, limit, len = 4),1)
      cex = scale * sqrt(abs(v))
      cex[cex == 0] <- 1
      legend("bottomleft", 
             legend = paste0(round(rev(v)), " #/km2"),
             pch = c(4, 21, 21, 21),
             pt.bg = c("black", "red3", "red3", "red3"),
             pt.cex = cex, bg = "white")
     # mtext(paste0(unique(data$year), " survey, instar ", as.roman(i), ", mean size = ", round(exp(r$mu_instar)[i-3], 1), " mm"),
   #         3, 0.5, cex = 1.5)
      dev.off()
   }
}

# Annual log-scale size differencesL
clg()
dev.new(width = 8.5, height = 11)
m <- kronecker(matrix(1:6, ncol = 1), matrix(1, ncol = 5, nrow = 5))
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







