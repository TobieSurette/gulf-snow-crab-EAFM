library(gulf.spatial)

vars <- c("date", "tow.number")
tows <- unique(b[c(vars, "tow")])

# Extract mixture parameters:
r <- obj$report()
mu <- r$mu
dimnames(mu) <- list(instar = names(mu_instars), tow = tows$tow)
mu_mature <- r$mu_mature
dimnames(mu_mature) <- list(instar = names(mu_instars), tow = tows$tow)
p <- r$p
dimnames(p) <- list(instar = names(mu_instars), tow = tows$tow)
p_mature <- r$p_mature
dimnames(p_mature) <- list(instar = names(mu_instars), tow = tows$tow)
sigma <- r$sigma
dimnames(sigma) <- list(instar = names(mu_instars), tow = tows$tow)

# Calculate instar counts:
#c(3.19, 5.12 7.65, 10.97, 15.32, 21.02, 28.48, 38.25, 50.73, 64.53, 79.79, 96.67, 115.34, 135.99)
#10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)

# Plot instar size anomalies grouped by year:
years <- sort(unique(data$year))
instars <- 4:9
for (j in 1:length(years)){
   file <- paste0("Female immature instar size anomaly maps ", years[j], ".pdf")
   pdf(file = file, width = 8.5, height = 11)
   m <- kronecker(matrix(1:6, ncol = 2), matrix(1, ncol = 5, nrow = 5))
   m <- rbind(0, 0, cbind( 0, m, 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))
   for (i in 1:length(instars)){
      scale <- 1.5 * sqrt((50 / mu_instars[as.character(instars[i])]))
      limit <- 2.5
      s <- read.scsset(years[j], valid = 1, survey = "regular")   
      
      # Instars:
      s$tow <- tows$tow[match(s[vars], tows[vars])] 
      ix <- !is.na(as.character(s$tow))
      s$value <- NA
      s$value[ix] <- exp(mu[as.character(instars[i]), as.character(s$tow[ix])]) - exp(r$mu_instar)[instars[i]-3]
      
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
      if (i == 6){
         legend("bottomleft", 
             legend = paste0(rev(v), " mm"),
             pch = c(21,21,21, 4, 21, 21, 21),
             pt.bg = c("blue", "blue", "blue", "black", "red", "red", "red"),
             pt.cex = cex, bg = "white")
      }
      if (i %in% 1:3) map.axis(2)
      if (i %in% c(3,6)) map.axis(1)
      
      text(par("usr")[1] + 0.85 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.85 * diff(par("usr")[3:4]),
           as.roman(instars[i]), cex = 1.5)
   }
   dev.off()
}

# Immatures
# Plot instar densities grouped by year:
n <- p * repvec(tows$n_imm, nrow = data$n_instar)
d <- 1000000 * n / repvec(tows$swept.area, nrow = data$n_instar)
years <- sort(unique(data$year))
instars <- 4:9
for (j in 1:length(years)){
   file <- paste0("Female immature instar density maps ", years[j], ".pdf")
   pdf(file = file, width = 8.5, height = 11)
   m <- kronecker(matrix(1:6, ncol = 2), matrix(1, ncol = 5, nrow = 5))
   m <- rbind(0, 0, cbind( 0, m, 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))
   for (i in 1:length(instars)){
      scale <- 0.02
      limit <- 20000
      s <- read.scsset(years[j], valid = 1, survey = "regular")   
      
      # Instars:
      s$tow <- tows$tow[match(s[vars], tows[vars])] 
      ix <- !is.na(as.character(s$tow))
      s$value <- NA
      s$value[ix] <- d[as.character(instars[i]), as.character(s$tow[ix])]
      
      map.new()
      #map("bathymetry")
      map("coast")
      ix <- s$value > 0
      points(lon(s)[ix], lat(s)[ix], pch = 21, cex = scale * sqrt(s$value[ix]), bg = "brown1")
      ix <- which(is.na(s$value)) 
      points(lon(s)[ix], lat(s)[ix], pch = "x", cex = 0.8, lwd = 2)
      box()
      v <- round(seq(0, limit, len = 5),1)
      cex = scale * sqrt(abs(v))
      cex[cex == 0] <- 1
      if (i == 6){
         legend("bottomleft", 
             legend = paste0(round(rev(v)), " #/km2"),
             pch = rev(c(4, 21, 21, 21, 21)),
             cex = 0.8,
             pt.bg = rev(c("black", rep("brown1", 4))),
             pt.cex = rev(cex), bg = "white")
      }
      if (i %in% 1:3) map.axis(2)
      if (i %in% c(3,6)) map.axis(1)
      
      text(par("usr")[1] + 0.85 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.85 * diff(par("usr")[3:4]),
           as.roman(instars[i]), cex = 1.5)
   }
   dev.off()
}

# Matures
n <- p_mature * repvec(tows$n, nrow = data$n_instar)
d <- 1000000 * n / repvec(tows$swept.area, nrow = data$n_instar)
# Plot instar densities grouped by year:
years <- sort(unique(data$year))
instars <- 5:9
for (j in 1:length(years)){
   file <- paste0("Female mature instar density maps ", years[j], ".pdf")
   pdf(file = file, width = 8.5, height = 11)
   m <- kronecker(matrix(1:6, ncol = 2), matrix(1, ncol = 5, nrow = 5))
   m <- rbind(0, 0, cbind( 0, m, 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))
   for (i in 1:length(instars)){
      scale <- 0.02
      limit <- 20000
      s <- read.scsset(years[j], valid = 1, survey = "regular")   
      
      # Instars:
      s$tow <- tows$tow[match(s[vars], tows[vars])] 
      ix <- !is.na(as.character(s$tow))
      s$value <- NA
      s$value[ix] <- d[as.character(instars[i]), as.character(s$tow[ix])]
      
      map.new()
      #map("bathymetry")
      map("coast")
      ix <- s$value > 0
      points(lon(s)[ix], lat(s)[ix], pch = 21, cex = scale * sqrt(s$value[ix]), bg = "brown1")
      ix <- which(is.na(s$value)) 
      points(lon(s)[ix], lat(s)[ix], pch = "x", cex = 0.8, lwd = 2)
      box()
      v <- round(seq(0, limit, len = 5),1)
      cex = scale * sqrt(abs(v))
      cex[cex == 0] <- 1
      if (i == 6){
         legend("bottomleft", 
             legend = paste0(round(rev(v)), " #/km2"),
             pch = rev(c(4, 21, 21, 21, 21)),
             cex = 0.8,
             pt.bg = rev(c("black", rep("brown1", 4))),
             pt.cex = rev(cex), bg = "white")
      }
      if (i %in% 1:3) map.axis(2)
      if (i %in% c(3,6)) map.axis(1)
      
      text(par("usr")[1] + 0.85 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.85 * diff(par("usr")[3:4]),
           as.roman(instars[i]), cex = 1.5)
   }
   dev.off()
}

# Annual log-scale size differences:
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







