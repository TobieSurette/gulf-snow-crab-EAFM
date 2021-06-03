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
sigma_mature <- r$sigma_mature
dimnames(sigma_mature) <- list(instar = names(mu_instars), tow = tows$tow)

# Calculate instar counts:
#c(3.19, 5.12 7.65, 10.97, 15.32, 21.02, 28.48, 38.25, 50.73, 64.53, 79.79, 96.67, 115.34, 135.99)
#10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)

# Plot instar size anomalies grouped by year:
years <- sort(unique(data$year))
instars <- 4:9
for (j in 1:length(years)){
   file <- paste0("maps/Female immature instar size anomaly maps ", years[j], ".pdf")
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
   file <- paste0("maps/Female immature instar density maps ", years[j], ".pdf")
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
instars <- 5:10
for (j in 1:length(years)){
   file <- paste0("maps/Female mature instar density maps ", years[j], ".pdf")
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

# Calculate instar frequencies:
n_imm <- t(p) * repvec(tows$n_imm, ncol = length(mu_instars))
n_mat <- t(p_mature) * repvec(tows$n_mat, ncol = length(mu_instars))
n_mat[, 7] / n_imm[, 6]
plot(t(n)[, 7])


z <- log(n_imm[, 6] / n_mat[, 6])
z[!is.finite(z)]   <- NA
z <- as.numeric(z)

clg()
file <- paste0("maps/Female ratio of immature instar IX to mature instar IX.pdf")
pdf(file = file, width = 8.5, height = 8.5)
m <- kronecker(matrix(1:4, ncol = 2), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, 0, cbind( 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (j in 1:length(years)){
   map.new()
   ix <- which((z > 0) & (year(tows) == years[j]))
   points(tows$longitude[ix], tows$latitude[ix], cex = 1.0 * sqrt(z[ix]), pch = 21, bg = "black")
   ix <- which((z < 0) & (year(tows) == years[j]))
   points(tows$longitude[ix], tows$latitude[ix], cex = 1.0 * sqrt(-z[ix]), pch = 21, bg = "red")
   map("coast")
   
   v <- seq(-3, 3)
   legend("bottomleft", 
          legend = c("20:1", "6:1", "3:1", "1:1", "1:3", "1:6", "1:20"),
          pch = 21,
          pt.bg = rev(c("red", "red", "red", "black", "black", "black", "black")),
          pt.cex = sqrt(abs(v)),
          bg = "white")

   
   text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.9 * diff(par("usr")[3:4]),
        years[j], cex = 1.25)
   box()
}
str <- "Black circles show >50% proportion of immatures, while red circles show > 50% proportions of matures" 
mtext(str, 1, at = par("usr")[1])
dev.off()
