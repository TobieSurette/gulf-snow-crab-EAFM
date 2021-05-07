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

c(10.5, 14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0)

years <- sort(unique(data$year))
instars <- 4:9
for (j in 1:length(years)){
   for (i in 1:length(instars)){
      scale <- 1.5 * sqrt((50 / mu_instars[as.character(instars[i])]))
      limit <- 2.5 #round(3 * (50 / mu_instars[as.character(instars[i])]), 1)
      s <- read.scsset(years[j], valid = 1, valid = 1)   
      
      # Instars:
      s$tow <- tows$tow[match(s[vars], tows[vars])] 
      ix <- !is.na(as.character(s$tow))
      s$value <- NA
      s$value[ix] <- exp(mu[as.character(instars[i]), as.character(s$tow[ix])]) - exp(r$mu_instar)[instars[i]-3]

      file <- paste0("maps/instar size/Female immature ", years[j], " - instar ", as.roman(instars[i]), ".pdf")
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
   


