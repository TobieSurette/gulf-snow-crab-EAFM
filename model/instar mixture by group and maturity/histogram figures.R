years <- sort(unique(data$year))
xlim = c(2.5, 4.25)

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

k <- 5
for (y in 1:length(years)){  
   t <- tows[year(tows$date) == years[y], ]
   groups <- t$tow[rev(order(t$n))]
   ix <- data$group %in% groups
   tmp <- lapply(data, function(x) if (length(x) == length(ix)) return(x[ix]) else return(x))

   clg()
   file <- paste0("Female histograms SCS maturity ", years[y], ".pdf")
   pdf(file = file, width = 8.5, height = 8.5)
   m <- kronecker(matrix(1:10, ncol = 2, byrow = TRUE), matrix(1, ncol = 5, nrow = 5))
   m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))

   for (i in 1:length(groups)){ 
      ix <- which(tmp$group == groups[i])

      if ((i == 1) | ((i %% k) == 0)) ylim <- c(0, 10*(floor(sum(tmp$f[ix]) / 200) + 1))
      
      # Immature portion
      if (length(ix) > 1){
         z <- rep(exp(tmp$x[ix][tmp$maturity[ix] == 0]), tmp$f[ix][tmp$maturity[ix] == 0])
         z <- z + rep(tmp$precision[ix][tmp$maturity[ix] == 0], tmp$f[ix][tmp$maturity[ix] == 0]) * (runif(length(z))-0.5)
         tab <- table(round(log(z), 2))        
         plot(xlim, ylim, type = "n",  xaxs = "i", xaxt = "n", yaxt = "n", yaxs = "i")
         grid()
         if (length(tab) > 1){
            gbarplot(tab, border = "grey50", xlim = xlim, ylim = ylim, width = 0.1, lwd = 0.5, add = TRUE)
            x0 <- seq(0, 5, by = 0.01)
            d <- rep(0, length(x0))
            for (j in 1:nrow(mu)){
               a <- 0.01 * sum(tmp$f[ix][tmp$maturity[ix] == 0]) * p[j,as.character(groups[i])] * dnorm(x0, mu[j,as.character(groups[i])], sigma[j,as.character(groups[i])])
               lines(x0, a, lwd = 1, lty = "dashed", col = "blue")
               d <- d + a
            }
            lines(x0, d, col = "blue", lwd = 0.5)
   
            v <- mu[,as.character(groups[i])]
            vline(v, col = "blue", label = paste0(round(exp(v),1), "mm"), at = d[match(round(v, 2), round(x0,2))] + 0.05 * diff(par("usr")[3:4]), cex = 0.8,
                  lwd = 0.5, lty = "dashed")
         }
         
         # Draw axes:
         if ((i %% k) <= 5 & (i != k)) axis(2)
         if ((i %% k) == 0) axis(1, at = seq(2.5, xlim[2], by = 0.5))
         if ((i %% k) == 5) axis(1, at = seq(2, xlim[2]-0.5, by = 0.5))
         if ((i %% k) == 0) mtext("log(cw)", 1, 2.5, at = 2.5)
         if ((i %% k) == 3)  mtext("Frequency", 2, 2.5, cex = 1.25, at = mean(par("usr")[3:4]))

         # Group label:
         label <- paste0(tows$date[which(tows$tow == groups[i])], ", tow #", tows$tow.number[tows$tow == groups[i]], " (n = ",  sum(tmp$f[ix][tmp$maturity[ix] == 0]), ")")
         text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), 
              par("usr")[3] + 0.9 * diff(par("usr")[3:4]), label, pos = 1, cex = 0.80)
         
         # Mature portion:
         plot(c(3.5, 4.5), ylim, type = "n",  xaxs = "i", xaxt = "n", yaxt = "n", yaxs = "i")
         grid()
         z <- rep(exp(tmp$x[ix][tmp$maturity[ix] == 1]), tmp$f[ix][tmp$maturity[ix] == 1])
         z <- z + rep(tmp$precision[ix][tmp$maturity[ix] == 1], tmp$f[ix][tmp$maturity[ix] == 1]) * (runif(length(z))-0.5)
         tab <- table(round(log(z), 2))            

         if (length(tab) > 1){
            gbarplot(tab, col = "darkolivegreen1", border = "darkolivegreen2", width = 0.1, lwd = 1, add = TRUE)
            x0 <- seq(0, 5, by = 0.01)
            d <- rep(0, length(x0))
            for (j in 1:nrow(mu)){
               a <- 0.01 * sum(tmp$f[ix][tmp$maturity[ix] == 1]) * p_mature[j,as.character(groups[i])] * dnorm(x0, mu_mature[j,as.character(groups[i])], sigma_mature[j,as.character(groups[i])])
               lines(x0, a, lwd = 1, lty = "dashed", col = "darkolivegreen")
               d <- d + a
            }
            lines(x0, d, col = "darkolivegreen", lwd = 1)
            v <- mu_mature[,as.character(groups[i])]
            vline(v, col = "darkolivegreen", label = paste0(round(exp(v),1), "mm"), at = d[match(round(v, 2), round(x0,2))] + 0.05 * diff(par("usr")[3:4]), cex = 0.8,
                  lwd = 0.5, lty = "dashed")
            
            label <- paste0(tows$date[which(tows$tow == groups[i])], ", tow #", tows$tow.number[tows$tow == groups[i]], " (n = ",  sum(tmp$f[ix][tmp$maturity[ix] == 1]), ")")
            text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), 
                 par("usr")[3] + 0.9 * diff(par("usr")[3:4]), label, pos = 1, cex = 0.80)
         }
         
         if ((i %% k) ==  0) axis(1)

         box()
      }
   }
   dev.off()
}
   