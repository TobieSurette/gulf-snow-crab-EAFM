xlim = c(2.5, 4.25)
ylim = c(0, 30)

r <- obj$report()
mu <- r$mu
dimnames(mu) <- list(instar = names(mu_instars), tow = tows$tow)
p <- r$p
dimnames(p) <- list(instar = names(mu_instars), tow = tows$tow)
sigma <- r$sigma
dimnames(sigma) <- list(instar = names(mu_instars), tow = tows$tow)

for (y in 1:length(years)){  
   t <- tows[year(tows$date) == years[y], ]
   groups <- t$tow[rev(order(t$n))]
   print(groups)
   ix <- data$group %in% groups
   tmp <- lapply(data, function(x) if (length(x) == length(ix)) return(x[ix]) else return(x))

   clg()
   file <- paste0("figures/Female immature histograms ", years[y], ".pdf")
   pdf(file = file, width = 8.5, height = 8.5)
   m <- kronecker(matrix(1:10, ncol = 2), matrix(1, ncol = 5, nrow = 5))
   m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))

   for (i in 1:length(groups)){ 
      ix <- which(tmp$group == groups[i])
      if (i %% 10 == 0){
         ylim <- c(0, 10*(floor(sum(tmp$f[ix]) / 300) + 1))
      }
      if (length(ix) > 1){
         # Plot output:
         z <- rep(exp(tmp$x[ix]), tmp$f[ix])
         z <- z + rep(tmp$precision[ix], tmp$f[ix]) * (runif(length(z))-0.5)
         tab <- table(round(log(z), 2))
         gbarplot(tab, 
                  border = "grey50", width = 0.1, xlim = xlim, 
                  xaxs = "i", xaxt = "n", yaxt = "n", ylim = ylim, lwd = 0.5)
         grid()
         x0 <- seq(0, 5, by = 0.01)
         d <- rep(0, length(x0))
         for (j in 1:nrow(mu)){
            a <- 0.01 * sum(tmp$f[ix]) * p[j,as.character(groups[i])] * dnorm(x0, mu[j,as.character(groups[i])], sigma[j,as.character(groups[i])])
            lines(x0, a, lwd = 0.5, lty = "dashed", col = "blue")
            d <- d + a
         }
         lines(x0, d, col = "blue", lwd = 0.5)
   
         v <- mu[,as.character(groups[i])]
         ix <- match(round(v, 2), round(x0,2))
         vline(v, col = "blue", label = paste0(round(exp(v),1), "mm"), at = d[ix] + 0.05 * diff(par("usr")[3:4]), cex = 0.8,
               lwd = 0.5, lty = "dashed")
         
         if ((i %% 10) <= 5 & (i != 10)) axis(2)
         if ((i %% 10) ==  5) axis(1, at = seq(2, xlim[2]-0.5, by = 0.5))
         if ((i %% 10) ==  0) axis(1, at = seq(2.5, xlim[2], by = 0.5))
         if ((i %% 10) ==  0) mtext("log(cw)", 1, 2.5, at = 2.5)
         if ((i %% 10) == 3)  mtext("Frequency", 2, 2.5, cex = 1.25, at = mean(par("usr")[3:4]))

         # Group label:
         ix <- which(tows$tow == groups[i])
         label <- paste0(tows$date[ix], ", tow #", tows$tow.number[ix], " (n = ",  tows$n[ix], ")")
         text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), 
              par("usr")[3] + 0.9 * diff(par("usr")[3:4]), label, pos = 1, cex = 0.80)
         box()
      }
   }
   dev.off()
}
   