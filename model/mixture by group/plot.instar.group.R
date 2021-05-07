plot.instar.group <- function(obj, data, groups, labels, xlim = c(2.0, 5), ylim = c(0, 200)){
   clg()
   dev.new(width = 8.5, height = 11)
   m <- kronecker(matrix(1:10, ncol = 2), matrix(1, ncol = 5, nrow = 5))
   m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))
   m <- s <- p <- NULL
   
   if (missing(groups)) groups <- sort(unique(data$groups))
   if (missing(labels)) labels <- groups

   for (i in 1:length(groups)){   
      ix <- which(data$group == groups[i])
      if (length(ix) > 1){
         # Plot output:
         
         z <- rep(exp(data$x[ix]), data$f[ix])
         z <- z + rep(data$precision[ix], data$f[ix]) * (runif(length(z))-0.5)
         tab <- table(round(log(z), 2))
         gbarplot(tab, 
                  border = "grey50", width = 0.01, xlim = xlim, 
                  xaxs = "i", xaxt = "n", yaxt = "n", ylim = ylim, lwd = 0.5)
         grid()
         x0 <- seq(0, 5, by = 0.01)
         d <- rep(0, length(x0))
         for (j in 1:length(obj$report()$mu_instar)){
            lines(x0, 0.01 * obj$report()$p[j,groups[i]+1] * sum(data$f[ix]) * dnorm(x0, obj$report()$mu[j,groups[i]+1], obj$report()$sigma[j,groups[i]+1]), lwd = 0.5, lty = "dashed", col = "blue")
            d <- d + .01 * obj$report()$p[j,groups[i]+1] * sum(data$f[ix]) * dnorm(x0, obj$report()$mu[j,groups[i]+1], obj$report()$sigma[j,groups[i]+1]) 
         }
         lines(x0, d, col = "blue", lwd = 0.5)
   
         mu <- obj$report()$mu[,groups[i]+1]
         ix <- match(round(mu, 2), round(x0,2))
         vline(mu, 
               col = "blue", 
               label = paste0(round(exp(mu),1), "mm"), 
               at = d[ix] + 0.05 * diff(par("usr")[3:4]),
               cex = 0.8,
               lwd = 0.5, lty = "dashed")
          
         if ((i %% 10) <=  5) axis(2)
         if ((i %% 10) ==  5) axis(1, at = seq(2, xlim[2]-0.5, by = 0.5))
         if ((i %% 10) ==  0) axis(1, at = seq(2.5, xlim[2], by = 0.5))
         if ((i %% 10) ==  0) mtext("log(cw)", 1, 2.5, at = 2.5)
         if ((i %% 10) == 3) mtext("Frequency", 2, 2.5, cex = 1.25, at = 0)
         #if (i == 16) mtext("ln(cw)", 1, 2.5, cex = 1.25)
      
         # Group label:
         text(par("usr")[1] + 0.1 * diff(par("usr")[1:2]), 
              par("usr")[3] + 0.85 * diff(par("usr")[3:4]), labels[i], cex = 1.25)
      
         box()
      }
   }
}
