pie.plot <- function(x, y, p, cex, reference = 0, col, ...){
   # PIE.PLOT - Pie plots.
   
   # Parse input arguments:
   if (length(x) != length(y)) stop("'x' and 'y' must have the same length.")
   if (!is.numeric(x) | !is.numeric(y)) stop("'x' and 'y' must be numeric.")
   if (nrow(p) != length(x)) stop("'p' must have the same number of entries as 'x' and 'y'.")
   
   # Size of pies:
   if (missing(cex)) cex <- apply(p, 1, sum)
   if (length(cex) == 1) cex <- cex * apply(p, 1, sum)
   
   # Standardize to proportions:
   p <- p / kronecker(apply(p, 1, sum), matrix(rep(1, ncol(p)), nrow = 1))

   # Parse colours:
   if (missing(col)) col <- colorRampPalette(c("black", "white"))(ncol(p))
   print(col)
   
   angle <- seq(pi / 2, 2.5*pi, len = 500)  
   circle <- cbind(sin(angle), cos(angle))
   a <- (1:length(angle)) / length(angle)
   for (i in 1:nrow(p)){
      cp <- c(0, cumsum(p[i,]))
      for (j in 1:ncol(p)){
         index <- (a >= cp[j]) & (a < cp[j+1])
         xx <- c(0, circle[index,1], 0)
         yy <- c(0, circle[index,2], 0)
         xx <- cex[i] * xx + x[i]
         yy <- cex[i] * yy + y[i]
         polygon(xx, yy, col = col[j], ...)
      }
   }
}
