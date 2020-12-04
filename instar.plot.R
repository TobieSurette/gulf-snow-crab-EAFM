plot.instar <- function(obj, data, xlim = c(0, 140), ylim = c(0, 50)){
  mu <- obj$report()$mu
  sigma <- obj$report()$sigma
  n_instar <- length(mu)
  lambda_imm <- obj$report()$lambda_imm
  lambda_mat <- obj$report()$lambda_mat 
  
  m <- kronecker(matrix(1:2), matrix(1, nrow = 5, ncol = 5))
  m <- rbind(0, cbind(0, m, 0), 0)
  layout(m)
  par(mar = c(0,0,0,0))
  
  # Immature:
  gbarplot(data$f_imm, data$x_imm, width = data$dx, border = "grey70", 
           grid = TRUE, xlim = xlim, ylim = ylim, xaxs = "i", xaxt = "n")
  grid()
  x0 <- seq(0, 140, len = 1000)
  d <- rep(0, length(x0))
  for (k in 1:(n_instar-1)){
    p <- pnorm(x0+data$dx/2, mu[k], sigma[k]) - pnorm(x0-data$dx/2, mu[k], sigma[k])
    p[is.na(p)] <- 0
    lines(x0, lambda_imm[k] * p, lwd = 1, lty = "dashed", col = "blue")
    d <- d + lambda_imm[k] * p 
  }
  lines(x0, d, col = "blue", lwd = 2)
  vline(mu, col = "red", lwd = 1, lty = "dashed")
  mtext("Frequency", 2, 2.25, cex = 1.25)
  axis(2)
  axis(3, at = mu, label = as.roman(4:(4+n_instar-1)))
  text(xlim[1] + 0.9 * diff(xlim), ylim[1] + 0.8 * diff(ylim), years, cex = 1.5)   
  box()
  mtext("Immature", 4, 1.5, cex = 1.25, srt = 180)
  
  # Mature:
  gbarplot(data$f_mat, data$x_mat, width = data$dx, border = "grey70", 
           grid = TRUE, xlim = xlim, ylim = ylim, xaxs = "i")
  grid()
  x0 <- seq(0, 140, len = 1000)
  d <- rep(0, length(x0))
  for (k in 1:(n_instar-1)){
    p <- pnorm(x0+data$dx/2, mu[k+1], sigma[k+1]) - pnorm(x0-data$dx/2, mu[k+1], sigma[k+1])
    lines(x0, lambda_mat[k+1] * p, lwd = 1, lty = "dashed", col = "blue")
    p[is.na(p)] <- 0
    d <- d + lambda_mat[k+1] * p 
  }
  lines(x0, d, col = "blue", lwd = 2)
  vline(mu, col = "red", lwd = 1, lty = "dashed")
  mtext("Frequency", 2, 2.25, cex = 1.25)
  axis(2)
  text(xlim[1] + 0.9 * diff(xlim), ylim[1] + 0.8 * diff(ylim), years + 1, cex = 1.5)
  box()
  mtext("Carapace width(mm)", 1, 2.5, cex = 1.25)
  axis(1)
  mtext("New-shelled mature", 4, 1.5, cex = 1.25, srt = 180)
}