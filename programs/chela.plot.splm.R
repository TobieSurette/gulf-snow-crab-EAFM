chela.plot.splm <- function(data, parameters, scale = 1){

   # Calculate predicted allometric means:
   mu_imm <- splm(log(data$cw), parameters$log_alpha_immature, parameters$log_beta_immature, parameters$transition_immature, parameters$log_precision_immature)
   mu_mat <- parameters$log_alpha_mature + parameters$log_beta_mature * log(data$cw)    
   
   # Determine maturity and ouliers:
   logit_pp <- splm(data$cw, parameters$eta_alpha, parameters$eta_beta, parameters$eta_transition, parameters$log_eta_precision)
   p_mature <- 1 / (1 + exp(-logit_pp))
   
   # Error parameters:
   sigma <- exp(parameters$log_sigma)
   sigma_outlier <- exp(parameters$log_sigma_outlier)
   p_outlier <- 1 / (1 + exp(-parameters$alpha_outlier))  
   
   # Component densities:
   dout <- (1-p_mature) * p_outlier * dnorm(log(data$y), mu_imm, sigma + sigma_outlier, FALSE) +
               p_mature * p_outlier * dnorm(log(data$y), mu_mat, sigma + sigma_outlier, FALSE) 
   dimm <- (1-p_mature) * (1-p_outlier) * dnorm(log(data$y), mu_imm, sigma, FALSE)                 
   dmat <-     p_mature * (1-p_outlier) * dnorm(log(data$y), mu_mat, sigma, FALSE)  

   # Classify observations:
   pmat <- dmat / (dimm + dmat + dout)
   pimm <- dimm / (dimm + dmat + dout)
   pout <- dout / (dimm + dmat + dout)  
    
   data$class <- c("mature", "immature", "outlier")[apply(data.frame(pmat, pimm, pout), 1, which.max)]
   
   windows()
   m <- rbind(matrix(1, nrow = 5, ncol = 5), matrix(2, nrow = 3, ncol = 5))
   m <- rbind(0, cbind(0, m, 0), 0)
   layout(m)
   par(mar = c(0, 0, 0, 0))

   xx <- seq(40, 140, len = 1000)
   plot(data$cw, data$y, type = "n", ylab = "Chela height (mm)", xlim = c(40, 140), ylim = c(0, 40), cex.lab = 1.25, xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   col <- colorRamp(c("white", "blue"))(0.25)[1, ] / 255
   col <- rgb(col[1], col[2], col[3])
   points(data$cw[data$class == "mature"], data$y[data$class == "mature"], pch = 21, bg = col, cex = 0.25)
   col <- colorRamp(c("white", "green"))(0.25)[1, ] / 255
   col <- rgb(col[1], col[2], col[3])
   points(data$cw[data$class == "immature"], data$y[data$class == "immature"], pch = 21, bg = col, cex = 0.25)
   points(data$cw[data$class == "outlier"], data$y[data$class == "outlier"], pch = "x", col = "red", lwd = 2, cex = 1)
   lines(xx, exp(splm(log(xx), parameters$log_alpha_immature, parameters$log_beta_immature, parameters$transition_immature, parameters$log_precision_immature)), lwd = 2, col = "green")
   lines(xx, exp(parameters$log_alpha_mature + parameters$log_beta_mature * log(xx)), lwd = 2, col = "blue")  
   mtext("Chela height (mm)", 2, 2.5, cex = 1.25)
   box()
   
   # Plot maturity proportions:
   data$rcw <- round(data$cw)
   res <- aggregate(list(n_imm = data$class), list(cw = data[["rcw"]]), function(x) sum(x == "immature"))
   res$n_mat <- aggregate(list(n = data$class), list(cw = data[["rcw"]]), function(x) sum(x == "mature"))$n
   res$n_out <- aggregate(list(n = data$class), list(cw = data[["rcw"]]), function(x) sum(x == "mature"))$n
   res$p_mat <-  res$n_mat / (res$n_imm + res$n_mat)
   plot(c(40, 140), c(0, 1.1), type = "n", xaxs = "i", yaxs = "i")
   grid()
   dbarplot(res$p_mat, res$cw, col = "grey90", width = 1, add = TRUE)
   
   # Model:
   pp <- 1 / (1 + exp(-splm(xx, parameters$eta_alpha, parameters$eta_beta, parameters$eta_transition, parameters$log_eta_precision)))   
   lines(xx, pp, col = "red", lwd = 2)
   box()
   mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
   mtext("Mature Proportion", 2, 2.5, cex = 1.25)
   
   # Residuals:
   windows()
   m <- rbind(matrix(1, nrow = 3, ncol = 5), matrix(2, nrow = 3, ncol = 5), matrix(3, nrow = 3, ncol = 5))
   m <- rbind(0, cbind(0, m, 0), 0)
   layout(m)
   par(mar = c(0, 0, 0, 0))

   # Calculations:
   r_imm <- (log(data$y[data$class == "immature"]) - mu_imm[data$class == "immature"]) / sigma
   r_mat <- (log(data$y[data$class == "mature"]) - mu_mat[data$class == "mature"]) / sigma   
   res$logit_p_model <- splm(res$cw, parameters$eta_alpha, parameters$eta_beta, parameters$eta_transition, parameters$log_eta_precision)
   res$logit_p <- log(res$p_mat / (1-res$p_mat)) 
   res$residuals_p <- res$logit_p - res$logit_p_model
   
   plot(data$cw[data$class == "immature"], r_imm, pch = 21, bg = "grey", cex = 0.25, xaxt = "n", xlim = c(40, 140), xaxs = "i")
   abline(0, 0, col = "red", lwd = 2)
   plot(data$cw[data$class == "mature"], r_mat, pch = 21, bg = "grey", cex = 0.25, xaxt = "n", xlim = c(40, 140), xaxs = "i")
   abline(0, 0, col = "red", lwd = 2)   
   plot(res$cw, res$residuals_p, pch = 21, bg = "grey", cex = 1, xlim = c(40, 140), xaxs = "i")
   abline(0, 0, col = "red", lwd = 2) 
   
   windows()
   m <- rbind(matrix(1, nrow = 3, ncol = 5), matrix(2, nrow = 3, ncol = 5))
   m <- rbind(0, cbind(0, m, 0), 0)
   layout(m)
   par(mar = c(0, 0, 0, 0))
   
   boxplot(r_imm ~ data$station[data$class == "immature"]) 
   boxplot(r_mat ~ data$station[data$class == "mature"])    
   
}
