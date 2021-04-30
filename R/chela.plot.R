chela.plot <- function(data, parameters, observer, scale = 1){
    # Model plot:
    p <- parameters
    
    # Function for subsetting data:
    subset.data <- function(x, index){
       flag <- is.data.frame(x)
       x <- as.list(x)
       vars <- names(x)[unlist(lapply(x, length)) == max(unlist(lapply(x, length)))]
       for (i in 1:length(vars)) x[[vars[i]]] <- x[[vars[i]]][index]
       if (flag) x <- as.data.frame(x)
       return(x)
    }
       
    # Subset data by observer:
    if (!missing(observer)){
       if (is.character(observer)){
          observer <- toupper(observer)          
          if (length(observer) == 1) observer <- grep(observer, data$observers) else observer <- which(observer %in% data$observers)
       }
       vars <- names(data)[unlist(lapply(data, length)) == max(unlist(lapply(data, length)))]
       index <- data$observer %in% observer
       data <- subset.data(data, index)
    }
    
    # Mean maturity functions:
    mu_imm <- function(data, p){
       v <- data$cw + p$mu_observer
       if ("observer" %in% names(data)) v <- v + p$observer_effect[data$observer]
       if ("trip" %in% names(data)) v <- v + p$trip_effect[data$trip]
       v <- p$log_alpha_immature + p$log_beta_immature * log(v) 
       return(v)
    }
    mu_mat <- function(data, p){
       v <- data$cw + p$mu_observer
       if ("observer" %in% names(data)) v <- v + p$observer_effect[data$observer]
       if ("trip" %in% names(data)) v <- v + p$trip_effect[data$trip]
       v <- p$log_alpha_mature + p$log_beta_mature * log(v) 
       return(v)
    }

    # Function to calculate maturity probability:
    p_mat <- function(data, p){
       if ("observer" %in% names(data)){
          v = p$alpha_p + p$observer_effect_p[data$observer] + (p$beta_p + p$observer_effect_beta_p[data$observer]) * data$cw 
       }else{
          v = p$alpha_p + p$beta_p * data$cw 
       }
       if ("grid" %in% names(data)) v = v + p$grid_effect_p[data$grid]  
       if ("week" %in% names(data)) v = v + p$week_effect_p[data$week]  
       v <- 1 / (1 + exp(-v))
       return(v)
    }       
                          
    # Calculate data maturities:
    p_mature = p_mat(data, p)
    p_immature = 1 - p_mature
    
    # Chela observation error:                                            
    sigma = exp(p$log_sigma + p$observer_error_effect[data$observer])
      
    # Outlier component proportions:
    logit_p_outlier = p$alpha_outlier + p$observer_outlier_effect[data$observer]
    p_outlier = 1 / (1 + exp(-logit_p_outlier))    

    # Outlier standard errors:
    sigma_outlier = exp(p$log_sigma_outlier)
     
    # Gaussian mixture compoennt densities:                                            
    dout <- (1-p_mature) * p_outlier * dnorm(log(data$y), mu_imm(data, p), sigma + sigma_outlier, FALSE) +
                p_mature * p_outlier * dnorm(log(data$y), mu_mat(data, p), sigma + sigma_outlier, FALSE) 
    dimm <- (1-p_mature) * (1-p_outlier) * dnorm(log(data$y), mu_imm(data, p), sigma, FALSE)                 
    dmat <-     p_mature * (1-p_outlier) * dnorm(log(data$y), mu_mat(data, p), sigma, FALSE)  

    # Classify observations:
    pmat <- dmat / (dimm + dmat + dout)
    pimm <- dimm / (dimm + dmat + dout)
    pout <- dout / (dimm + dmat + dout)    
    data$class <- c("mature", "immature", "outlier")[apply(data.frame(pmat, pimm, pout), 1, which.max)]
  
    # Prepare bubble plot:   
    d <- as.data.frame(data[c("cw", "y", "class")])
    res <- aggregate(list(n_imm = data$class), by = data[c("cw", "y")], function(x) sum(x == "immature"))
    res$n_mat <- aggregate(list(n = data$class), by = data[c("cw", "y")], function(x) sum(x == "mature"))$n
    res$n_out <- aggregate(list(n = data$class), by = data[c("cw", "y")], function(x) sum(x == "outlier"))$n
    res$n <- aggregate(list(n = data$class), by = data[c("cw", "y")], length)$n
    res$max <- apply(res[c("n_imm", "n_mat", "n_out")], 1, which.max)
    
    # Chela plot:
    windows()
    m <- rbind(matrix(1, nrow = 5, ncol = 5), matrix(2, nrow = 3, ncol = 5))
    m <- rbind(0, cbind(0, m, 0), 0)
    layout(m)
    par(mar = c(0, 0, 0, 0))
    
    # Chela plot:
    plot(log(res$cw), log(res$y), type = "n", xlab = "", ylab = "", ylim = c(2.6, 3.7), xlim = c(4.3, 4.9), xaxs = "i", yaxs = "i", xaxt = "n")
    grid()
    points(log(res$cw[res$max == 1]), log(res$y[res$max == 1]), cex = 5 * scale * (nrow(res) / sum(res$n)) * sqrt(res$n[res$max == 1]), pch = 21, bg = "green")
    points(log(res$cw[res$max == 2]), log(res$y[res$max == 2]), cex = 5 * scale * (nrow(res) / sum(res$n)) * sqrt(res$n[res$max == 2]), pch = 21, bg = "blue")
    points(log(res$cw[res$max == 3]), log(res$y[res$max == 3]), cex = 1, pch = "x", col = "red", lwd = 2)
    mtext("Log(chela height)", 2, 2.5, cex = 1.25)
    if (length(unique(data$observer)) == 1) mtext(data$observers[unique(data$observer)], 3, 1.5, cex = 1.25)
       
    cw <- seq(80, 140, len = 1000)
    lines(log(cw), mu_mat(list(cw = cw, observer = observer), p), col = "blue", lwd = 3)
    lines(log(cw), mu_imm(list(cw = cw, observer = observer), p), col = "green", lwd = 3)
    
    s = exp(p$log_sigma + p$observer_error_effect[observer])
    lines(log(cw), mu_mat(list(cw = cw, observer = observer), p) - 1.96 * s, lty = "dashed", col = "blue", lwd = 2)
    lines(log(cw), mu_mat(list(cw = cw, observer = observer), p) + 1.96 * s, lty = "dashed", col = "blue", lwd = 2)
    lines(log(cw), mu_imm(list(cw = cw, observer = observer), p) - 1.96 * s, lty = "dashed", col = "green", lwd = 2)
    lines(log(cw), mu_imm(list(cw = cw, observer = observer), p) + 1.96 * s, lty = "dashed", col = "green", lwd = 2)  
    legend("topleft", legend = c("Mature", "Immature"), col = c("blue", "green"), lwd = 2, bg = "white", cex = 1.5)
    box()
    
    # Maturity proportions plot:
    plot(par("usr")[1:2], c(0, 1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")  
    grid()      
    res <- aggregate(list(n_imm = data$class), by = data[c("cw")], function(x) sum(x == "immature"))
    res$n_mat <- aggregate(list(n = data$class), by = data[c("cw")], function(x) sum(x == "mature"))$n
    res$n_out <- aggregate(list(n = data$class), by = data[c("cw")], function(x) sum(x == "outlier"))$n
    res$n <- aggregate(list(n = data$class), by = data[c("cw")], length)$n
    res$max <- apply(res[c("n_imm", "n_mat", "n_out")], 1, which.max)
    res$p <- res$n_mat / (res$n_imm + res$n_mat)
    points(log(res$cw), res$p, cex = 1 * sqrt(res$n_imm + res$n_mat), pch = 21, bg = "red")
    points(log(jitter(data$cw, amount = 0.5)), p_mat(data, p), pch = 21, bg = "grey", cex = 0.5)
    lines(log(cw), p_mat(list(cw = cw, observer = observer), p), lwd = 3)
    mtext("Log(carapace width)", 1, 2.5, cex = 1.25) 
    mtext("Mature proportion", 2, 2.5, cex = 1.25)
    box()
  
}
