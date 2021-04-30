library(gulf.data)
library(gulf.graphics)

years <- 1990:2020
sex = 1
path <- "results/figures/diagnostics/"

x0 <- 10:140
p <- matrix(NA, nrow = length(x0), ncol = length(years))
dimnames(p) <- list(cw = x0, year = years)
pm <- p
pars <- NULL
for (i in 1:length(years)){
   # Read data:
   b <- read.scsbio(years[i], survey = "regular", sex = 1)
   b <- b[which(is.new.shell(b)), ]
   b <- b[!is.na(b$carapace.width), ]
   b <- b[-which(is.na(b$chela.height) & (b$carapace.width >= 30)), ]

   # Prepare variables:
   x <- b$carapace.width
   y <- b$chela.height
   z <- rep(NA, nrow(b))
   z[which(x < 30)] <- 0 # Crab smaller than 30mm are considered immature.

   # Fit morphometric model:
   theta <- fit.morphometry.scsbio(x, y, z, sex = sex, discrete = years[i] < 1998)
   
   # Compile model parameters:
   pars <- rbind(pars, theta)
   
   # Extract maturity proportions:
   pm[, i] <- morphometry.scsbio(x0, theta = theta, discrete = years[i] < 1998)$p_mature
   v <- morphometry.scsbio(x, y, theta = theta, discrete = years[i] < 1998)$p_mature_posterior
   tmp <- aggregate(list(k = v >= 0.5), by = list(cw = round(x)), sum) 
   tmp$n <- aggregate(list(n = v >= 0.5), by = list(cw = round(x)), length)$n
   tmp$p <- tmp$k / tmp$n
   tmp <- tmp[tmp$cw %in% rownames(p), ]
   p[as.character(tmp$cw), i] <- tmp$p
   
   # Output figure:
   file <- paste0("sGSL male snow crab morphometry diagnostics ", years[i])
   gdevice("pdf", file = paste0(path, file), height = 11, width = 8.5)
   plot.morphometry.scsbio(x, y, theta, xlim = c(10, 140), title = years[i], discrete = years[i] < 1998)
   dev.off()
}
rownames(pars) <- years

# Write empirical probability table:
file <- paste0("results/tables/", "sGSL male maturity proportions - empirical.csv")
cw <- t(t(as.numeric(rownames(p))))
colnames(cw) <- "carapace.width"
write.csv(cbind(cw, p), file = file, row.names = FALSE)

# Write model probability table:
file <- paste0("results/tables/", "sGSL male maturity proportions - model.csv")
cw <- t(t(as.numeric(rownames(pm))))
colnames(cw) <- "carapace.width"
write.csv(cbind(cw, pm), file = file, row.names = FALSE)

# Write parameter estimates:
file <- paste0("results/tables/", "sGSL male maturity model parameters.csv")
tmp <- t(t(years))
colnames(tmp) <- "year"
write.csv(cbind(tmp, pars), file = file, row.names = FALSE)
