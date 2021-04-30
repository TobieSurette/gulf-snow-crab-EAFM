library(gulf.data)
library(gulf.graphics)
library(mgcv)

years <- 1990:2020
sex <- 2
path <- "results/figures/"
output <- "pdf"
if (sex == 1){
   xlim <- c(10, 140) 
   range <- c(-2, 2)
}else{
   xlim <- c(30, 90)
   range <- c(-4, 4)
}
if (sex == 1) x0 <- 10:140 else x0 <- 20:90

# Compile maturity statistics:
logit.pm <- matrix(0, nrow = length(x0), ncol = length(years))
dimnames(logit.pm) <- list(cw = x0, year = years)
logit.pm.sd <- logit.pm
p <- logit.pm
s <- list(mu = NULL, sigma = NULL, n = NULL)
for (i in 1:length(years)){
   print(years[i])

   # Load data:
   b <- read.scsbio(years[i], survey = "regular", sex = sex)
   b  <- b[which(is.new.shell(b)), ]
   b$maturity <- is.mature(b)
   b <- b[!is.na(b$maturity), ]
   b <- b[!is.na(b$carapace.width) & (b$carapace.width <= xlim[2]), ]

   # Model proportions
   m <- gam(maturity ~ carapace.width + s(carapace.width), family = binomial, data = b)
   tmp <- predict(m, newdata = list(carapace.width = as.numeric(row.names(p))), se.fit = TRUE)
   logit.pm[,i] <- tmp$fit
   logit.pm.sd[,i] <- tmp$se.fit
   
   # Empirical proportions:
   tmp <- aggregate(list(k = b$maturity), by = list(cw = round(b$carapace.width)), sum, na.rm = TRUE)
   tmp$n <- aggregate(list(n = b$maturity), by = list(cw = round(b$carapace.width)), function(x) sum(!is.na(x)))$n   
   tmp$p <- tmp$k / tmp$n
   tmp <- tmp[tmp$cw %in% x0, ]
   p[as.character(tmp$cw), i] <- tmp$p
}  

# Calculate average global proportions:
pm <- 1 / (1 + exp(-logit.pm))
pm[logit.pm.sd > 3] <- NA

# Write empirical probability table:
file <- paste0("results/tables/", "sGSL female maturity proportions - empirical.csv")
cw <- t(t(as.numeric(rownames(p))))
colnames(cw) <- "carapace.width"
write.csv(cbind(cw, p), file = file, row.names = FALSE)

# Write model probability table:
file <- paste0("results/tables/", "sGSL female maturity proportions - model.csv")
cw <- t(t(as.numeric(rownames(pm))))
colnames(cw) <- "carapace.width"
write.csv(cbind(cw, pm), file = file, row.names = FALSE)


