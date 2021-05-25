
j <- 0:5
first_mature_instar   <- 5

# Instar size parameters:
mu_instar_0           <- 2.35       # Instar IV size.
log_increment        <- log(0.345)  # Log-scale mean parameter associated with instar growth increments.
log_increment_delta  <- log(0.010)  # Quadratic instar size increment effect.
log_increment_mature <- log(0.20)   # Instar size moult-to-maturity growth scaling. (new)
log_delta_maturation <- log(0.11)   # Instar size selectivity maturation bias.      (new)

# Instar error parameters:
log_sigma_0           <- -2         # Instar IV error. (new)
log_sigma_increment   <- 0          # Instar error incremental effect.  (new)
log_sigma_maturation  <- -0.25      # Instar error effect associated with maturation.  (new)

# Growth equations:
mu = mu_instar_0 + j * exp(log_increment) - 0.5 * (j-1) * j * exp(log_increment_delta);
names(mu) <- as.roman(j + 4)
mu_mature <- mu
mu["IX"]        <- mu["VIII"] +                # Previous immature instar size.
                   exp(log_increment) -        # General size increment for immatures.
                   exp(log_increment_delta) -  # Decreasing increment effect for immatures.
                   exp(log_delta_maturation)   # Penalty for size-selectivity in maturation.
mu_mature["IX"] <- mu["VIII"] +                # Previous immature instar size.
                   exp(log_increment_mature) + # General size increment for matures.
                   exp(log_delta_maturation)   # Bonus for size-selectivity in maturation.   
mu_mature["X"]  <- mu["IX"] +                  # Previous immature instar size.
                   exp(log_increment_mature)   # General size increment for matures. 

# Instar error equations:
log_sigma              <- log_sigma_0 + j * log_sigma_increment;
names(log_sigma)       <- as.roman(j + 4)
log_sigma["IX"]        <- log_sigma_mature["VIII"] + log_sigma_maturation
log_sigma_mature       <- log_sigma
log_sigma_mature["IX"] <- log_sigma["IX"]
log_sigma_mature["X"]  <- log_sigma["IX"] 

cat("Immature:\n")
print(round(exp(mu), 1))
print(round(exp(log_sigma),2))

cat("Mature:\n")
print(round(exp(mu_mature),1))
print(round(exp(log_sigma_mature),2))
