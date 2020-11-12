library(TMB)

Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("instar.cpp")
dyn.load(dynlib("instar"))

# Define input data:
data <- list(x = c(rnorm(100, 20, 2.5), 
                   rnorm(100, 30, 2.5),
                   rnorm(100, 40, 2.5),
                   rnorm(100, 50, 2.5),
                   rnorm(100, 60, 2.5),
                   rnorm(100, 70, 2.5),
                   rnorm(100, 80, 2.5),
                   rnorm(100, 90, 2.5)))

# Define initial parameters:
parameters <- list(log_mu_0 = log(20), 
                   log_sigma = log(2.5),
                   log_mu_inc = log(10),
                   log_sigma_inc = -2,
                   log_inc = rep(log(10), 7),
                   log_p = rep(0, 7))

obj <- MakeADFun(data, parameters, DLL = "instar", random = "log_inc")

obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt

rep <- sdreport(obj)
summary(rep, "random")
summary(rep, "fixed")


           