parameters.cpp <- function(x){
   # Return model parameters declared in a TMB C++ code file.
   
   y <- readLines(x)                       # Read source code.
   y <- y[grep("PARAMETER", y)]            # Find lines with parameter declarations.
   y <- gsub("//.+", "", y)                # Remove comments.
   y <- gsub(";| ", "", y)                 # Remove spacing and separation characters.
   y <- gsub("PARAMETER[_A-Z]*", "", y)    # Remove parameter identifier.
   y <- gsub("[(|)]", "", y)               # Remove parentheses.
   y <- y[y != ""]
   
   return(y)
}

data.cpp <- function(x){
   # Returns the data declared in a TMB C++ code file.
   
   y <- readLines(x)                       # Read source code.
   y <- y[grep("DATA", y)]                 # Find lines with parameter declarations.
   y <- gsub("//.+", "", y)                # Remove comments.
   y <- gsub(";| ", "", y)                 # Remove spacing and separation characters.
   y <- gsub("DATA[_A-Z]*", "", y)         # Remove parameter identifier.
   y <- gsub("[(|)]", "", y)               # Remove parentheses.
   y <- y[y != ""]
      
   return(y)
}

keywords.cpp <- function(x){
   # Extract keywords from a TMB code file.
   
   y <- readLines(x)  
   y <- gsub("//.+", "", y)                # Remove comments.
   y <- gsub("[(]|[)]|[*]|[+]|[-]|[=]|[/]|[,]|[{]|[}]", " ", y)
   y <- gsub("[;]|[[]|[]]", " ", y)
   y <- gsub("<.*>", " ", y)
   y <- gsub("[<]|[>]", " ", y)
   y <- gsub("[:]|[#]", " ", y)
   y <- gsub("PARAMETER.* $", " ", y)
   y <- gsub("DATA.* $", " ", y)
   y <- gsub("[ ]+", " ", y)  
   y <- gsub("[.]", " ", y)
   y <- gsub("[0-9]+ ", " ", y)
   y <- paste(y, collapse = " ")
   y <- strsplit(y, " ")
   y <- unique(y[[1]])
   y <- setdiff(y, c("", "template", "Type", "include", "objective_function", "operator", "int", "return", "size",  
                     "matrix", "vector", "exp", "fill", "log", "sum", "pow", "if", "for", "else", "i", "j", "k"))
                     
   return(y)
}

extract.tmb <- function(x, variable = "estimate"){
   # EXTRACT - Extracts variables from the output of an 'sdreport' object.
   
   # Parse input arguments:
   variable <- match.arg(tolower(variable), c("estimate", "sd"))
   if (variable == "estimate") var.col <- 1
   if (variable == "sd") var.col <- 2 

   names <- unique(rownames(x))
   str <- names
   str <- gsub("_effect", "", str)
   vars <- unique(unlist(strsplit(str, "_")))
 
   # Calculate number of elements in random effect:
   temp <- expand.grid(vars, vars)
   n <- matrix(NA, nrow = nrow(temp), ncol = ncol(temp))
   for (i in 1:nrow(temp)){
      for (j in 1:2){
         a <- paste(temp[i,j], "effect", sep = "_")
         if (a %in% names) n[i,j] <- sum(rownames(x) == a)
      }
   }
   n <- data.frame(variable = paste(temp[,1], temp[,2], sep = "_"), n1 = n[,1], n2 = n[,2])
   n$variable <- as.character(n$variable)
   n[temp[,1] == temp[,2], 1] <- as.character(temp[temp[,1] == temp[,2],1])
   n[temp[,1] == temp[,2], "n2"] <- 1
   n$variable <- paste(n$variable, "effect", sep = "_")
   
   # Extract information:
   v <- list()
   for (i in 1:length(str)){
      if (length(strsplit(str[i], "_")[[1]]) == 1){
         v[[i]] <- as.numeric(x[rownames(x) == names[i], var.col])
      }
      if (length(strsplit(str[i], "_")[[1]]) == 2){
         v[[i]] <- as.numeric(x[rownames(x) == names[i], var.col])
         index <- which(names[i] == n$variable)
         if (length(index) > 0){
            if (prod(as.numeric(n[index, c("n1", "n2")])) == length(v[[i]])){
               dim(v[[i]]) <- as.numeric(n[index, c("n1", "n2")])
            }
         }
      }
   }
   names(v) <- unique(rownames(x))
   return(v)
}

indent.str <- function(x, n = 3){
   # Indent a vector of character strings:
   str <- paste(rep(" ", 3), collapse = "")
   return(paste0(str, x))
}

fit.tmb <- function(x, optim = FALSE){
   # Fit a TMB model:

   opt <- nlminb(x$par, x$fn, x$gr)
   if (optim) opt <- optim(opt$par, x$fn, x$gr, control = list(maxit = 5000, trace = 3))
   x$par <- opt$par
   rep  <- sdreport(x)
   fixed <- summary(rep, "fixed")
   random <- summary(rep, "random")
   
   # Define results:
   res <- list(obj = x, fixed = fixed, random = random, output = cbind(rep$value, rep$sd))

   return(res)
}

control.parameters <- function(x, active, inactive){
   # Used to activate/disactivate parameters to be optimized in the 'map' argument of the 'MakeADFun function.
   if (!missing(active)){
      active <- as.character(active)
      active <- active[active %in% names(x)]
      active <- unique(active)
   }
   if (!missing(inactive)){
      inactive <- as.character(inactive)
      inactive <- inactive[inactive %in% names(x)]
      inactive <- unique(inactive)
   }
   
   if (!missing(active)) if (missing(inactive)) inactive <- setdiff(names(x), active)
   if (missing(inactive)) return(x)
   for (i in 1:length(inactive)) x[[inactive[i]]] <- as.factor(NA * x[[inactive[i]]])
   
   return(x[names(x) %in% inactive])
}

update.parameters <- function(x, fixed, random){
   # Update fixed parameters:
   vars <- unique(rownames(fixed))
   for (i in 1:length(vars)){
      d <- dim(x[[vars[i]]])
      x[[vars[i]]] <- as.numeric(fixed[rownames(fixed) == vars[i], 1])
      dim(x[[vars[i]]]) <- d
   }
   # Update random parameters:
   if (!missing(random)){
      vars <- unique(rownames(random))
      for (i in 1:length(vars)){
         d <- dim(x[[vars[i]]])
         x[[vars[i]]] <- as.numeric(random[rownames(random) == vars[i], 1])
         dim(x[[vars[i]]]) <- d
      }
   }
     
   return(x) 
}

parse.parameters <- function(x, data){
   # Define dimensions of single variables:
   vars <- gsub("_effect", "", names(parameters))
   index <- unlist(lapply(strsplit(vars, "_"), length)) == 1
   vars <- gsub("length", "len", vars)
   d <- unlist(lapply(x[index], length))
   names(d) <- vars[index]
   
   v <- list()
   for (i in 1:length(x)){
      vars <- gsub("_effect", "", names(x[i]))
      vars <- unlist(strsplit(vars, "_")[[1]])
      vars <- gsub("length", "len", vars) 
      v[[i]] <- x[[i]]
      if ((length(vars) > 1) & is.null(dim(x[[i]]))){
         if (all(vars %in% names(d))) if (prod(d[vars]) == length(x[[i]])) dim(v[[i]]) <- d[vars]
      }
   }
   names(v) <- names(x)
   return(v)
}

comment.str <- function(x, comment){
   # Add a comment to a character string:
   
   if (missing(comment)){
      comment <- x 
      x <- rep("", length(x))
   }
   
   comment[comment != ""] <- paste0("// ", comment[comment != ""])
   return(paste0(x, c("", " ")[2-(nchar(x) == 0)], comment))
}

strip.comments.str <- function(x){
   # Remove comments from character strings:
   
   index <- grep("\\", x, fixed = TRUE)
   if (length(index)){
      for (i in 1:length(index)){
         x[index[i]] <- strsplit(x[index[i]], "\\", fixed = TRUE)[[1]][1]
      }
   }
   
   return(x)
}

square.str <- function(x, strip.trailing = TRUE){
   # Adds spaces to a vector of character strings so that their lengths are standardized.
   
   # Strip trailing blank spaces:
   if (strip.trailing) x <- gsub(" *$", "", x) 

   # Add trailing spaces to even out string length:
   for (i in 1:length(x)){
      x[i] <- paste0(x[i], paste(rep(" ", max(nchar(x)) - nchar(x[i])), collapse = ""))
   }
   
   return(x)
}

half.cauchy.prior <- function(variable, mu = 0, scale = 1, log.like.variable = "res"){
   # Priors over random effect standard errors (sigma ~ Half-Cauchy(0,5)):

   v <- paste0(log.like.variable, " -= log(2 * ", scale, ") - log(pi) - log(", scale, "*", scale, " + (", variable, " * ", variable, "))")
   
   return(v)
}

pad.str <- function(x, n, char = " "){
   v <- paste0(paste(rep(char, n), collapse = ""), x)    
   return(v)    
}

# List of prior forms:
#
# Half-Cauchy(0, 5)
# N(mu, sigma)
# LN(mu, sigma)
# Gam(alpha, beta)
# Beta(a,b)
# Dirichlet(alpha)

# TMB object properties:
# - CPP code
# - Input data frame.
# - Data input for compilation
# - List of model variables along with:
#     - Variable names.
#     - Factor level labels for vectors.
#     - Dimension labels for matrices.
#     - Functions for redimensioning matrices and arrays from TMB output.

#v$code 
#v$data
#v$variables
#v$link
#v$family

#offset(log(distance) - log(1.75))

gp.prior <- function(x, model = "exponential", metric = "Euclidean", discrete, 
                     coordinate.variable = "X", effect.variable, 
                     calculate.distance = TRUE, data, print = TRUE, log.like.var = "res", nugget = FALSE){ 
   # 
   
   if (is.data.frame(x)){
      if ((ncol(x) == 1) & is.factor(x[,1])){
         coordinate.variable <- names(x[1])
         if (missing(effect.variable )) effect.variable <- names(x[1])
      }
   }      
   
   # Parse input variables:
   model <- match.arg(tolower(model), c("exponential", "gaussian", "matern", "rational.quadratic"))
   
   # Initialize code variable:
   v <- NULL
   
   # Variable shortcuts:
   dvar <- paste0("Distance_", coordinate.variable)   
   cvar <- coordinate.variable
   
   # Distance matrix:
   if (calculate.distance){
      v <- NULL
      v <- c(v, comment.str(paste0("DATA_MATRIX(", cvar, ");"), "Coordinate matrix."))
      
      v <- c(v, "")
      v <- c(v, comment.str("Calculate distance matrix:"))  
      v <- c(v, paste0("Type<matrix> ", dvar, "(", cvar, ".rows(), ", cvar, ".rows());"))
      
      # Calculate Euclidean distance:
      v <- c(v, paste0("for (int i = 1; i < ", cvar, ".rows(); i++){"))
      v <- c(v, paste0("   for (int j = 1; j < i; j++){"))
      v <- c(v, paste0("      ", dvar, "(i,j) = 0;")) 
      v <- c(v, paste0("      ", "for (int k = 1; k < ", cvar, ".cols(); k++){"))
      v <- c(v, paste0("         ", dvar, "(i,j) = ", dvar, "(i,j) + (", cvar, "(i,k) - ", cvar, "(j,k)) * (", cvar, "(i,k) - ", cvar, "(j,k));"))
      v <- c(v, paste0("      ", "}"))
      v <- c(v, paste0("      ", dvar, "(i,j) = sqrt(", dvar, "(i,j));")) 
      v <- c(v, paste0("      ", dvar, "(j,i) = ", dvar, "(i,j);"))     
      v <- c(v, paste0("   ", "}"))
      v <- c(v, paste0("", "}"))  
   }   
      
   # PARAMETER DECLARATIONS:
   v <- c(v, "")
   v <- c(v, paste0("PARAMETER_VECTOR(", effect.variable, "_effect);"))
   range.var <- paste0("range_", effect.variable)
   scale.var <- paste0("scale_", effect.variable)
   epsilon.var <- paste0("epsilon_", effect.variable)
   v <- c(v, paste0("PARAMETER(", range.var, ");"))
   v <- c(v, paste0("PARAMETER(", scale.var, ");"))
   if (nugget) v <- c(v, paste0("PARAMETER(", epsilon.var, ");"))
   if (model == "matern"){
      kappa.var <- paste0("kappa_", effect.variable)
      v <- c(v, paste0("PARAMETER(", kappa.var, ");"))
   }
   if (model == "rational.quadratic"){
      alpha.var <- paste0("alpha_", effect.variable)
      v <- c(v, paste0("PARAMETER(", alpha.var, ");"))      
   }
      
   # Covariance definition:
   v <- c(v, "")
   v <- c(v, comment.str("Define covariance matrix:"))
   v <- c(v, paste0("Type<matrix> Sigma_", effect.variable, "(", cvar, ".rows(), ", cvar, ".rows());"))
   v <- c(v, paste0("for (int i = 1; i < ", cvar, ".rows(); i++){"))
   v <- c(v, paste0("   for (int j = 1; j < i; j++){"))
   
   v <- c(v, paste0("      Sigma_", effect.variable, "(i,j) = ", scale.var, " * "))
   
   # Define covariance model:
   if (model == "exponential") v[length(v)] <- paste0(v[length(v)], "exp(-", dvar, "(i,j) / range_", effect.variable, ");")
   if (model == "gaussian")    v[length(v)] <- paste0(v[length(v)], "exp(-", dvar, "(i,j)*", dvar, "(i,j) / ", range.var, ");")
   if (model == "matern")      v[length(v)] <- paste0(v[length(v)], "matern(", dvar, "(i,j), ", range.var, ", ", kappa.var, ");")
   if (model == "rational.quadratic"){ 
      v[length(v)] <- paste0(v[length(v)], "(1 + (", dvar, "(i,j) * ", dvar, "(i,j)) / (2 * ", alpha.var, " * ", range.var, " * ", range.var, "))^(-", alpha.var, ");")
   }
   
   v <- c(v, paste0("   }"))
   if (nugget) v <- c(v, paste0("   Sigma_", effect.variable, "(i,j) = Sigma_", effect.variable, "(i,j) + epsilon_", effect.variable, ";"))
   v <- c(v, paste0("}"))
    
   v <- c(v, "")
   v <-c(v, paste0("Type ", log.like.var, " += density::MVNORM_t<Type>(Sigma_", effect.variable, ")(", effect.variable, "_effect);"))
     
   if (print) cat.str(v)
   
   return(v)   
} 

as.tmb <- function(x, data, link = "identity", family = "Gaussian", loglike.var = "res", response.var = "obs"){
   # 'x' : Model formula.
   # 'data' : Data frame containing the formula variables.    
   
   # Parse input arguments:
   link <- match.arg(tolower(link), c("identity", "log", "logit"))
   family <- match.arg(tolower(family), c("gaussian", "normal", "log-normal", "poisson", "negbin"))
   
   # Extract model terms:
   if ("formula" %in% class(x)){
         terms <- attributes(terms(x))$term.labels
      str <- as.character(x)
      if (length(str) == 3) response.var <- str[2] else response.var <- "obs" 
      
      # Check for offset terms:
      str <- rownames(attributes(terms(x))$factors)
      str <- str[grep("^offset", str)]
      if (length(str) > 0){
         str <- gsub("offset", "", str)
         if ((substr(str, 1, 1) == "(") & (substr(str, nchar(str), nchar(str)) == ")")) str <- substr(str, 2, nchar(str)-1)
         offset.str <- str
         # Need to extract data variables within the offset function!
         attributes(terms(as.formula(paste("~", str))))
      }else{
         offset.str <- NULL          
      }
   }
  
   # 'x' is a character vector of variable terms:
   if (is.character(x)) terms <- x
   
   # Generate code for each term:
   v <- list()
   for (i in 1:length(terms)){
      v[[i]] <- parse.variable(terms[i], loglike.var = loglike.var)      
   }
   
   # Combine code:
   v <- catenate.list(v)
   
   # Add response variable declarations:
   v$data <- c(paste0("DATA_INT(n_obs);"), v$data) 
   v$data <- c(paste0("DATA_VECTOR(", response.var, ");"), v$data)
   
   # Complete code for linear terms:
   if (!is.null(offset.str)) v$linear.term <- c(v$linear.term, offset.str)
   nabla.str <- "nabla = "
   nabla.str <- paste0(nabla.str, v$linear.term[1])
   if (length(v$linear.term) > 1){
         for (i in 2:length(v$linear.term)){
         nabla.str[i-1] <- paste0(nabla.str[i-1], " + ")
         nabla.str[i] <- paste0("        ", v$linear.term[i])      
      }
   }
   nabla.str[length(nabla.str)] <- paste0(nabla.str[length(nabla.str)], ";")
   
   # Add inverse link function transform:
   if (link == "identity") mu.str <- "mu = nabla;"
   if (link == "log")      mu.str <- "mu = exp(nabla);"
   if (link == "logit")    mu.str <- "mu = exp(nabla) / (1 + exp(nabla));"
   
   # Loop over obvservations:
   model.str <- "for (i = 0; i < n_obs; i++){"
   model.str <- c(model.str, pad.str(nabla.str, 3))
   model.str <- c(model.str, pad.str(mu.str, 3))
   
   # Observational model:
   if (family == "gaussian"){
      model.str <- c(model.str, paste0("   ", loglike.var, " -= dnorm(", response.var, "[i], mu, sigma)"))    
   }
   
   # End loop:
   model.str <- c(model.str, "}")    
   
   v$linear.term <- model.str     
   
   return(v)
}
     
print.tmb <- function(x){
   for (i in 1:length(x)){
         for (j in 1:length(x[[i]])){
            cat(paste(x[[i]][j], "\n"))       
         }
         cat("\n")    
   }        
}    

catenate.list <- function(x){
   # Combines code generated by 'parse.variable'.    
   
   # Generate list of common list field names:
   names <- unique(unlist(lapply(x, names)))
      
   # Create result variable:
   v <- vector(mode = "list", length(names))
   names(v) <- names

   # Catenate code: 
   for (i in 1:length(x)){
         vars <- intersect(names(x[[i]]), names(v))
         for (j in 1:length(vars)){
              v[[vars[j]]] <- c(v[[vars[j]]], x[[i]][[vars[j]]])
         }    
   }
   
   # Remove redundant code:
   
   return(v)
}    

cat.str <- function(x){
   for (i in 1:length(x)){
      cat(paste0(x[i], "\n"))
   }
}

parse.variable <- function(x, type = "factor", distribution = "Gaussian", sigma.prior, 
                           loglike.var = "res", as.list = TRUE, indent = 0, byrow = TRUE){
   
   # Arguments:
   # 'x' : Character vector stating the names of variables, a factor variable or an integer vector whose values will be treated 
   #       as factor labels. If the character string is of the form "a:b", then the variable is interpreted as an interaction 
   #       term
   # 'type' : Variable type (i.e. 'integer' or 'factor'). Only used when 'x' is a character vector.
   # 'distribution' : The distributional family to be used as a prior of the generated random effect.
   # 'log.like.var' : Character string stating the name of the variable which stores the negative log-likelihood value.
   # 'as.list' : A logical value stating whether to return the code as a list or a character vector.
   # 'indent' : A positive integer stating the number of spaces to indent the generated code. 
   
   #if (!is.character(x)){
   #      if (is.list(x) 
   #}
   
   # Check input arguments:
   if (!is.character(x)) stop("'x' must be a character string.")
   type <- match.arg(tolower(type), c("integer", "real", "factor"))
   distribution <- match.arg(tolower(distribution), c("gaussian"))
   
   # Determine if variable is an interaction term:
   if (length(grep(":", x)) > 0) interaction <- TRUE else interaction <- FALSE 
   
   # Modify variable name:
   if (interaction){
      x <- gsub(":", "_", x)
      variables <- unlist(strsplit(x, "_"))
   }
   
   # Initialize result variable:
   code <- list(data = NULL, 
                parameters = NULL,
                constants = NULL, 
                transforms = NULL, 
                priors = NULL, 
                linear.term = NULL)
 
   # Define variable names:
   n.str <- paste0("n_", x)
   if (interaction){
      n.str[2] <- paste0("n_", variables[1])
      n.str[3] <- paste0("n_", variables[2])       
   }    
   effect.str <- paste0(x, "_effect")
   log.sigma.str <- paste0("log_sigma_", x)
   sigma.str <- paste0("sigma_", x)
       
   # Data declaration:
   if (!interaction){
      if (type == "integer") code$data <- c(code$data, paste0("DATA_INT(", x, ");"))   
      if (type == "real")    code$data <- c(code$data, paste0("DATA_VECTOR(", x, ");"))
      if (type == "factor")  code$data <- c(code$data, paste0("DATA_FACTOR(", x, ");"))
   }
   
   # Parameter declaration:
   if (type == "factor") code$parameters <- c(code$parameters, paste0("PARAMETER_VECTOR(", effect.str, ");"))
   code$parameters <- c(code$parameters, paste0("PARAMETER(", log.sigma.str, ");"))
  
   # Constants and parameter transforms:
   if (type == "factor") code$constants <- c(code$constants, paste0("int ", n.str[1], " = ", effect.str, ".size();"))    
   if (type == "factor") code$transforms <- c(code$transforms, paste0("Type ", sigma.str, " = exp(", log.sigma.str, ");"))
 
   # Prior declaration:
   if (distribution == "gaussian") 
      code$priors <- c(code$priors, paste0(loglike.var, " -= dnorm(", effect.str, ", Type(0), ", sigma.str, ");"))
         
   # Observation term:
   # e.g. year_effect[year[i]] 
   if (!interaction){
         code$linear.term <- c(code$linear.term, paste0(effect.str, "[", x, "[i]]"))
   }else{
         # year_region_effect : year x region matrix (row x column)
         # region[i] * n_year + year[i]    # by column
         # year[i] * n_region + region[i]  # by row
         if (byrow){
              code$linear.term <- c(code$linear.term, 
                                    paste0(effect.str, "[", variables[1], "[i] * ", n.str[3], " + ", variables[2], "[i]]"))               
         }else{
              code$linear.term <- c(code$linear.term, 
                                    paste0(effect.str, "[", variables[2], "[i] * ", n.str[2], " + ", variables[1], "[i]]"))

         }    
   }
   
   # Collapse to character vector:
   if (!as.list) code = unlist(code) 
   
   return(code)
}

