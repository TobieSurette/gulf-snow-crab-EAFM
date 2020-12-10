# Function to update parameter list from TMB

update.parameters <- function(x, obj, fixed, random){
   if (!missing(obj)){
      if (all(c("par", "fn") %in% names(obj))){
         rep <- sdreport(obj)
         fixed <- summary(rep, "fixed")[, 1]
         random <- summary(rep, "random")[, 1]
      }
   }

   # Update fixed parameters:
   vars <- unique(names(fixed))
   for (i in 1:length(vars)){
      d <- dim(x[[vars[i]]])
      v <- fixed[which(names(fixed) == vars[i])]
      if (length(v) == length(x[[vars[i]]])){
         x[[vars[i]]] <- as.numeric(v)
      }else{
         print(paste0(vars[i], " problem."))
         dim(x[[vars[i]]]) <- d
      } 
   }
   
   # Update random parameters:
   vars <- unique(names(random))
   if (length(vars) > 0){
      for (i in 1:length(vars)){
         d <- dim(x[[vars[i]]])
         v <- random[which(names(random) == vars[i])]
         if (length(v) == length(x[[vars[i]]])){
            print(v)
            x[[vars[i]]] <- as.numeric(v)
         }else{
            print(paste0(vars[i], " problem."))
            dim(x[[vars[i]]]) <- d
         } 
      }
   }
   
   return(x) 
}

update.map <- function(map, free.variables, fixed.variables){
   if (!missing(free.variables)){
      free.variables <- free.variables[which(free.variables %in% names(map))]
      if (length(free.variables) > 0){
         for (i in 1:length(free.variables)){
            map[[free.variables[i]]] <- factor(1:length(map[[free.variables[i]]]))
         }
      }
   }
   
   if (!missing(fixed.variables)){
      fixed.variables <- fixed.variables[which(fixed.variables %in% names(map))]
      if (length(fixed.variables) > 0){
         for (i in 1:length(fixed.variables)){
            map[[fixed.variables[i]]] <- factor(rep(NA, length(map[[free.variables[i]]])))
         }
      }
   }

   return(map)
}
