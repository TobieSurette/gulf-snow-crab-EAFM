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
   if (length(vars) > 0){
      for (i in 1:length(vars)){
         v <- as.numeric(fixed[which(names(fixed) == vars[i])])
         if ((length(v) == length(x[[vars[i]]])) | (length(v) == 1)){
            dim(v) <- dim(x[[vars[i]]])
            x[[vars[i]]] <- v
         }else{
            print(paste0(vars[i], " problem."))
         } 
      }
   }
   
   # Update random parameters:
   vars <- unique(names(random))
   if (length(vars) > 0){
      for (i in 1:length(vars)){
         v <- as.numeric(random[which(names(random) == vars[i])])
         if ((length(v) == length(x[[vars[i]]])) | (length(v) == 1)){
            dim(v) <- dim(x[[vars[i]]])
            x[[vars[i]]] <- v
         }else{
            print(paste0(vars[i], " problem."))
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
