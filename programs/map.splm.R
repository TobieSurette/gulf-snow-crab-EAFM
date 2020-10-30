map.splm <- function(data, p, scale = 3){

   index <- grep("station_effect", names(p))
   
   if (length(index) > 0){
      vars <- names(p)[index]
      for (i in 1:length(vars)){
         res <- unique(as.data.frame(data[c("station", "longitude", "latitude", "sampler")]))
         res <- res[order(res$station), ]       
         if (!all(p[[vars[i]]] == 0)){
            res <- res[res$station %in% which(p[[vars[i]]] != 0), ] 
            tmp <- data.frame(station = 1:length(p[[vars[i]]]), station_effect = p[[vars[i]]])
            res$station_effect <- tmp$station_effect[match(res$station, tmp$station)]
         
            windows(width = 10, height = 10)
            gulf.map(sea = TRUE)
            index <- res$station_effect >= 0
            size <- 1.5 * sqrt(res$station_effect[index]) / max(sqrt(res$station_effect[index]))
            points(res$longitude[index], res$latitude[index], cex = scale * size, pch = c(21:24)[as.numeric(res$sampler)], bg = "grey90")

            index <- res$station_effect < 0
            size <- 1.5 * sqrt(-res$station_effect[index]) / max(sqrt(-res$station_effect[index]))
            points(res$longitude[index], res$latitude[index], cex = scale * size, pch = c(21:24)[as.numeric(res$sampler)], bg = "black")     
            title(main = vars[i])    
         }
      }
   }
}
