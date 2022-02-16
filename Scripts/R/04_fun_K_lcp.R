library(raster)
library(gdistance)
library(yenpathy)
library(rgdal)
library(sf)
library(magrittr)

# Function for low cost paths ---------------------------------------------
gen_top_paths <- function(tr, resist,  numpath, bufdist, orig, goal){
  originCells <- raster::cellFromXY(tr, orig)
  goalCells <- raster::cellFromXY(tr, goal)
  indexOrigin <- originCells 
  indexGoal <- goalCells 
  result.list <- vector("list",numpath)
  update.res.list <- vector("list",numpath)
  for(z in 1:numpath){
    if(z == 1){
      y <- transitionMatrix(tr)
      if(isSymmetric(y)) {
        mode <- "undirected"
      }else{
        mode <- "directed"
      }
      adjacencyGraph <- igraph::graph.adjacency(y, mode=mode, weighted=TRUE)
      E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
      shortestPaths <- get.shortest.paths(adjacencyGraph,
                                          indexOrigin, indexGoal)$vpath
      result <- tr
      transitionMatrix(result) <- Matrix::Matrix(0, ncol=ncell(tr), nrow=ncell(tr))
      for(i in 1:length(shortestPaths)){
        sPVector <- shortestPaths[[i]]
        adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
        adj <- rbind(adj,cbind(adj[,2], adj[,1]))
        transitionMatrix(result)[adj] <- 1/length(shortestPaths) + transitionMatrix(result)[adj]}
      result.rast <- raster(result)
      result.buf <- raster::buffer(result.rast, width = bufdist, doEdge=TRUE)
      result.list[[z]] <- result.buf
      result.buf.2 <- raster::buffer(result.rast, width = 2 * bufdist, doEdge=TRUE)
      update <- tr
      adj <- raster::adjacent(result.buf.2, Which(!is.na(result.buf.2), cells = TRUE), directions=16)
      transitionMatrix(update)[adj] <- cellStats(raster(tr), min)
      update.res.list[[z]] <- update
    }else{
      y <- transitionMatrix(update.res.list[[z-1]])
      if(isSymmetric(y)) {
        mode <- "undirected"
      }else{
        mode <- "directed"
      }
      adjacencyGraph <- igraph::graph.adjacency(y, mode=mode, weighted=TRUE)
      E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
      shortestPaths <- get.shortest.paths(adjacencyGraph,
                                          indexOrigin, indexGoal)$vpath
      result <- tr
      transitionMatrix(result) <- Matrix::Matrix(0, ncol=ncell(tr), nrow=ncell(tr))
      for(i in 1:length(shortestPaths)){
        sPVector <- shortestPaths[[i]]
        adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
        adj <- rbind(adj,cbind(adj[,2], adj[,1]))
        transitionMatrix(result)[adj] <- 1/length(shortestPaths) + transitionMatrix(result)[adj]}
      result.rast <- raster(result)
      result.buf <- raster::buffer(result.rast, width = bufdist, doEdge=TRUE)
      result.list[[z]] <- result.buf
      result.buf.2 <- raster::buffer(result.rast, width = 2 * bufdist, doEdge=TRUE)
      update <- update.res.list[[z-1]] 
      adj <- raster::adjacent(result.buf.2, Which(!is.na(result.buf.2), cells = TRUE), directions=16)
      transitionMatrix(update)[adj] <- cellStats(raster(tr), min)
      update.res.list[[z]] <- update
    }
  }
  all.res <- list(result.list, update.res.list)
  return(all.res)
}

