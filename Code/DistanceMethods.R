############################ Dissimilarity Distance methods #################################

### Calling all needed functions

source("CommonFunc.R")
LoadLibraries()
if (!exists(c("Data_Overlap","Data_WellSeparated", "Data_Real"))){
  source("LoadData.R")
}

############## Dissimilarity distances #############################

### K-means
distance.list <- c("Euclidean", "Manhattan", "correlation")

Distances.Kmeans <- function (data, k) {
  dist.Kmeans <- list(Euclidean = list(), Manhattan = list(), Pearson = list())
  # Euclidean
  dist.Kmeans$Euclidean <- K_means(data, k, "k-means++", "Euclidean")
  # Manhattan
  dist.Kmeans$Manhattan <- K_means(data, k, "k-means++", "Manhattan")
  # correlation
  dist.Kmeans$Pearson <- K_means(data, k, "k-means++", "correlation")
  return(dist.Kmeans)
}  

Distance.KM.Experiment <- function(Data_Overlap, Data_WellSeparated, Data_Real){
  distance.KM.DT.Overlap <- Distances.Kmeans(Data_Overlap$dataSet, 5)
  print("K-means's Sum of Square Error on Overlapping Artificial dataset")
  print(ldply(distance.KM.DT.Overlap, function (x) { return(x$WSS)}))
  
  to.plot.Euc <- cbind(Data_Overlap$dataSet, cluster = distance.KM.DT.Overlap$Euclidean$cluster )
  g1 <- ggplot(to.plot.Euc, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3)  +
    ggtitle("Overlapping Artificial Dataset with Euclidean Distance") + theme_classic() + scale_color_hue("Clusters")
  
  to.plot.Man <- cbind(Data_Overlap$dataSet, cluster = distance.KM.DT.Overlap$Manhattan$cluster )
  g2 <- ggplot(to.plot.Man, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3)  +
    ggtitle("Overlapping Artificial Dataset with Manhattan Distance") + theme_classic() + scale_color_hue("Clusters")
  
  to.plot.Pear <- cbind(Data_Overlap$dataSet, cluster = distance.KM.DT.Overlap$Pearson$cluster )
  g3 <- ggplot(to.plot.Pear, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("Overlapping Artificial Dataset with Pearson Correlation") + theme_classic() + scale_color_hue("Clusters")
  
  multiplot(g1,g2,g3, cols = 3)
  
  distance.KM.DT.WellSeparated <- Distances.Kmeans(Data_WellSeparated$dataSet, 5)
  print("K-means's Sum of Square Error on Well-separated Artificial dataset")
  print(ldply(distance.KM.DT.WellSeparated, function (x) { return(x$WSS)}))
  
  to.plot.Euc <- cbind(Data_WellSeparated$dataSet, cluster = distance.KM.DT.WellSeparated$Euclidean$cluster )
  g1 <- ggplot(to.plot.Euc, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3)  +
    ggtitle("Well-separated Artificial Dataset with Euclidean Distance") + theme_classic() + scale_color_hue("Clusters")
  
  to.plot.Man <- cbind(Data_WellSeparated$dataSet, cluster = distance.KM.DT.WellSeparated$Manhattan$cluster )
  g2 <- ggplot(to.plot.Man, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3)  +
    ggtitle("Well-separated Artificial Dataset with Manhattan Distance") + theme_classic() + scale_color_hue("Clusters")
  
  to.plot.Pear <- cbind(Data_WellSeparated$dataSet, cluster = distance.KM.DT.WellSeparated$Pearson$cluster )
  g3 <- ggplot(to.plot.Pear, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3)  +
    ggtitle("Well-separated Artificial Dataset with Pearson Correlation") + theme_classic() + scale_color_hue("Clusters")
  
  multiplot(g1,g2,g3, cols = 3)
  
  distance.KM.DT.Real <- Distances.Kmeans(Data_Real, 5)
  print("K-means's Sum of Square Error on Real dataset")
  print(ldply(distance.KM.DT.Real, function (x) { return(x$WSS)}))
  cluster <- t(ldply(distance.KM.DT.Real, function(x) {return(x$cluster)}))
  cluster <- cluster[-1,]
  colnames(cluster) <- c("Euclidean", "Manhattan", "correlation")
  data.to.plot <- cbind(id = c(1:nrow(Data_Real)), Data_Real, cluster)
  data.to.plot <- melt(data.to.plot, id.vars = c("id", "Euclidean", "Manhattan", "correlation"))
  
  plotEuc <- list()
  for (i in 1:5){
    temp.to.plot <- data.to.plot[Euclidean == i,]
    plotEuc[[i]] <- ggplot(temp.to.plot, aes(x=variable, y = value, group = id)) + 
      geom_line(colour = "royalblue") + xlab("time") + ylab("Gene expression") +
      stat_summary(fun.y = mean, aes(group = 1), geom = "line", size = 3, colour = "red")
  }
  print("Plot of K-means with Euclidean distance on Real Dataset")
  multiplot(plotlist = plotEuc, cols = 5)
  
  
  plotMan <- list()
  for (i in 1:5){
    temp.to.plot <- data.to.plot[Manhattan == i,]
    plotMan[[i]] <- ggplot(temp.to.plot, aes(x=variable, y = value, group = id)) + 
      geom_line(colour = "royalblue") + xlab("time") + ylab("Gene expression") +
      stat_summary(fun.y = mean, aes(group = 1), geom = "line", size = 3, colour = "red")
  }
  print("Plot of K-means with Manhanttan distance on Real Dataset")
  multiplot(plotlist = plotMan, cols = 5)
  
  plot.Cor <- list()
  for (i in 1:5){
    temp.to.plot <- data.to.plot[correlation == i,]
    plot.Cor[[i]] <- ggplot(temp.to.plot, aes(x=variable, y = value, group = id)) + 
      geom_line(colour = "royalblue") + xlab("time") + ylab("Gene expression") +
      stat_summary(fun.y = mean, aes(group = 1), geom = "line", size = 3, colour = "red")
  }
  print("Plot of K-means with Pearson correlation distance on Real Dataset")
  multiplot(plotlist = plot.Cor, cols = 5)
} 

### Hierarchical clustering

## Dissimilarity distance
distance.HC <- function(data, k, distance.list) {
  llply(distance.list, function(y) {
    hc.result <- HC(data, k, y, "complete")
    return(hc.result)
  })
}

Distance.HC.Experiment <- function(Data_Overlap, Data_WellSeparated, Data_Real){
  distance.list <- c("Euclidean", "Manhattan", "correlation")
  
  distance.HC.Overlap <- distance.HC(Data_Overlap$dataSet, 5, distance.list)
  names(distance.HC.Overlap) <- distance.list
  print("Hierarchical's Sum of Square Error on Overlapping artificial dataset")
  print(ldply(distance.HC.Overlap, function(x){
    return(x$stat[[1]]$within.cluster.ss)
  }))
  
  distance.HC.WellSeparated <- distance.HC(Data_WellSeparated$dataSet, 5, distance.list)
  names(distance.HC.WellSeparated) <- distance.list
  print("Hierarchical's Sum of Square Error on Well-separated artificial dataset")
  print(ldply(distance.HC.WellSeparated, function(x){
    return(x$stat[[1]]$within.cluster.ss)
  }))
  
  distance.HC.Real <- distance.HC(Data_Real, 5, distance.list)
  names(distance.HC.Real) <- distance.list
  print("Hierarchical's Sum of Square Error on Real dataset")
  print(ldply(distance.HC.Real, function(x){
    return(x$stat[[1]]$within.cluster.ss)
  }))
}

## Linkages distance
linkage.HC <- function(data, k, distance, linkage.list){
  llply(linkage.list, function(y) {
    hc.result <- HC(data, k, distance, y)
    return(hc.result)
  })
}

Linkage.HC.Experiment <- function(Data_Overlap, Data_WellSeparated, Data_Real){
  par(mfrow = c(1,3))
  linkage.list <- c("single", "complete", "average")
  linkage.HC.Overlap <- linkage.HC(Data_Overlap$dataSet, 5,"Euclidean", linkage.list)
  names(linkage.HC.Overlap) <- linkage.list
  print("Hierarchical's Sum of Square Error on Overlapping artificial dataset")
  print(ldply(linkage.HC.Overlap, function (x){
    return(x$stat[[1]]$within.cluster.ss)
  }))
  print("Dendrogram for Overlapping artificial data set")
  plot(linkage.HC.Overlap$single$hierchical, main = "Single linakge")
  tryCatch(rect.hclust(linkage.HC.Overlap, k = 5, border = "red"), error = function(e) NULL)
  plot(linkage.HC.Overlap$complete$hierchical, main = "Complete linakge")
  tryCatch(rect.hclust(linkage.HC.Overlap, k = 5, border = "red"), error = function(e) NULL)
  plot(linkage.HC.Overlap$average$hierchical, main = "Average linakge")
  tryCatch(rect.hclust(linkage.HC.Overlap, k = 5, border = "red"), error = function(e) NULL)
  
  linkage.HC.WellSeparated <- linkage.HC(Data_WellSeparated$dataSet, 5,"Euclidean", linkage.list)
  names(linkage.HC.WellSeparated) <- linkage.list
  print("Hierarchical's Sum of Square Error on Well-separated artificial dataset")
  print(ldply(linkage.HC.WellSeparated, function (x){
    return(x$stat[[1]]$within.cluster.ss)
  }))
  print("Dendrogram for Well-separated artificial data set")
  par(mfrow = c(1,3))
  plot(linkage.HC.WellSeparated$single$hierchical, main = "Single linakge")
  tryCatch(rect.hclust(linkage.HC.WellSeparated, k = 5, border = "red"), error = function(e) NULL)
  plot(linkage.HC.WellSeparated$complete$hierchical, main = "Complete linakge")
  tryCatch(rect.hclust(linkage.HC.WellSeparated, k = 5, border = "red"), error = function(e) NULL)
  plot(linkage.HC.WellSeparated$average$hierchical, main = "Average linakge")
  tryCatch(rect.hclust(linkage.HC.WellSeparated, k = 5, border = "red"), error = function(e) NULL)
  
  linkage.HC.Real <- linkage.HC(Data_Real, 5,"correlation", linkage.list)
  names(linkage.HC.Real) <- linkage.list
  print("Hierarchical's Sum of Square Error on Real dataset")
  print(ldply(linkage.HC.Real, function (x){
    return(x$stat[[1]]$within.cluster.ss)
  }))
  print("Dendrogram for Real data set")
  plot(linkage.HC.Real$single$hierchical, main = "Single linakge")
  tryCatch(rect.hclust(linkage.HC.Real, k = 5, border = "red"), error = function(e) NULL)
  plot(linkage.HC.Real$complete$hierchical, main = "Complete linakge")
  tryCatch(rect.hclust(linkage.HC.Real, k = 5, border = "red"), error = function(e) NULL)
  plot(linkage.HC.Real$average$hierchical, main = "Average linakge")
  tryCatch(rect.hclust(linkage.HC.Real, k = 5, border = "red"), error = function(e) NULL)
}

### Dissimilarity distancce
## With K-means
Distance.KM.Experiment(Data_Overlap, Data_WellSeparated, Data_Real)
## With Agglomerative Hierarchicial 
Distance.HC.Experiment(Data_Overlap, Data_WellSeparated, Data_Real)
## Linkages experiment with Agglomerative Hierarchical
Linkage.HC.Experiment(Data_Overlap, Data_WellSeparated, Data_Real)
