LoadLibraries = function(){
  library(data.table)
  library(MixSim)
  library(som)
  library(plyr) 
  library(dplyr)
  library(pdist)
  library(MASS)
  library(corpcor)
  library(mvtnorm)
  library(proxy)
  library(cluster)
  library(clusterSim)
  library(reshape2)
  library(clValid)
  library(fpc)
  library(rbenchmark)
  library(RmixmodCombi)
  library(ggplot2)
  print("The libraries have been loaded.") }

LoadLibraries()

############# Create Artificial dataset#################################
# observations: Number of Observations in the dataset
# attributes: Number of Dimensions of the datase
# clusters: Number of Clusters
# overlap: Average overlap value , should be larger than 0
# Outliers: Number of Outlying observations

Overlapped_AD <- function(observations, attributes, clusters, overlap, Outliers){
  set.seed(1112)
  ms <- MixSim(BarOmega = overlap, MaxOmega = 0.15, K = clusters, p = attributes, int = c(-10,10))
  sigmak <- alply(ms$S, 3)
  dataSet <- simdataset(observations,ms$Pi,ms$Mu,ms$S, n.out = Outliers)
  
  data <- list(dataSet = as.data.table(dataSet$X), cluster = dataSet$id, w = ms$Pi, covData = sigmak, meanData = ms$Mu)
  
  return (data)
}
# Example 4: Dendrogram for Hierarchical clustering
example4 <- function(){
  simple.data = Overlapped_AD(150,2,3,0.05,0)
  print("Example 4")
  par(mfrow = c(1,4))
  linkages = c("single", "complete", "average", "centroid")
  distance <- proxy::dist(simple.data$dataSet, method = "Manhattan")
  for (i in 1:length(linkages)){
    simple.data.HC <- hclust(distance, method = linkages[i])
    plot(simple.data.HC, main = paste(linkages[i], " linkage"))
  }
}
example4()