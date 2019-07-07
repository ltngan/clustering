
############################ FUNCTION ##################################
source("CommonFunc.R")
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

WellSeparated_AD <- function(observations, attributes, clusters){
  set.seed(1112)
  dataSet <- data.table()
  meanData <- list()
  covData <- list()
  
  for (i in 1:clusters){
    
    # Covariance matrix
    sigma <- matrix(rnorm(attributes*attributes), ncol = attributes)
    sigma <- crossprod(sigma) #+ diag(rep(0.1, attributes))
    
    # Mean vector
    mu <- runif(attributes, min = -10, max = 30)
    
    
    # Generate obs number of sample
    obs_eachCluster = floor(observations/clusters)
    tempSet <- MASS::mvrnorm(obs_eachCluster,mu, sigma)
    meanData[[i]] <- colMeans(tempSet)
    covData[[i]] <- cov(tempSet)
    dataSet <- rbind(dataSet, as.data.table(tempSet))
    if (i == 1){
      clusterIndex <- rep(i, obs_eachCluster)
    } else {
      clusterIndex <- c(clusterIndex, rep(i, obs_eachCluster))
    }
  }
  data <- list(dataSet = dataSet, cluster = clusterIndex, covData = covData, meanData = meanData)
  return (data)
}

############ Plot Artificial Dataset #############################

GeneratingDataset <- function(Data_Overlap, Data_Overlap_ExtraCluster, Data_WellSeparated){
  gplot <- list()
  ## Overlapping Dataset
  to.plot.Ovr <- cbind(Data_Overlap$dataSet, cluster = Data_Overlap$cluster)
  print(paste0("Created an overlapping artificial data set having ", nrow(Data_Overlap$dataSet), " observations and ", ncol(Data_Overlap$dataSet), " attributes"))
  gplot[[1]] <- ggplot(to.plot.Ovr, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("Overlapping Artificial Dataset") + theme_classic() + scale_color_hue("Clusters")
  
  ## Overlapping Dataset with Outliers  
  to.plot.Ovr.ex <- cbind(Data_Overlap_ExtraCluster$dataSet, cluster = Data_Overlap_ExtraCluster$cluster)
  print(paste0("Created an overlapping artificial data set with extra cluster of outliers having ", nrow(Data_Overlap_ExtraCluster$dataSet), " observations and ", ncol(Data_Overlap_ExtraCluster$dataSet), " attributes"))
  gplot[[2]] <- ggplot(to.plot.Ovr.ex, aes(V1, V2, color = factor(cluster))) +  geom_point(size = 3)  +
    ggtitle("Overlapping Artificial Dataset with Outliers") + theme_classic() + scale_color_hue("Clusters")
  
  ## Well separating Dataset
  to.plot.WS <- cbind(Data_WellSeparated$dataSet, cluster = Data_WellSeparated$cluster)
  print(paste0("Created a well-separated artificial data set having ", nrow(Data_WellSeparated$dataSet), " observations and ", ncol(Data_WellSeparated$dataSet), " attributes"))
  gplot[[3]] <- ggplot(to.plot.WS, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("Well Separating Artificial Dataset ") + theme_classic() + scale_color_hue("Clusters")
  
  ## Plot dataset
  multiplot(plotlist = gplot, cols = 3)
}

### Load all the libraries
LoadLibraries()

### Generating artificial datasets
## Overlapping dataset
Data_Overlap <- Overlapped_AD(500, 2, 5, 0.05, 0)
## Overlapping dataset with Outliers  
Data_Overlap_ExtraCluster <- Overlapped_AD(500, 2, 5, 0.05, 100)
## Well separating dataset
Data_WellSeparated <- WellSeparated_AD(500, 2, 5)
## Plot Artificial dataset
GeneratingDataset(Data_Overlap, Data_Overlap_ExtraCluster, Data_WellSeparated)

### Load the real dataset
## Remember to change the working directory to where the data is
Data_Real <- fread("CW_preprocessed_microarray.csv")
Data_Real_class <- Data_Real$V1
Data_Real[, V1 := NULL]
