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

########## Multiple graphs on one page for ggplot2 #####################
# Code is credited to cookbook-R website
# 27.	Chang, Winston. "Graphs with Ggplot2." R Graphics Cookbook. Beijing: O'Reilly, 2013.
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

############# K-means Clustering #######################################
# Data: a data.table format
# k: number of desired cluster   
# initial_method: initialisation method, either 'random' or 'k-means++'
# dist_method: dissimilarity distance method, such as "Euclidean", "Manhattan", "correlation"

## K_means function
K_means <- function(data, k, initial_method, dist_method) {
  if (class(data)[1] != "data.table"){
    data <- as.data.table(data)
  }
  iteration = FALSE
  ite = 0
  WSS = 0
  cluster <- data.frame(rep(0,nrow(data)))
  ## Create initial cluster assignments for the observations
  
  ## Choose a random initial points
  means <- initialisation(initial_method, data, k, dist_method)
  
  ## Iterate until the cluster assignments stop changing
  while (!iteration) {
    ite = ite + 1
    
    ## Calculate distance
    distance <- proxy::dist(data, means, method = dist_method, convert_similarities = FALSE)
    
    distanceDT <- as.matrix(distance)^2
    
    # Assign points to nearest cluster centroid
    minDistance <- apply(distanceDT, 1, which.min)
    uniqueCluster <- unique(minDistance)
    
    biTable <- distanceDT
    biTable[1:nrow(biTable), ] = 0
    
    for (i in 1:length(uniqueCluster)) {
      index <- uniqueCluster[i]
      tempInd <- which(minDistance == index)
      # Create binary index for each point
      biTable[tempInd, index] <- 1
      # Calculate new k value
      means[index,] <- data[tempInd, sapply(.SD, mean)]
    }
    
    cluster[,ite] = as.vector(minDistance)
    # Calculate distortion value
    WSS[ite] <- sum(distanceDT * biTable)
    if (ite > 1){
      iteration = WSS[ite]== WSS[ite - 1]
    } 
    if (ite > 20 && iteration == FALSE) {
      iteration = WSS[ite] < WSS[ite-1] && WSS[ite] == WSS[ite - 2]
    }
    # minimium distance index
    cluster_summary <- plyr::count(minDistance)
    
  }
  ## Store result to the list
  k_means_result <- list(
    summary = cluster_summary, 
    cluster = cluster, 
    means = means, 
    WSS = WSS[length(WSS)],
    iteration = ite
  )
  return(k_means_result)
}

### Initialisation methods
initialisation <- function(method, data, k, dist_method) {
  if (method == "random") {
    # Randomly recreate a matrix of k x d
    means <- replicate(ncol(data), rnorm(k))
  } else if (method == "Forgy") {
    # Randomly assigned observations as mean
    random_index = sample(1:nrow(data), k)
    means <- as.matrix(data[random_index,])
  } else if (method == "k-means++") {
    # K-means++
    means <- matrix(rep(0, ncol(data)*k), ncol = ncol(data), nrow = k)
    random_index = sample(1:nrow(data), 1)
    means[1,] = as.matrix(data[random_index,])
    for (i in 2:k){
      distanceMatrix <- proxy::dist(data, means, method = dist_method, convert_similarities = FALSE)
      distanceMatrix <- as.matrix(distanceMatrix)
      minDistance <- apply(distanceMatrix, 1, min)^2
      prob <- minDistance/sum(minDistance)
      if (any(is.na(prob))) {
        random_index <- sample(1:nrow(data), 1)
      } else {
        random_index <- sample(1:nrow(data), 1, prob = prob)
      }
      means[i,] <- as.matrix(data[random_index,])
    }
  } else {
    print("Please specify the initialisation method")
  }
  return(means)
}

# Example 1: K-mean iteration process
example1 <- function(){
  simple.data = Overlapped_AD(150,2,3,0.05,0)
  simple.data.KMN <- K_means(simple.data$dataSet, 3, "random", "Manhattan")
  dd = simple.data.KMN$cluster
  colnames(dd) <- paste0("V", c(1:ncol(dd)))
  temp <- list()
  for (i in 1: ncol(dd)){
    factor.to.plot = cbind(simple.data$dataSet, cluster = dd[,i])
    mean.vector <-  factor.to.plot %>% group_by(cluster) %>% summarise(mean(V1),mean(V2))
    setnames(mean.vector,c("cluster", "V1","V2"))
    temp[[i]] <- ggplot(factor.to.plot, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) +
      geom_point(data = mean.vector, aes(V1,V2), size = 5, shape = 24, colour = "black", fill = "black") +
      theme_classic() + scale_color_hue(guide = "none") + ggtitle(paste("Iteration ", i))
  }          
  
  print("Example 1")
  multiplot(plotlist = temp, cols = min(5,length(temp)))
}
example1()
