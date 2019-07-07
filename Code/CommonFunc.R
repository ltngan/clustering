############################ FUNCTION ##################################
LoadLibraries = function(){
  library(data.table)
  library(MixSim)
  library(mclust)
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
  ## Create initial cluster assignments for the observations
  means <- initialisation(initial_method, data, k, dist_method)
  
  ## Iterate until the cluster assignments stop changing
  while (!iteration) {
    ite = ite + 1
    ## Calculate distance
    while (!iteration) {
      ite = ite + 1
      ## Calculate distance
      distance <- proxy::dist(data, means, method = dist_method, convert_similarities = FALSE)
      distanceDT <- as.matrix(distance)^2
      
      # Assign points to nearest cluster centroid
      minDistance <- apply(distanceDT, 1, which.min)
      uniqueCluster <- unique(minDistance)
      
      #Create Binary table 1 if observation is assigned to that cluster, 0 if not
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
      
      # Calculate distortion value
      WSS[ite] <- sum(distanceDT * biTable)
      if (ite > 1){
        iteration = WSS[ite]== WSS[ite - 1]
      } 
      if (ite > 20 && iteration == FALSE) {
        iteration = WSS[ite] < WSS[ite-1] && WSS[ite] == WSS[ite - 2]
      }
      
      cluster_summary <- plyr::count(minDistance)
    }
    # Assign points to nearest cluster centroid
    minDistance <- apply(distanceDT, 1, which.min)
    uniqueCluster <- unique(minDistance)
    
    #Create Binary table 1 if observation is assigned to that cluster, 0 if not
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
    
    # Calculate distortion value
    WSS[ite] <- sum(distanceDT * biTable)
    if (ite > 1){
      iteration = WSS[ite]== WSS[ite - 1]
    } 
    if (ite > 20 && iteration == FALSE) {
      iteration = WSS[ite] < WSS[ite-1] && WSS[ite] == WSS[ite - 2]
    }
    
    cluster_summary <- plyr::count(minDistance)
  }
  # Record the result
  k_means_result <- list(
    summary = cluster_summary, 
    cluster = minDistance, 
    means = means, 
    WSS = WSS[length(WSS)],
    iteration = ite
  )
  return(k_means_result)
}

## K_medoids function
K_medoids <- function(data, k, initial_method, dist_method) {
  if (class(data)[1] != "data.table"){
    data <- as.data.table(data)
  }
  iteration = FALSE
  ite = 0
  WSS = 0
  ## Create initial cluster assignments for the observations
  if (initial_method == "random") {
    print("Please select other initialisation methods (Forgy or k-means++), this is not sufficient")
    break
  }
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
      temp_cluster <- data[tempInd,]
      temp_distance <-as.matrix(proxy::dist(temp_cluster, method = dist_method, convert_similarities = FALSE))
      avg_distance <- apply(temp_distance, 1, mean)
      means[index,] <- as.matrix(temp_cluster[which.min(avg_distance)[1], ])[1,]
    }
    
    # Calculate distortion value
    WSS[ite] <- sum(distanceDT * biTable)
    if (ite > 1){
      iteration = WSS[ite]== WSS[ite - 1]
    } 
    if (ite > 20 && iteration == FALSE) {
      iteration = WSS[ite] < WSS[ite-1] && WSS[ite] == WSS[ite - 2]
    }
    
    cluster_summary <- plyr::count(minDistance)
  }
  # Record result
  k_means_result <- list(
    summary = cluster_summary, 
    cluster = minDistance, 
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

############# Expectation Maximisation for Mixture of Gaussian##########
# data: a data.table format
# k: number of desired cluster
# initial_method: initialisation method, either 'random' or 'k-means'
# dist_method: dissimilarity distance method, such as "Euclidean", "Manhattan", "correlation" 

EM_main <- function(data, k, initial_method, dist_method) {
  if (class(data)[1] != "data.table"){
    data = as.data.table(data)
  }
  if (k > 1){
    initial <- initialisation.EM(data, k, initial_method, dist_method)
    w <- initial$w
    gauss <- initial$gauss
    em_table <- EM(k, data, w, gauss, initial_method, dist_method)
    return(em_table)
  } else {
    cluster <- list(cluster = rep(1, nrow(data)))
    return (cluster)
  }
}

## Support functions for EM for MOG
EM <- function (k, data, w, gauss,initial_method, dist_method) {
  atts <- ncol(data)
  obs <- nrow(data)
  log_list <- NULL
  iteration = FALSE
  ite = 0
  
  while (!iteration) {
    # Initial parameters
    if (ite > 0 && ln == 0){
      print("Reassigned initialise parameters")
      initial <- initialisation.EM(data, k, initial_method, dist_method)
      w <- initial$w
      gauss <- initial$gauss
    }
    # E-Step
    gamma <- e.step(data, k, w, gauss)
    
    if (length(which(is.na(rowSums(gamma)))) > 0 ){
      gamma <- gamma[-which(is.na(rowSums(gamma))),.SD]
    } else if (length(which(rowSums(gamma) == 0)) > 0){
      gamma <- gamma[-which(rowSums(gamma) == 0),.SD]
    } else {
      gamma <- gamma
    }
    # M-Step
    # Update effective number of points assigned to cluster
    Nk <- colSums(gamma)
    
    # Update mixing coefficient
    w <- Nk/nrow(gamma)
    
    # Update mean
    NewMean <- gamma[, laply(.SD, function(y) {
      newMeanVec = t(colSums(as.matrix(y) * data))
      return (newMeanVec)
    })]/Nk
    
    
    # Update covariance matrix
    gamma = as.data.frame(gamma)
    
    distanceCov <- list()
    for (i in 1:k){
      newMean <- NewMean[i,]
      temp <- matrix(rep(0,atts * atts), ncol = atts)
      diff = as.matrix(data - newMean)
      diff_cross <- alply(diff,1,crossprodList)
      for (j in 1:obs) {
        dd = gamma[j,i] * diff_cross[[j]]
        temp = temp + dd
      }
      distanceCov[[i]] = 1/Nk[i] * temp
    } 
    
    ## If a Gaussian has no assigned observation, find a new initial parameter
    if (any(is.na(NewMean))) {
      temp.index <- which(is.na(NewMean[,1]))
      for (j in 1:length(temp.index)){
        newRandomGaussian <- GenerateGaus(obs, atts)
        NewMean[temp.index[j],] <- newRandomGaussian$meank
        distanceCov[[temp.index[j]]] <- newRandomGaussian$sigmak
      }
    }
    
    # Evaluate the log likelihood  
    logLikelihood = data.table(rep(0,obs))
    gauss = list()
    for (i in 1:k) {
      gauss[[i]] <- list(meank = NewMean[i,], sigmak = distanceCov[[i]])
      temp <- data[,pdf(.SD,w[i],gauss[[i]]), by = .I]
      if (any(temp == Inf)) {
        temp.inf.index <- which(temp == Inf)
        for (j in 1: length(temp.inf.index)){
          temp[temp.inf.index[j],] = 0
        }
      }
      logLikelihood <- cbind(logLikelihood, temp)
    }
    logLikelihood <- logLikelihood[, V1 := NULL]
    setnames(logLikelihood, paste0("V", c(1:k)))
    
    sum_logLikelihood <- rowSums(logLikelihood)
    if(length(which(sum_logLikelihood == 0)) > 0) {
      ln <- sum(log(sum_logLikelihood[-which(sum_logLikelihood == 0)]))
    } else {
      ln <- sum(-log(sum_logLikelihood))
    }
    log_list <- c(log_list,ln)
    ite = ite + 1
    if (ite > 1 && sum(log_list) != 0){
      LL_diff = log_list[ite - 1] - log_list[ite]
      iteration = LL_diff < 10^-8
    }
    if (ite > 1000 && ln != 0) {
      print("It's not the best solution but cannot find any better initial points")
      break
    } else if (ite > 1000 && ln == 0){
      print ("IteratiTake too long to find initial points. Break the code")
      break
    }
  }
  # Assign points to the highest probability
  maxProb <- apply(logLikelihood, 1, which.max)
  N <- plyr::count(maxProb)
  
  MOG <- list(numberOfCluster = nrow(N), clusterSummary = N, cluster = maxProb, coef = w, mean = NewMean, cov = distanceCov, logLikelihood = ln, iteration = ite)
  return(MOG)
}  

GenerateGaus <- function(observations, attributes){
  # Covariance matrix
  sigma <- matrix(rnorm(attributes*attributes), ncol = attributes)
  sigma <- crossprod(sigma)+ diag(rep(0.1, attributes))
  # Mean vector
  mu <- rnorm(attributes)
  # List of initialise covariance and mean
  data <- list(sigmak = sigma, meank = mu)
  return (data)
}

initialisation.EM <- function(data, k, initial_method, dist_method){
  if (initial_method == "k-means") {
    # Initialise mixing coefficient vector
    # Generate N random numbers, compute their sum, divide each one by the sum
    temp <- K_means(data, k, "k-means++", dist_method)
    exclude.index <- which(temp$summary$freq == 1)
    
    # Initialise Covariance matrix and mean
    Gaussian <- list()
    for (i in 1:k) {
      if (length(exclude.index) > 0 && any(i == exclude.index)){
        Gaussian[[i]] <- GenerateGaus(nrow(data), ncol(data))
      } else {
        meank <- temp$means[i,]
        sigmak <- cov(data %>% filter(temp$cluster == i))
        Gaussian[[i]] <- list(meank = meank, sigmak = sigmak)
      }
    }
    # Initial mixing coefficient
    w <- temp$summary$freq
    w <- w/sum(w)
  } else if (initial_method == "random") {
    # Initialise mixing coefficient vector
    # Generate N random numbers, compute their sum, divide each one by the sum
    w <- runif(k, 0, 1)
    w <- w/sum(w)
    
    # Initialise Covariance matrix and mean
    Gaussian <- list()
    for (i in 1:k) {
      Gaussian[[i]] <- GenerateGaus(nrow(data), ncol(data))
    }
  } else {
    print("Please specify initialisation method (random or k-means)")
    break
  }
  return(list(w = w, gauss = Gaussian))
}

e.step <- function(data,k, w, gauss) {
  gamma <- data.table(rep(0,nrow(data)))
  
  for (i in 1:k) {
    temp <- data[,pdf(.SD, w[i], gauss[[i]]), by = .I]
    gamma <- cbind(gamma, temp)
  }
  gamma <- gamma[, V1 := NULL]
  setnames(gamma, paste0("V", c(1:k)))
  
  sumPDF <- rowSums(gamma)
  gamma <- gamma[,.SD]/sumPDF
  return(gamma)
}

pdf <- function(x, wk, gauss) {
  res = wk*mvtnorm::dmvnorm(x, gauss$meank, gauss$sigmak)
  return (res)
}

crossprodList <- function(x) {
  x = as.matrix(t(x))
  d = crossprod(x)
  return(d)
}

############# Hierarchical Clustering ##################################
## data: a data.table format
## k: a number of a vector of number of desirred cluster
## dist_method: choice of dissimilarity method "Euclidean", "Manhattan", "correlation"
## linkage: choice of linkage is  "single", "complete", "average"

HC <- function (data, k, dist_method, linkage){
  cluster <- NULL
  stat.values <- list()
  pairwise_data <- proxy::dist(data, method = dist_method, convert_similarities = FALSE)
  hc_temp <- hclust(pairwise_data^2, method = linkage)
  for (i in 1:length(k)) {3
    cluster = rbind(cluster,cutree(hc_temp, k = k[i]))
    stat.values[[i]] <- cluster.stats(pairwise_data, cluster[i,]) 
  }
  return(list(hierchical = hc_temp, cluster = cluster, stat = stat.values))
}

############# Combine K-mean, Hierarchical and EM in 1 run ##############
collect.func <- function (data, kmin, kmax, dist_method){
  kmeans.list <- list()
  hc.list <- HC(data, c(kmin:kmax), dist_method, "complete")
  em.list <- list()
  
  for (k in kmin:kmax){
    kmeans.list[[k]] <- K_means(data, k, "k-means++", dist_method)
    em.list[[k]] <- EM_main(data, k, "k-means", dist_method)
  }
  return(list(KM = kmeans.list,HC = hc.list, EM = em.list))
}


