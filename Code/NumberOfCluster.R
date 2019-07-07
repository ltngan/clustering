############################ Number of Clusters methods #################################

source("CommonFunc.R")
LoadLibraries()
############ Initialisation ######################################
if (!exists(c("Data_Overlap","Data_WellSeparated", "Data_Real"))){
  source("LoadData.R")
}

### Number of Cluster

# 1. Elbow Method

elbow.func <- function (elbow.kmeans, elbow.hc, elbow.em, kmin, kmax){
  elbow.kmeans <- ldply(elbow.kmeans, function(x) {return(x$WSS)})
  elbow.kmeans <- cbind(elbow.kmeans, id = c(kmin:kmax))
  elbow.hc <- laply(elbow.hc$stat, function(x) {
    return(x$within.cluster.ss)
  })
  elbow.hc <- cbind(elbow.hc, id = c(kmin:kmax))
  elbow.em <- ldply(elbow.em, function(x) {return(x$logLikelihood)})
  elbow.em <- cbind(elbow.em, id = c(kmin:kmax))
  
  g <- list(ggplot(elbow.kmeans, aes(x = id, y = V1)) + geom_point(shape = 8) + 
              geom_line(colour = "blue", size = 1) + xlab("Number of Clusters") + ylab("Sum of Square Error") + 
              ggtitle("K-means") + theme_classic(),
            ggplot(as.data.frame(elbow.hc), aes(x = id, y = elbow.hc)) + geom_point(shape = 8) + 
              geom_line(colour = "blue", size = 1) + xlab("Number of Clusters") + ylab("Sum of Square Error") + 
              ggtitle("Hierarchical") + theme_classic(),
            ggplot(elbow.em, aes(x = id, y = V1)) + geom_point(shape = 8) + 
              geom_line(colour = "blue", size = 1) + xlab("Number of Clusters") + ylab("Log Likelihood") + 
              ggtitle("EM for MOG") + theme_classic()
  )
  
  print("Elbow methods on Artificial dataset")
  multiplot(plotlist = g, cols = 3)
}

# 2. Gap Statistic

Gap.Stat.Artificial <- function(Data_Overlap, Data_WellSeparated){
  dataSet.DF <- as.data.frame(Data_Overlap$dataSet)
  dataSet2.DF <- as.data.frame(Data_WellSeparated$dataSet)
  mets <- eval(formals(maxSE)$method)
  
  # K-means
  print("Gap statistic for K-means")
  gskmn1 <- clusGap(dataSet.DF,FUN = K_means, initial_method = "k-means++", dist_method = "Euclidean", K.max = 15)
  gstab = data.frame(gskmn1$Tab, k = 1:nrow(gskmn1$Tab))
  p1 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p1 = p1 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  p1 = p1 + ggtitle("Gap statistic for overlapping dataset: K-means++")
  
  gskmn2 <- clusGap(dataSet2.DF,FUN = K_means, initial_method = "k-means++", dist_method = "Euclidean", K.max = 15)
  gstab = data.frame(gskmn2$Tab, k = 1:nrow(gskmn2$Tab))
  p2 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p2 = p2 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  p2 = p2 + ggtitle("Gap statistic for well-separated dataset: K-means++")
  
  multiplot(p1, p2, cols = 2)
  
  # Hierarchical Clustering
  print("Gap statistic on Hierarchical clustering")
  gshc1 <- clusGap(dataSet.DF,FUN = HC, dist_method = "Euclidean", linkage = "complete", K.max = 15)
  #  sapply(c(1/4, 1,2,4), function(SEf)
  #    sapply(mets, function(M) maxSE(gshc1$Tab[,"gap"], gshc1$Tab[,"SE.sim"], method = M, SE.factor = SEf)))
  
  gstab = data.frame(gshc1$Tab, k = 1:nrow(gshc1$Tab))
  p1 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p1 = p1 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  p1 = p1 + ggtitle("Gap statistic for overlapping dataset: Hierarchical")
  
  gshc2 <- clusGap(dataSet2.DF,FUN = HC, dist_method = "Euclidean", linkage = "complete", K.max = 15)
  gstab = data.frame(gshc2$Tab, k = 1:nrow(gshc2$Tab))
  p2 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p2 = p2 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  p2 = p2 + ggtitle("Gap statistic for well-separated dataset: Hierarchical")
  
  multiplot(p1, p2, cols = 2)
  
  # Expectation maximisation
  yn <- quote(var1 <- readline('Gap statistic for EM for MOG will take approximate 2 hours. \n Do you want to continue (Y/N)? '))
  eval(yn)
  if(var1 == "Y"){
    gsem1 <- clusGap(dataSet.DF, FUN = EM_main, initial_method = "k-means", dist_method = "Euclidean", K.max = 15)
    #  sapply(c(1/4, 1,2,4), function(SEf)
    #    sapply(mets, function(M) maxSE(gsem1$Tab[,"gap"], gsem1$Tab[,"SE.sim"], method = M, SE.factor = SEf)))
    
    gstab = data.frame(gsem1$Tab, k = 1:nrow(gsem1$Tab))
    p1 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
    p1 = p1 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
    p1 = p1 + ggtitle("Gap statistic for overlapping dataset: EM for MOG")
    
    gsem2 <- clusGap(dataSet2.DF,  FUN = EM_main, initial_method = "k-means", dist_method = "Euclidean", K.max = 15)
    gstab = data.frame(gsem2$Tab, k = 1:nrow(gsem2$Tab))
    p2 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
    p2 = p2 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
    p2 = p2 + ggtitle("Gap statistic for well-separated dataset: EM for MOG")
    
    multiplot(p1, p2, cols = 2)
  }
}

Gap.Stat.Real <- function(Data_Real){
  warning("This will execute for more than 2 hours")
  df.real <- as.data.frame(Data_Real)
  gskmn.real <- clusGap(df.real,FUN = K_means, initial_method = "k-means++", dist_method = "correlation", K.max = 15, B = 50)
  gstab = data.frame(gskmn.real$Tab, k = 1:nrow(gskmn.real$Tab))
  p1 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p1 = p1 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  p1 = p1 + ggtitle("Gap statistic for real dataset: K-means++")
  
  gshc.real <- clusGap(df.real,FUN = HC, dist_method = "correlation", linkage = "complete", K.max = 15, B = 50)
  gstab = data.frame(gshc.real$Tab, k = 1:nrow(gshc.real$Tab))
  p2 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  p2 = p2 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  p2 = p2 + ggtitle("Gap statistic for real dataset: Hierarchical")
  
  #gsem.real <- clusGap(df.real,FUN = EM_main, initial_method = "k-means", dist_method = "correlation", K.max = 15, B = 50)
  #gstab = data.frame(gsem.real$Tab, k = 1:nrow(gsem.real$Tab))
  #p3 = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size = 5)
  #p3 = p3 + geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), color = "red")
  #p3 = p3 + ggtitle("Gap statistic for real dataset: EM for MOG")
  print("Gap statistic for real data set")
  multiplot(p1,p2, cols = 2)
}

# 3. BIC & Combining mixture component for clustering
# K-means
BIC.KMN.Func <- function(data, KMNresult, kmin, kmax) {
  N <- nrow(data)
  d <- ncol(data)
  BIC.KMN.Pen <- c(kmin:kmax)*(d+1/2*(d+1))*N
  
  BIC.KMN <- ldply(KMNresult, function (x) {
    logLikelihood = data.table(rep(0, nrow(data)))
    gauss = list()
    tryCatch(
      {k = nrow(x$summary)
      for(i in 1:k) {
        gauss[[i]] <- list(meank = x$means[i,], sigmak = cov(data%>%filter(x$cluster == i)))
        temp <- data[,pdf(.SD,1,gauss[[i]]), by = .I]
        if (any(temp == Inf)) {
          temp.inf.index <- which(temp == Inf)
          for (j in 1: length(temp.inf.index)){
            temp[temp.inf.index[j],] = 0
          }
        }
        logLikelihood <- cbind(logLikelihood, temp)
      }
      logLikelihood <- logLikelihood[, V1 := NULL]
      setnames(logLikelihood, paste0("V", c(1:nrow(x$summary))))
      sum.log <- sum(-1/2*x$summary$freq*colSums(logLikelihood))
      return(sum.log )}, error = function(e) NULL)
  })
  
  BIC.KMN <- BIC.KMN - BIC.KMN.Pen
  BIC.KMN.BestK <- c(kmin:kmax)[which.max(BIC.KMN[,"V1"])]
  print("BIC for K-means return the optimal  number of cluster is: ")
  print(BIC.KMN.BestK)
  return(BIC.KMN.BestK)
}

# EM for MOG
BIC.EM.func <- function(data, kmin, kmax, dist_method){
  BIC_EM <- Mclust(data, G = c(kmin:kmax), initialization = hclust(proxy::dist(data, method = dist_method, convert_similarities = FALSE)))
  print(paste("BIC for EM returns the optimal number of components is ", BIC_EM$G))
  return(BIC_EM)
}

BIC.CMC.Artificial <- function(Data_Overlap, Data_Overlap_ExtraCluster){
  to.plot.OrgAD <- cbind(Data_Overlap$dataSet, cluster = Data_Overlap$cluster)
  g <- ggplot(to.plot.OrgAD, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("Artificial Dataset") + theme_classic() + scale_color_hue("Clusters")
  
  #K-means
  print("Overlapping artificial data set")
  BIC.KMN.Overlap <- BIC.KMN.Func(Data_Overlap$dataSet, full.Overlap.3.15K.Euc$KM, 3, 15)
  to.plot.BIC.Ovr <- cbind(Data_Overlap$dataSet, cluster = full.Overlap.3.15K.Euc$KM[[BIC.KMN.Overlap]]$cluster)
  g1 <- ggplot(to.plot.BIC.Ovr, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("Kmeans with the best BIC") + theme_classic() + scale_color_hue("Clusters")
  
  print("Overlapping artificial data set with outliner")
  BIC.KMN.Overlap.Outlier <- BIC.KMN.Func(Data_Overlap_ExtraCluster$dataSet, full.Overlap.3.15K.Euc$KM, 3, 15)
  
  # EM for MOG
  print("Overlapping artificial data set")
  BIC.EM.Overlap <- BIC.EM.func(Data_Overlap$dataSet, 3, 15, "Euclidean")
  to.plot.BIC.EM.Ovr <- cbind(Data_Overlap$dataSet, cluster = BIC.EM.Overlap$classification)
  g2 <- ggplot(to.plot.BIC.EM.Ovr, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("EM for MOG with the best BIC") + theme_classic() + scale_color_hue("Clusters")
  
  print("Overlapping artificial data set with outliner")
  BIC.EM.Overlap.Outlier <- BIC.EM.func(Data_Overlap_ExtraCluster$dataSet, 3, 15, "Euclidean")
  to.plot.BIC.EM.Outlier <- cbind(Data_Overlap_ExtraCluster$dataSet, cluster = BIC.EM.Overlap.Outlier$classification)
  g3 <- ggplot(to.plot.BIC.EM.Outlier, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + 
    ggtitle("EM for MOG with the best BIC") + theme_classic() + scale_color_hue("Clusters")
  
  print("K-means and EM for MOG’s best results with BIC on Overlapping Artificial Data set")
  multiplot(g,g1,g2, cols = 3)
  
  cmc.Overlap.AD <- mixmodCombi(Data_Overlap$dataSet, nbCluster = c(3:15))
  comb.cluster <- cmc.Overlap.AD@mixmodOutput@bestResult@partition
  count(comb.cluster)
  to.plot.CMC.Ovr <- cbind(Data_Overlap$dataSet, cluster = comb.cluster)
  g4 <- ggplot(to.plot.CMC.Ovr, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + ggtitle("Combining Mixture Component Clustering") + theme_classic()
  
  print("BIC vs CMC on Overlapping artificial data set")
  
  pairwise.distance <- proxy::dist(Data_Overlap$dataSet, method = "Euclidean")
  compare <- cbind(
    rbind(BIC.EM.Overlap$loglik, cluster.stats(pairwise.distance,BIC.EM.Overlap$classification)$avg.silwidth),
    rbind(cmc.Overlap.AD@mixmodOutput@bestResult@likelihood, cluster.stats(pairwise.distance,cmc.Overlap.AD@mixmodOutput@bestResult@partition)$avg.silwidth))
  colnames(compare) <- c("BIC", "CMC")
  rownames(compare) <- c("Log likelihood", "Silhouette Width")
  
  print(compare)
  multiplot(g, g2, g4, cols = 3)
  
  cmc.Overlap.Extra.AD <- mixmodCombi(Data_Overlap_ExtraCluster$dataSet, nbCluster = c(3:15))
  cmc.Overlap.Extra.AD@mixmodOutput@bestResult@likelihood
  comb.cluster <- cmc.Overlap.Extra.AD@mixmodOutput@bestResult@partition
  count(comb.cluster)
  to.plot.CMC.Ovr.Outlier <- cbind(Data_Overlap_ExtraCluster$dataSet, cluster = comb.cluster)
  g5 <- ggplot(to.plot.CMC.Ovr.Outlier, aes(V1, V2, color = factor(cluster))) + geom_point(size = 3) + ggtitle("Combining Mixture Component Clustering") + theme_classic()
  
  print("BIC vs CMC on Overlapping artificial data set with Outliers")
  pairwise.distance <- proxy::dist(Data_Overlap_ExtraCluster$dataSet, method = "Euclidean")
  compare <- cbind(
    rbind(BIC.EM.Overlap.Outlier$loglik, cluster.stats(pairwise.distance,BIC.EM.Overlap.Outlier$classification)$avg.silwidth),
    rbind(cmc.Overlap.Extra.AD@mixmodOutput@bestResult@likelihood, cluster.stats(pairwise.distance,cmc.Overlap.Extra.AD@mixmodOutput@bestResult@partition)$avg.silwidth))
  colnames(compare) <- c("BIC", "CMC")
  rownames(compare) <- c("Log likelihood", "Silhouette Width")
  
  multiplot(g, g3, g5, cols = 3)
  
}

BIC.CMC.Real <- function(Data_Real, kmin, kmax){
  cmc.Real <- mixmodCombi(Data_Real, nbCluster = c(kmin:kmax))
  cmc.Real@mixmodOutput@bestResult@likelihood
  comb.cluster <- cmc.Real@mixmodOutput@bestResult@partition
  print("CMC returns optimal number of clusters on Real data set: ")
  print(count(comb.cluster))
}

## Slow execution
## K-means, Hierarchical, EM for MOG for K = 3...15
if (!exists(c("full.Overlap.3.15K.Euc", "full.Overlap.Outlier.3.15K.Euc"))){
  print("Find K-means with K = 3, .., 15 on Overlapping Artificial Data set")
  full.Overlap.3.15K.Euc <- collect.func(Data_Overlap$dataSet, 3,15,"Euclidean")
  print("Find K-means with K = 3, .., 15 on Overlapping Artificial Data set with Outliers")
  full.Overlap.Outlier.3.15K.Euc <- collect.func(Data_Overlap_ExtraCluster$dataSet, 3,15,"Euclidean")
}

## Experiment on Artificial data set
## Elbow method
print("Elbow methods")
elbow.func(full.Overlap.3.15K.Euc$KM,full.Overlap.3.15K.Euc$HC,full.Overlap.3.15K.Euc$EM, 3, 15)
## Gap statistic -- Really slow, due to the bootstrapping step -- Consider berfore running
print("Gap statistic methods")
Gap.Stat.Artificial(Data_Overlap, Data_WellSeparated)
## BIC & CMC
print("BIC and CMC methods")
BIC.CMC.Artificial(Data_Overlap, Data_Overlap_ExtraCluster)

yn <- quote(var1 <- readline('It will take approximate 2 hours to execute on the Real data set \n Do you want to continue (Y/N)? '))
eval(yn)
if(var1 == "Y"){
  full.Real.3.15k.Pear <- collect.func(Data_Real, 3, 15, "correlation")
  ## Experiment on Real data set
  ## Elbow method
  print("Elbow methods")
  elbow.func(full.Real.3.15k.Pear$KM, full.Real.3.15k.Pear$HC, full.Real.3.15k.Pear$EM, 3, 15)
  ## Gap statistic -- Really slow, it will take more than 2 hours, due to the bootstrapping step -- Consider berfore running
  print("Gap statistic methods")
  Gap.Stat.Real(Data_Real)
  ## BIC & CMC
  print("BIC and CMC methods")
  BIC.KMN.Func(Data_Real, full.Real.3.15k.Pear$KM, 3, 15)
  BIC.EM.func(Data_Real, 3, 15, "correlation")
  BIC.CMC.Real(Data_Real, 3, 15)
}



