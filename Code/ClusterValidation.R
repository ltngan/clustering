############################ Cluster Validation #################################

### Calling all needed functions

### Calling all needed functions
source("CommonFunc.R")
LoadLibraries()
############ Initialisation ######################################
if (!exists(c("Data_Overlap","Data_WellSeparated", "Data_Real"))){
  source("LoadData.R")
}
if(!exists(c("full.Overlap.3.15K.Euc"))){
  print("Find K-means, Hierarchicial and EM for MOG with K = 3, .., 15 on Overlapping Artificial Data set")
  full.Overlap.3.15K.Euc <- collect.func(Data_Overlap$dataSet, 3,15,"Euclidean")
  
}
############################# Evaluation cluster ############################################

eval.AD <- function(Data_Overlap, full.Overlap.3.15K.Euc){
  ## Artificial Dataset k = 5, distance = Euclidean
  AD.pairwise_data <- proxy::dist(Data_Overlap$dataSet, method = "Euclidean")
  ## Best BIC for EM for MOG
  BIC.EM.func <- function(data, kmin, kmax, dist_method){
    BIC_EM <- Mclust(data, G = c(kmin:kmax), initialization = hclust(proxy::dist(data, method = dist_method, convert_similarities = FALSE)))
    return(BIC_EM)
  }
  BIC.EM.Overlap <- BIC.EM.func(Data_Overlap$dataSet, 3, 15, "Euclidean")
  # Internal Criterion & External Criterion
  AD.KMN.best <- cluster.stats(AD.pairwise_data, full.Overlap.3.15K.Euc$KM[[5]]$cluster, Data_Overlap$cluster)
  AD.EM.best <- cluster.stats(AD.pairwise_data,BIC.EM.Overlap$classification, Data_Overlap$cluster)
  AD.HC.best <- cluster.stats(AD.pairwise_data, full.Overlap.3.15K.Euc$HC$cluster[3,], Data_Overlap$cluster)
  
  summary.func <- function(data){
    return(rbind(SSE = data$within.cluster.ss, 
                 Silhouette =  data$avg.silwidth, 
                 PearsonGamma = data$pearsongamma, 
                 Entropy = data$entropy, 
                 RandIndex = data$corrected.rand))
  }
  print("K-means with K = 5 on Overlapping Artificial data set")
  print(summary.func(AD.KMN.best))
  
  print("Hierarchical with K = 5 on Overlapping Artificial data set")
  print(summary.func(AD.HC.best))
  
  print("EM for MOG with K = 5 on Overlapping Artificial data set")
  print(summary.func(AD.EM.best))
}

eval.Real <- function(Data_Real, full.Real.3.15k.Pear){
  ## Real Dataset k = 10, distance = "correlation"
  Real.pairwise.data <- proxy::dist(Data_Real, method = "correlation")
  ## Best BIC for EM for MOG
  BIC.EM.func <- function(data, kmin, kmax, dist_method){
    BIC_EM <- Mclust(data, G = c(kmin:kmax), initialization = hclust(proxy::dist(data, method = dist_method, convert_similarities = FALSE)))
    return(BIC_EM)
  }
  BIC.EC.Real <- BIC.EM.func(Data_Real, 3, 15, "Euclidean")
  Real.KMN.best <- cluster.stats(Real.pairwise.data, full.Real.3.15k.Pear$KM[[8]]$cluster)
  Real.EM.best <- cluster.stats(Real.pairwise.data,BIC.EC.Real$classification)
  Real.HC.best <- cluster.stats(Real.pairwise.data, full.Real.3.15k.Pear$HC$cluster[8,])
  print("K-means with K = 5 on Overlapping Artificial data set")
  print(summary.func(Real.KMN.best))
  
  print("Hierarchical with K = 5 on Overlapping Artificial data set")
  print(summary.func(Real.HC.best))
  
  print("EM for MOG with K = 5 on Overlapping Artificial data set")
  print(summary.func(Real.EM.best))
}

### Cluster validation
print("Cluster validation on Overlapping artificial data set")
eval.AD(Data_Overlap, full.Overlap.3.15K.Euc)

yn <- quote(var1 <- readline('It will take approximate 2 hours to execute on the Real data set \n Do you want to continue (Y/N)? '))
eval(yn)
if(var1 == "Y"){
  if (!exists("full.Real.3.15k.Pear")){
    print("Find K-means, Hierarchicial and EM for MOG with K = 3, .., 15 on real data set")
    full.Real.3.15k.Pear <- collect.func(Data_Real, 3, 15, "correlation")
  }
  print("Cluster validation on real data set")
  eval.Real(Data_Real, full.Real.3.15k.Pear)
}

