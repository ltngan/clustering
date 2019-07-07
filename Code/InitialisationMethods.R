############################ Initialisation methods #################################

### Calling all needed functions

source("CommonFunc.R")
LoadLibraries()
############ Initialisation ######################################
if (!exists(c("Data_Overlap","Data_WellSeparated", "Data_Real"))){
  source("LoadData.R")
}

### Initialisation with K-means
initialisation.Kmeans <- function (data, k, dist_method) {
  ini.SSE <- data.frame(random = 0, forgy = 0, kmeanPP = 0)
  for (i in 1:20){
    # Random
    temp <- K_means(data, k, "random", dist_method)$WSS
    # Forgy
    temp <- c(temp, K_means(data, k, "Forgy", dist_method)$WSS)
    # K-means++
    temp <- c(temp, K_means(data, k, "k-means++", dist_method)$WSS)
    ini.SSE[i, ] <- temp
  }
  ini.SSE <- melt(ini.SSE)
  return(ini.SSE)
}  

initialisation.KM.Experiment <- function(Data_Overlap, Data_WellSeparated, Data_Real){
  ini.SSE.DT.Overlap <- initialisation.Kmeans(Data_Overlap$dataSet, 5, "Euclidean")  
  avgSSE <- ini.SSE.DT.Overlap %>% group_by(variable) %>% summarise(mean(value))
  print("Average SSE for Overlapping artificial data")
  print(avgSSE)
  g1 <- ggplot(ini.SSE.DT.Overlap, aes(x = rep((1:20),3), y = value, group = variable, fill = variable)) + 
    geom_bar(stat="identity", position="dodge", width=.5) + ggtitle("Ovelapping Artificial Dataset") + xlab("Trial time") + ylab("Sum of Square Error")
  
  ini.SSE.DT.WellSeparated <- initialisation.Kmeans(Data_Overlap_ExtraCluster$dataSet, 5, "Euclidean") 
  avgSSE <- ini.SSE.DT.WellSeparated %>% group_by(variable) %>% summarise(mean(value))
  print("Average SSE for Well-separated artificial data")
  print(avgSSE)
  g2 <- ggplot(ini.SSE.DT.WellSeparated, aes(x = rep((1:20),3), y = value, group = variable, fill = variable)) + 
    geom_bar(stat="identity", position="dodge", width=.5) + ggtitle("Well-Separated Artificial Dataset") + xlab("Trial time") + ylab("Sum of Square Error")
 
  yn <- quote(var1 <- readline('It will take more than 10 mininutes to execute on the Real data set \n Do you want to continue (Y/N)? '))
  eval(yn)
  if(var1 == "Y"){
    ini.SSE.DT.Real <- initialisation.Kmeans(Data_Real, 15, "Euclidean")
    avgSSE <- ini.SSE.DT.Real %>% group_by(variable) %>% summarise(mean(value))
    print("Average SSE for Real data")
    print(avgSSE)
    g3 <- ggplot(ini.SSE.DT.Real, aes(x = rep((1:20),3), y = value, group = variable, fill = variable)) + 
      geom_bar(stat="identity", position="dodge", width=.5) + ggtitle("Real Dataset") + xlab("Trial time") + ylab("Sum of Square Error")
    
    multiplot(g1, g2, g3, cols = 1)
  } else {
    multiplot(g1, g2, cols = 1)
  }
  
}


### Initialisation with EM for MOG
initialisation.EM.func <- function(data, k, dist_method) {
  ini.LLH <- data.frame(random = 0, kmean = 0)
  for (i in 1:20) {
    temp <- EM_main(data, k, "random", dist_method)$logLikelihood
    temp <- c(temp,EM_main(data, k, "k-means", dist_method)$logLikelihood)
    ini.LLH[i, ] <- temp
  }
  ini.LLH <- melt(ini.LLH)
  return(ini.LLH)
}

initialisation.EM.func.Real <- function(data, k, dist_method) {
  ini.LLH <- data.frame(kmean = 0)
  for (i in 1:20) {
    temp <- EM_main(data, k, "k-means", dist_method)$logLikelihood
    ini.LLH[i, ] <- temp
  }
  ini.LLH <- melt(ini.LLH)
  return(ini.LLH)
}
initialisation.EM.Experiment <- function(Data_Overlap, Data_WellSeparated, Data_Real){
  print("Starting EM for MOG experiment on initialisation parameters")
  ini.LLH.DT.Overlap <- initialisation.EM.func(Data_Overlap$dataSet, 5, "Euclidean")  
  avgLLH <- ini.LLH.DT.Overlap %>% group_by(variable) %>% summarise(mean(value))
  print("Average Log Likelihood for Overlapping artificial data")
  print(avgLLH)
  g1 <- ggplot(ini.LLH.DT.Overlap, aes(x = rep((1:20),2), y = value, group = variable, fill = variable)) + 
    geom_bar(stat="identity", position="dodge", width=.5) + ggtitle("Ovelapping Artificial Dataset") + xlab("Trial time") + ylab("Log Likelihood")
  
  ini.LLH.DT.WellSeparated <- initialisation.EM.func(Data_Overlap_ExtraCluster$dataSet, 5, "Euclidean") 
  avgLLH <- ini.LLH.DT.WellSeparated %>% group_by(variable) %>% summarise(mean(value))
  print("Average Log Likelihood for Well-separated artificial data")
  print(avgLLH)
  g2 <- ggplot(ini.LLH.DT.WellSeparated, aes(x = rep((1:20),2), y = value, group = variable, fill = variable)) + 
    geom_bar(stat="identity", position="dodge", width=.5) +
    ggtitle("Well-Separated Artificial Dataset") + xlab("Trial time") + ylab("Log Likelihood")
  
  yn <- quote(var1 <- readline('It will take more than 10 mininutes to execute on the Real data set \n Do you want to continue (Y/N)? '))
  eval(yn)
  if(var1 == "Y"){
    ini.LLH.DT.Real <- initialisation.EM.func.Real(Data_Real, 15, "Euclidean")
    avgLLH <- ini.LLH.DT.Real %>% group_by(variable) %>% summarise(mean(value))
    print("Average Log Likelihood for real data")
    print(avgLLH)
    g3 <- ggplot(ini.LLH.DT.Real, aes(x = c(1:20), y = value)) + 
      geom_bar(stat="identity", position="dodge", width=.5, fill = "#00CCCC") + ggtitle("Real Dataset") + xlab("Trial time") + ylab("Log Likelihood")
    
    multiplot(g1, g2, g3, cols = 1)
  } else {
    multiplot(g1, g2, cols = 1)
  }
}


### Initialisation method experiment 
## With K-means
initialisation.KM.Experiment(Data_Overlap, Data_WellSeparated, Data_Real)

