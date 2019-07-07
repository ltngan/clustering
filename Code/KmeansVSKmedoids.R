############################ Initialisation methods #################################

### Calling all needed functions

source("CommonFunc.R")
LoadLibraries()
############ Initialisation ######################################
if (!exists(c("Data_Overlap","Data_WellSeparated", "Data_Real"))){
  source("LoadData.R")
}

print("The next step will take more than 3 hours to finish")
yn <- quote(var1 <- readline('It will take more than 3 hours to run K-medoids on the Real data set \n Do you want to continue (Y/N)? '))
eval(yn)
if(var1 == "Y"){
############################# K-means vs K-medoids #########################
  kmnVSkmd <- function(Data_Real, full.Real.3.15k.Pear){
    Real.pairwise.data <- proxy::dist(Data_Real, method = "correlation")
    print("Executing K-means with K = 8")
    Real.KMN <- K_means(Data_Real, 8, "k-means++", "correlation")
    Real.KMN.best <- cluster.stats(Real.pairwise.data,Real.KMN$cluster)
    
    print("Executing K-medoids with K = 8")
    Real.KMD <- K_medoids(Data_Real, 8, "k-means++", "correlation")
    Real.KMD.best <- cluster.stats(Real.pairwise.data,Real.KMD$cluster)
    
    data.to.plot <- cbind(id = c(1:nrow(Data_Real)), Data_Real, kmedoids = Real.KMD$cluster, kmeans = Real.KMN$cluster)
    data.to.plot <- melt(data.to.plot, id.vars = c("id", "kmedoids", "kmeans"))
    
    plotKmn <- list()
    for (i in 1:8){
      temp.to.plot <- data.to.plot[kmeans == i,]
      plotKmn[[i]] <- ggplot(temp.to.plot, aes(x=variable, y = value, group = id)) + 
        geom_line(colour = "royalblue") + xlab("time") + ylab("Gene expression") +
        stat_summary(fun.y = mean, aes(group = 1), geom = "line", size = 3, colour = "red")
    }
    print("Best K-means cluster")
    multiplot(plotlist = plotKmn, cols = 8)
    
    plotKmd <- list()
    for (i in 1:8){
      temp.to.plot <- data.to.plot[kmedoids == i,]
      plotKmd[[i]] <- ggplot(temp.to.plot, aes(x=variable, y = value, group = id)) + 
        geom_line(colour = "royalblue") + xlab("time") + ylab("Gene expression") +
        stat_summary(fun.y = mean, aes(group = 1), geom = "line", size = 3, colour = "red")
    }
    print("Best K-medoids cluster")
    multiplot(plotlist = plotKmd, cols = 8)
  }
  
  ### K-means vs K-medoids
  kmnVSkmd(Data_Real)
}
