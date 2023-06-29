#Purpose: data generation for two-level models with count outcome and continuous x and W.
#####################Data simulation######################
start <- Sys.time()
#setwd("G:/My Drive/Classes/Simulation/Mini project/Models/Simulation data")

library(MASS)
library(lme4)
#install.packages("writexl")

set.seed(77841)
reps<-1
#Manipulated variables
Cluster_list <- c(30,50,100)
Gamma11_list <- c(0.2, 0.45)
T00_list <- c(0.004, 0.0015, 0.001)

#Cluster_list <- c(30)
#Gamma11_list <- c(0.09)
#T00_list <- c(0.004)
#Fixed variables
Gamma00_list <- c(2) 
Gamma01_list <- c(0.5)
Gamma10_list <- c(0.5)
ClusterSize_list <- c(20)


Conditions <-
  length(Cluster_list) * length(Gamma11_list)* length(T00_list) * length(Gamma00_list)* length(Gamma01_list)*length(Gamma10_list)*length(ClusterSize_list)

design <- expand.grid("Cluster" = Cluster_list, "Gamma11" = Gamma11_list, "T00" = T00_list, "Gamma00" = Gamma00_list, "Gamma10" = Gamma10_list, "ClusterSize" = ClusterSize_list)

write.csv(design,"design.csv", row.names=F)

Data<- list()

design <- numeric()
for (Cluster in Cluster_list) {
  for (Gamma11 in Gamma11_list){
    for (T00 in T00_list) {
      for (Gamma00 in Gamma00_list) {
        for (Gamma01 in Gamma01_list){
          for(Gamma10 in Gamma10_list){
            for (ClusterSize in ClusterSize_list){
              design <- rbind(design,matrix(data=rep(c(Cluster, Gamma11, T00, Gamma00, Gamma01,Gamma10, ClusterSize),reps),nrow = reps,byrow = T))
            }
          }
        }
      }
    }
  }
}


colnames(design) <- c("Cluster","Gamma11", "T00", "Gamma00", "Gamma01","Gamma10","Cluster Size") 


Data_generation <- function (Cluster, Gamma11, T00, Gamma00, Gamma01, Gamma10, ClusterSize) #function of all design factors
  {   w <- rnorm(n=Cluster, 0, 1)
        x <- rnorm(n=Cluster*ClusterSize, 0, 1)
        u0 <- rnorm(n=Cluster,0, T00) #The variance of T00 from the list of T00
         lv1_2 <- matrix(NA, nrow=Cluster*ClusterSize, ncol=6)
      
      case <- 0
      for (j in 1:Cluster) {
        for (i in 1:ClusterSize) {
          case <- case +1
          lamda <- exp(Gamma00+Gamma01*w+ Gamma10*x + Gamma11*w*x + u0)
          Y <- rpois(1, lambda = lamda[case])
          lv1_2[case,] <- c(j, w[j], x[case], u0[j],lamda[case], Y)}
  }
  data <- as.data.frame(lv1_2)
  names(data)[1] <- "Cluster_ID"
  names(data)[2] <- "w"
  names(data)[3] <- "x"
  names(data)[4] <- "u0"
  names(data)[5] <- "lamda"
  names(data)[6] <- "Y"
  
  
  return(data)
}


for (k in 1:(Conditions*reps)) {
  
  Cluster<- design[k,1]
  Gamma11 <- design[k,2]
  T00 <- design[k,3]
  Gamma00 <- design[k,4]
  Gamma01 <- design[k,5]
  Gamma10 <- design[k,6]
  ClusterSize <- design[k,7]
  
  Data[[k]] <- Data_generation (Cluster, Gamma11, T00, Gamma00, Gamma01, Gamma10, ClusterSize)
  
  cat("Just completed", k, "of", Conditions*reps, "Replications", "(","remaining", Conditions*reps-k,")", "\n")
}

save(Data,file=paste0("Data mini project.Rdata"))


end <- Sys.time()
end-start
#####################End of Data simulation######################



#####################Performance measures########################
start1 <- Sys.time()
setwd("dirname(rstudioapi::getActiveDocumentContext()$path)")
library(nlme)
library(Matrix)
library(lme4)
install.packages("doParallel")
library(doParallel)
#install.packages("Rtools")

registerDoParallel(cores=8) 
load(file=paste0("Data mini project.Rdata"))

#container for simulation
Sim <- foreach(k=1:(reps*Conditions),.combine = "rbind",.packages = c("lme4"), .inorder=T) %dopar% {
  
  Cluster<- design[k,1]
  Gamma11 <- design[k,2]
  T00 <- design[k,3]
  Gamma00 <- design[k,4]
  Gamma01 <- design[k,5]
  Gamma10 <- design[k,6]
  ClusterSize <- design[k,7]
  
  
  data_all<-Data[[k]]
 #data_all<-paste0("Data mini project.Rdata")
  raw <- rep(NA,33) #placeholder of 33

  #Q2: check model align with design  
  error1 <- try (model <- glmer(Y ~ 1+ x + w + w*x + (1 | Cluster_ID),   #check model specification
                 data=data_all, 
                 family=poisson))
  
  #model <- glmer(Y ~ 1+ x + w + w*x + (1 | Cluster_ID), 
                                #data=data_all, 
                                #family=poisson)
summary(model)
    
  if(class(error1) == "try-error") {
    raw <- c(rep(NA, 33))}
  else {
    raw <- c(AIC(model),BIC(model),
             summary(model)$coefficients[,1], #ES: intercept, x, w, x*w
             summary(model)$coefficients[,2], #standard errors
             summary(model)$coefficients[,3], #z scores
             summary(model)$coefficients[,4]) # p value 
  }
  return(c(design[k,],raw))
}
  
#To produce the result summary for the model fitting
AIC <- as.vector((by(Sim[,8],ceiling(1:(Conditions*reps)/reps), function(x) mean(x,na.rm=T))))
BIC <- as.vector((by(Sim[,9],ceiling(1:(Conditions*reps)/reps), function(x) mean(x,na.rm=T))))

bias.ES <- as.vector(by((Sim[,10]-Sim[,4]), ceiling(1:(Conditions*reps)/reps),function(x) mean(x,na.rm=T))) #b00 bias
MSE.ES <-  as.vector(by((Sim[,10]-Sim[,4])^2, ceiling(1:(Conditions*reps)/reps),function(x) mean(x,na.rm=T))) #MSE of b00
cp.ES <- as.vector(by((Sim[,10]-1.96*Sim[,13])<Sim[,4] & (Sim[,10]+1.96*Sim[,13])>Sim[,4], ceiling(1:(Conditions*reps)/reps),function(x) mean(x,na.rm=T))) #coverage
power.ES <- as.vector(by(Sim[,13]<0.05, ceiling(1:(Conditions*reps)/reps),function(x) mean(x,na.rm=T)))



## Convergence ##
convergence.rate <- 1-as.vector(by(Sim[,18], ceiling(1:(Conditions*reps)/reps),function(x) sum(is.na(x))))/reps




##save results ##
Results<-cbind(unique(design),
               AIC,BIC,bias.ES,MSE.ES,cp.ES,power,convergence.rate)



colnames(Sim) <- c("Cluster","Gamma11","T00","Gamma00",
                   " Gamma01","Gamma10","ClusterSize",
                   
                   "AIC","BIC")



write.csv(Sim, file="Results_raw.csv")
write.csv(Results, file="Results.csv")
save.image(file="Results.RData")
end <- Sys.time()
end - start 








