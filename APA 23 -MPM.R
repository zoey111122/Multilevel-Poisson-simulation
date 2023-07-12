#Purpose: data generation for two-level models with count outcome and continuous x and W.
#####################Step 1:Data simulation######################
library(MASS)
library(lme4)
install.packages("writexl")

set.seed(789654)
reps <- 1000
# Manipulated variables
Cluster_list <- c(5,10,15,20)
Gamma11_list <- c(0.2,0.45)
T00_list <- c(0.004, 0.0015, 0.001)
# Fixed variables
Gamma00_list <- c(2) 
Gamma01_list <- c(0.5)
Gamma10_list <- c(0.5)
ClusterSize_list <- c(5,10,20)

Conditions <- expand.grid("Cluster" = Cluster_list, "Gamma11" = Gamma11_list, "T00" = T00_list, "Gamma00" = Gamma00_list, "Gamma10" = Gamma10_list, "ClusterSize" = ClusterSize_list)
Conditions$number <- 1:nrow(Conditions)

Data <- list()
design <- numeric()

for (Cluster in Cluster_list) {
  for (Gamma11 in Gamma11_list){
    for (T00 in T00_list) {
      for (Gamma00 in Gamma00_list) {
        for (Gamma01 in Gamma01_list){
          for(Gamma10 in Gamma10_list){
            for (ClusterSize in ClusterSize_list){
              design <- rbind(design, matrix(data = rep(c(Cluster, Gamma11, T00, Gamma00, Gamma01, Gamma10, ClusterSize), reps), nrow = reps, byrow = T))
            }
          }
        }
      }
    }
  }
}

colnames(design) <- c("Cluster","Gamma11", "T00", "Gamma00", "Gamma01", "Gamma10", "Cluster Size")

Data_generation <- function(Cluster, Gamma11, T00, Gamma00, Gamma01, Gamma10, ClusterSize) {
  w <- rnorm(n = Cluster, 0, 1)
  x <- rnorm(n = Cluster * ClusterSize, 0, 1)
  u0 <- rnorm(n = Cluster, 0, T00)
  lv1_2 <- matrix(NA, nrow = Cluster * ClusterSize, ncol = 6)  # Update the number of columns
  case <- 0
  for (j in 1:Cluster) {
    for (i in 1:ClusterSize) {
      case <- case + 1
      lamda <- exp(Gamma00 + Gamma01 * w[j] + Gamma10 * x[case] + Gamma11 * w[j] * x[case] + u0[j])
      if (length(lamda) > 1) {
        lamda <- mean(lamda)
      }
      Y <- rpois(1, lambda = lamda)
      lv1_2[case, ] <- c(j, w[j], x[case], u0[j], lamda, Y)  # Remove Y_ceiling from the matrix
    }
  }
  data_1 <- as.data.frame(lv1_2)
  names(data_1) <- c("Cluster_ID", "w", "x", "u0", "lamda", "Y")  # Update column names
  return(data_1)
}

# Create a list of lists for the data
Data <- vector("list", nrow(Conditions))
names(Data) <- paste0("Condition_", 1:nrow(Conditions))

for (i in 1:nrow(Conditions)) {
  Data[[i]] <- list()  # Create an empty list for each condition
  
  for (k in 1:reps) {  
    Cluster <- design[(i-1)*reps + k,1]
    Gamma11 <- design[(i-1)*reps + k,2]
    T00 <- design[(i-1)*reps + k,3]
    Gamma00 <- design[(i-1)*reps + k,4]
    Gamma01 <- design[(i-1)*reps + k,5]
    Gamma10 <- design[(i-1)*reps + k,6]
    ClusterSize <- design[(i-1)*reps + k,7]
    
    Data[[i]][[k]] <- Data_generation(Cluster, Gamma11, T00, Gamma00, Gamma01, Gamma10, ClusterSize)
    
    cat("Just completed", k, "of", reps, "Replications for Condition", i, "(", "remaining", reps - k, ")\n")
  }
}

#####################Step 1:End of Data simulation######################



#####################Step 2:Performance measures########################
#install.packages("rstudioapi")
library(rstudioapi)
#setwd("dirname(rstudioapi::getActiveDocumentContext()$path)")
library(nlme)
library(Matrix)
library(lme4)
#install.packages("doParallel")
library(doParallel)
#install.packages("Rtools")
# load(file=paste0("Condition 1.Rdata"))
#registerDoParallel(cores=8) 

library(lme4)

# Initialize MLM list of arrays
MLM <- vector("list", length = nrow(Conditions))
for (i in 1:nrow(Conditions)) {
  MLM[[i]] <- array(dim=c(reps,10))
}

non_converge <- 0

# Initialize result list
for (i in 1:nrow(Conditions)) {
  MLM[[i]] <- array(dim = c(reps, 10))
  for (J in 1:reps) {
    error1 <- try (
      model <- glmer(Y ~ 1 + x + w + w * x + (1 | Cluster_ID),
                     data = Data[[i]][[J]],
                     family = poisson),
      silent = TRUE
    )
    if (class(error1) == "try-error") {
      non_converge <- non_converge + 1
    } else {
      summary_model <- summary(model)
      # To produce the result summary for the model fitting
      par.est.int <- summary_model$coefficients[1, 1]   # ES estimate for the intercept
      par.est.x <- summary_model$coefficients[2, 1]  # ES estimate for the x
      par.est.w <- summary_model$coefficients[3, 1]  # ES estimate for the w
      par.est.xw <- summary_model$coefficients[4, 1]  # ES estimate for the x*w
      
      p.value.int <- summary_model$coefficients[1, 4]    # p-value for the intercept
      p.value.x <- summary_model$coefficients[2, 4]   # p-value for the x
      p.value.w <- summary_model$coefficients[3, 4]   # p-value for the w
      p.value.xw <- summary_model$coefficients[4, 4]   # p-value for the x*w
      
      MLM[[i]][J, 1] <- par.est.int       # ES estimate for the intercept
      MLM[[i]][J, 2] <- par.est.x      # ES estimate for the x
      MLM[[i]][J, 3] <- par.est.w     # ES estimate for the w
      MLM[[i]][J, 4] <- par.est.xw     # ES estimate for the xw interaction
      MLM[[i]][J, 5] <- p.value.int       # p-value for the intercept
      MLM[[i]][J, 6] <- p.value.x      # p-value for the x
      MLM[[i]][J, 7] <- p.value.w     # p-value for the w
      MLM[[i]][J, 8] <- p.value.xw     # p-value for the xw interaction
    }
    cat("Just completed", J, "of", reps, "simulation for Condition", i, "(", "remaining", reps - J, ")", "\n")
  }
}

#############Step 2:End of performance analysis########################




#############Step 3: Results summary####################################
#save results
MLM.ES.int <- vector("numeric", length = nrow(Conditions))
MLM.ES.x <- vector("numeric", length = nrow(Conditions))
MLM.ES.w <- vector("numeric", length = nrow(Conditions))
MLM.ES.xw <- vector("numeric", length = nrow(Conditions))
power.MLM.x <- vector("numeric", length = nrow(Conditions))
power.MLM.w <- vector("numeric", length = nrow(Conditions))
power.MLM.xw <- vector("numeric", length = nrow(Conditions))
converge.rate <- vector("numeric", length = nrow(Conditions))

for(i in 1:nrow(Conditions)){
  MLM.ES.int[i] <- mean(MLM[[i]][,1], na.rm = T)
  MLM.ES.x[i] <- mean(MLM[[i]][,2], na.rm = T)
  MLM.ES.w[i] <- mean(MLM[[i]][,3], na.rm = T)
  MLM.ES.xw[i] <- mean(MLM[[i]][,4], na.rm = T)
  power.MLM.x[i] <- mean(MLM[[i]][,6] < 0.05, na.rm = T)
  power.MLM.w[i] <- mean(MLM[[i]][,7] < 0.05, na.rm = T)
  power.MLM.xw[i] <- mean(MLM[[i]][,8] < 0.05, na.rm = T)
  converge.rate[i] <- (reps - non_converge)/reps
}


results_3 <- t(rbind(Conditions$Cluster, Conditions$Gamma11, Conditions$T00,Conditions$Gamma00, Conditions$Gamma10, Conditions$ClusterSize,
                     Conditions$number, MLM.ES.int, MLM.ES.x, MLM.ES.w, MLM.ES.xw, power.MLM.x, power.MLM.w, power.MLM.xw,converge.rate))
results_3 <- as.data.frame(results_3)

write.csv(results_3, file = "Results_samll_sample.csv",row.names=F)




