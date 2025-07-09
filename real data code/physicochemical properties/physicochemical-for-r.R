
rm(list=ls())
library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)

cores <- 48#                               
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("main functions.R")
})
tau=c(1/6, 2/6, 3/6, 4/6, 5/6)
data <- read.csv("CASP.csv",head=TRUE)
x <- as.matrix(data[1:45730,2:10])
y <- as.matrix(data[1:45730,1])
x <- apply(x, 2, as.numeric)
y <- apply(y, 2, as.numeric)
set.seed(123)
train_indices <- sample(nrow(x), 40000)
x_train <- x[train_indices, ]#dim(x_train)
y_train <- y[train_indices, ]
x_test <- x[-train_indices, ]
y_test <- y[-train_indices, ]
means <- apply(x_train, 2, mean)
sds <- apply(x_train, 2, sd)
x_train <- scale(x_train, center = means, scale = sds)
x_test <- scale(x_test, center = means, scale = sds)

n <- N <- dim(x_train)[1]
p <- dim(x)[2]
H <- 200   
r2 <- 200 
r3 <- 300 
r4 <- 400 
r5 <- 500 
library(cqrReg)
beta_full <- cqr.fit(x_train,y_train,tau,method="ip")$beta

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  M <- 10
  
  pe_full <- MSE(y_test, x_test%*%beta_full)
  
  beta_U_r2 <- optsample(y_train, x_train, r2, 'U',tau=tau)
  pe_U2_r2 <- MSE(y_test, x_test%*%beta_U_r2)
  beta_L_r2 <- optsample(y_train, x_train, r2, 'L',tau=tau)
  pe_L2_r2 <- MSE(y_test, x_test%*%beta_L_r2)
  beta_per_pois_r2 <- sumwcqrfit(N,p,y_train, x_train, M, r2,method='pois',tau=tau)
  pe_per_pois2_r2 <- MSE(y_test, x_test%*%beta_per_pois_r2)#
  beta_ls_per_pois_r2 <- sumwlsfit(N,p,y_train, x_train, M, r2)
  pe_ls_per_pois2_r2 <- MSE(y_test, x_test%*%beta_ls_per_pois_r2)
  
  beta_U_r3 <- optsample(y_train, x_train, r3, 'U',tau=tau)
  pe_U2_r3 <- MSE(y_test, x_test%*%beta_U_r3)
  beta_L_r3 <- optsample(y_train, x_train, r3, 'L',tau=tau)
  pe_L2_r3 <- MSE(y_test, x_test%*%beta_L_r3)
  beta_per_pois_r3 <- sumwcqrfit(N,p,y_train, x_train, M, r3,method='pois',tau=tau)
  pe_per_pois2_r3 <- MSE(y_test, x_test%*%beta_per_pois_r3)#
  beta_ls_per_pois_r3 <- sumwlsfit(N,p,y_train, x_train, M, r3)
  pe_ls_per_pois2_r3 <- MSE(y_test, x_test%*%beta_ls_per_pois_r3)
  
  beta_U_r4 <- optsample(y_train, x_train, r4, 'U',tau=tau)
  pe_U2_r4 <- MSE(y_test, x_test%*%beta_U_r4)
  beta_L_r4 <- optsample(y_train, x_train, r4, 'L',tau=tau)
  pe_L2_r4 <- MSE(y_test, x_test%*%beta_L_r4)
  beta_per_pois_r4 <- sumwcqrfit(N,p,y_train, x_train, M, r4,method='pois',tau=tau)
  pe_per_pois2_r4 <- MSE(y_test, x_test%*%beta_per_pois_r4)#
  beta_ls_per_pois_r4 <- sumwlsfit(N,p,y_train, x_train, M, r4)
  pe_ls_per_pois2_r4 <- MSE(y_test, x_test%*%beta_ls_per_pois_r4)
  
  beta_U_r5 <- optsample(y_train, x_train, r5, 'U',tau=tau)
  pe_U2_r5 <- MSE(y_test, x_test%*%beta_U_r5)
  beta_L_r5 <- optsample(y_train, x_train, r5, 'L',tau=tau)
  pe_L2_r5 <- MSE(y_test, x_test%*%beta_L_r5)
  beta_per_pois_r5 <- sumwcqrfit(N,p,y_train, x_train, M, r5,method='pois',tau=tau)
  pe_per_pois2_r5 <- MSE(y_test, x_test%*%beta_per_pois_r5)#
  beta_ls_per_pois_r5 <- sumwlsfit(N,p,y_train, x_train, M, r5)
  pe_ls_per_pois2_r5 <- MSE(y_test, x_test%*%beta_ls_per_pois_r5)
  
  
  c(pe_full, pe_U2_r2, pe_L2_r2, pe_per_pois2_r2, pe_ls_per_pois2_r2,
    pe_full, pe_U2_r3, pe_L2_r3, pe_per_pois2_r3, pe_ls_per_pois2_r3,
    pe_full, pe_U2_r4, pe_L2_r4, pe_per_pois2_r4, pe_ls_per_pois2_r4,
    pe_full, pe_U2_r5, pe_L2_r5, pe_per_pois2_r5, pe_ls_per_pois2_r5) 
} 

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("pe_full", "pe_U2_r2", "pe_L2_r2", "pe_per_pois2_r2", "pe_ls_per_pois2_r2",
                         "pe_full", "pe_U2_r3", "pe_L2_r3", "pe_per_pois2_r3", "pe_ls_per_pois2_r3",
                         "pe_full", "pe_U2_r4", "pe_L2_r4", "pe_per_pois2_r4", "pe_ls_per_pois2_r4",
                         "pe_full", "pe_U2_r5", "pe_L2_r5", "pe_per_pois2_r5", "pe_ls_per_pois2_r5")
result_df
#write.csv(result_df,"E:\\CQR\\real_data_analysis\\Physicochemical\\pe_for-r.csv")










rm(list=ls())
library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)

cores <- 48#                               
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("main functions.R")
})
tau5=c(1/6, 2/6, 3/6, 4/6, 5/6)
data <- read.csv("CASP.csv",head=TRUE)
x <- as.matrix(data[1:45730,2:10])
y <- as.matrix(data[1:45730,1])
x <- apply(x, 2, as.numeric)
y <- apply(y, 2, as.numeric)
set.seed(123)
train_indices <- sample(nrow(x), 40000)
x_train <- x[train_indices, ]#dim(x_train)
y_train <- y[train_indices, ]
x_test <- x[-train_indices, ]
y_test <- y[-train_indices, ]
means <- apply(x_train, 2, mean)
sds <- apply(x_train, 2, sd)
x_train <- scale(x_train, center = means, scale = sds)
x_test <- scale(x_test, center = means, scale = sds)


n <- N <- dim(x_train)[1]
p <- dim(x_train)[2]
H <- 200    
r <- 300 
result <- foreach(h = 1:H, 
                  .packages = c("mvtnorm", "Rfit", "cqrReg"),  
                  .combine = "rbind") %dopar% {
                    
                    M <- 10
                    beta_est <- sumwcqrfit(N, p, y_train, x_train, M, r, method = 'pois', tau = tau5)
                    
                    if(is.null(names(beta_est))) {
                      names(beta_est) <- paste0("beta", 1:p)
                    }
                    
                    beta_est
                  }

stopCluster(cl)
result_df <- as.data.frame(result)
result_df$Iteration <- 1:H
#write.csv(result_df,"E:\\CQR\\real_data_analysis\\Physicochemical\\estimator.csv")








