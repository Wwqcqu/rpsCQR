
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
tau1=1/2
tau3=c(1/4, 2/4, 3/4)
tau5=c(1/6, 2/6, 3/6, 4/6, 5/6)
tau7=c(1/8,2/8,3/8,4/8,5/8,6/8,7/8)
tau9=c(1/10,2/10,3/10,4/10,5/10,6/10,7/10,8/10,9/10)
data <- read.csv("superconductivty_train.csv",head=TRUE)
x <- as.matrix(data[1:21263,1:81])
y <- as.matrix(data[1:21263,82])

y=log(y)
x <- apply(x, 2, as.numeric)
y <- apply(y, 2, as.numeric)

set.seed(123)
train_indices <- sample(nrow(x), 20000)
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
H <- 200    ##
r <- 300 
library(cqrReg)


result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  M <- 10
  
  beta_full_K1 <- cqr.fit(x_train,y_train,tau1,method="ip")$beta
  pe_full_K1 <- MSE(y_test, x_test%*%beta_full_K1)
  beta_U_K1 <- optsample(y_train, x_train, r, 'U',tau=tau1)
  pe_U2_K1 <- MSE(y_test, x_test%*%beta_U_K1)
  beta_L_K1 <- optsample(y_train, x_train, r, 'L',tau=tau1)
  pe_L2_K1 <- MSE(y_test, x_test%*%beta_L_K1)
  beta_per_pois_K1 <- sumwcqrfit(N,p,y_train, x_train, M, r,method='pois',tau=tau1)
  pe_per_pois2_K1 <- MSE(y_test, x_test%*%beta_per_pois_K1)#
  beta_ls_per_pois <- sumwlsfit(N,p,y_train, x_train, M, r)
  pe_ls_per_pois2 <- MSE(y_test, x_test%*%beta_ls_per_pois)
  
  beta_full_K3 <- cqr.fit(x_train,y_train,tau3,method="ip")$beta
  pe_full_K3 <- MSE(y_test, x_test%*%beta_full_K3)
  beta_U_K3 <- optsample(y_train, x_train, r, 'U',tau=tau3)
  pe_U2_K3 <- MSE(y_test, x_test%*%beta_U_K3)
  beta_L_K3 <- optsample(y_train, x_train, r, 'L',tau=tau3)
  pe_L2_K3 <- MSE(y_test, x_test%*%beta_L_K3)
  beta_per_pois_K3 <- sumwcqrfit(N,p,y_train, x_train, M, r,method='pois',tau=tau3)
  pe_per_pois2_K3 <- MSE(y_test, x_test%*%beta_per_pois_K3)#
  
  beta_full_K5 <- cqr.fit(x_train,y_train,tau5,method="ip")$beta
  pe_full_K5 <- MSE(y_test, x_test%*%beta_full_K5)
  beta_U_K5 <- optsample(y_train, x_train, r, 'U',tau=tau5)
  pe_U2_K5 <- MSE(y_test, x_test%*%beta_U_K5)
  beta_L_K5 <- optsample(y_train, x_train, r, 'L',tau=tau5)
  pe_L2_K5 <- MSE(y_test, x_test%*%beta_L_K5)
  beta_per_pois_K5 <- sumwcqrfit(N,p,y_train, x_train, M, r,method='pois',tau=tau5)
  pe_per_pois2_K5 <- MSE(y_test, x_test%*%beta_per_pois_K5)#
  
  beta_full_K7 <- cqr.fit(x_train,y_train,tau7,method="ip")$beta
  pe_full_K7 <- MSE(y_test, x_test%*%beta_full_K7)
  beta_U_K7 <- optsample(y_train, x_train, r, 'U',tau=tau7)
  pe_U2_K7 <- MSE(y_test, x_test%*%beta_U_K7)
  beta_L_K7 <- optsample(y_train, x_train, r, 'L',tau=tau7)
  pe_L2_K7 <- MSE(y_test, x_test%*%beta_L_K7)
  beta_per_pois_K7 <- sumwcqrfit(N,p,y_train, x_train, M, r,method='pois',tau=tau7)
  pe_per_pois2_K7 <- MSE(y_test, x_test%*%beta_per_pois_K7)#
  
  beta_full_K9 <- cqr.fit(x_train,y_train,tau9,method="ip")$beta
  pe_full_K9 <- MSE(y_test, x_test%*%beta_full_K9)
  beta_U_K9 <- optsample(y_train, x_train, r, 'U',tau=tau9)
  pe_U2_K9 <- MSE(y_test, x_test%*%beta_U_K9)
  beta_L_K9 <- optsample(y_train, x_train, r, 'L',tau=tau9)
  pe_L2_K9 <- MSE(y_test, x_test%*%beta_L_K9)
  beta_per_pois_K9 <- sumwcqrfit(N,p,y_train, x_train, M, r,method='pois',tau=tau9)
  pe_per_pois2_K9 <- MSE(y_test, x_test%*%beta_per_pois_K9)#
  
  
  c(pe_full_K1,pe_U2_K1, pe_L2_K1,pe_per_pois2_K1,pe_ls_per_pois2,
    pe_full_K3,pe_U2_K3, pe_L2_K3,pe_per_pois2_K3,pe_ls_per_pois2,
    pe_full_K5,pe_U2_K5, pe_L2_K5,pe_per_pois2_K5,pe_ls_per_pois2,
    pe_full_K7,pe_U2_K7, pe_L2_K7,pe_per_pois2_K7,pe_ls_per_pois2,
    pe_full_K9,pe_U2_K9, pe_L2_K9,pe_per_pois2_K9,pe_ls_per_pois2) 
} 

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("pe_full_K1","pe_U2_K1", "pe_L2_K1","pe_per_pois2_K1","pe_ls_per_pois2",
                         "pe_full_K3","pe_U2_K3", "pe_L2_K3","pe_per_pois2_K3","pe_ls_per_pois2",
                         "pe_full_K5","pe_U2_K5", "pe_L2_K5","pe_per_pois2_K5","pe_ls_per_pois2",
                         "pe_full_K7","pe_U2_K7", "pe_L2_K7","pe_per_pois2_K7","pe_ls_per_pois2",
                         "pe_full_K9","pe_U2_K9", "pe_L2_K9","pe_per_pois2_K9","pe_ls_per_pois2")
result_df
log10(result_df)
#write.csv(result_df,"E:\\CQR\\real_data_analysis\\superconductivty\\pe-r300-for-K.csv")





