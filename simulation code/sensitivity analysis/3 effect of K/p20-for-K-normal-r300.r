rm(list=ls())

library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)

cores <- 48
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("main functions.R")
})
tau1=1/2
tau5=c(1/6, 2/6, 3/6, 4/6, 5/6)
tau10=c(1/11,2/11,3/11,4/11,5/11,6/11,7/11,8/11,9/11,10/11)
beta_0 <- c(0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2)
#beta_0 <- c(0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,
            #0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2)

n <- N <- 10000
p <- length(beta_0)

H <- 200    ##
r <- 300

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "mean0-Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  
  beta_full_K1 <- cqr.fit(X=x,y,tau1,method="ip")$beta
  mse_full_K1 <- MSE(beta_full_K1, beta_0)
  M <- 10
  beta_U_K1 <- optsample(y, x, r, 'U',tau=tau1)
  mse_U2_K1 <- MSE(beta_0, beta_U_K1)
  beta_L_K1 <- optsample(y, x, r, 'L',tau=tau1)
  mse_L2_K1 <- MSE(beta_0, beta_L_K1)
  beta_per_pois_K1 <- sumwcqrfit(N,p,y, x, M, r,method='pois',tau=tau1)
  mse_per_pois2_K1 <- MSE(beta_0, beta_per_pois_K1)#20秒
  beta_U_F_K1 <- Unisample(y, x, r,tau1)
  mse_U_F2_K1 <- MSE(beta_0, beta_U_F_K1)
  beta_lev_K1 <- lev_subsample(N,p,y, x, r,tau1)###
  mse_lev2_K1 <- MSE(beta_0, beta_lev_K1)
  beta_ls_per_pois <- sumwlsfit(N,p,y, x, M, r)
  mse_ls_per_pois2 <- MSE(beta_0, beta_ls_per_pois)
  
  beta_full_K5 <- cqr.fit(X=x,y,tau5,method="ip")$beta
  mse_full_K5 <- MSE(beta_full_K5, beta_0)
  beta_U_K5 <- optsample(y, x, r, 'U',tau=tau5)
  mse_U2_K5 <- MSE(beta_0, beta_U_K5)
  beta_L_K5 <- optsample(y, x, r, 'L',tau=tau5)
  mse_L2_K5 <- MSE(beta_0, beta_L_K5)
  beta_per_pois_K5 <- sumwcqrfit(N,p,y, x, M, r,method='pois',tau=tau5)
  mse_per_pois2_K5 <- MSE(beta_0, beta_per_pois_K5)#20秒
  beta_U_F_K5 <- Unisample(y, x, r,tau5)
  mse_U_F2_K5 <- MSE(beta_0, beta_U_F_K5)
  beta_lev_K5 <- lev_subsample(N,p,y, x, r,tau5)###
  mse_lev2_K5 <- MSE(beta_0, beta_lev_K5)
  
  beta_full_K10 <- cqr.fit(X=x,y,tau10,method="ip")$beta
  mse_full_K10 <- MSE(beta_full_K10, beta_0)
  beta_U_K10 <- optsample(y, x, r, 'U',tau=tau10)
  mse_U2_K10 <- MSE(beta_0, beta_U_K10)
  beta_L_K10 <- optsample(y, x, r, 'L',tau=tau10)
  mse_L2_K10 <- MSE(beta_0, beta_L_K10)
  beta_per_pois_K10 <- sumwcqrfit(N,p,y, x, M, r,method='pois',tau=tau10)
  mse_per_pois2_K10 <- MSE(beta_0, beta_per_pois_K10)#20秒
  beta_U_F_K10 <- Unisample(y, x, r,tau10)
  mse_U_F2_K10 <- MSE(beta_0, beta_U_F_K10)
  beta_lev_K10 <- lev_subsample(N,p,y, x, r,tau10)###
  mse_lev2_K10 <- MSE(beta_0, beta_lev_K10)
  
  #beta_rr_per_pois <- sumwrfit(N,p,y, x, M, r)
  #mse_rr_per_pois2 <- MSE(beta_0, beta_rr_per_pois)
  
    c(mse_full_K1,mse_U2_K1, mse_U_F2_K1,mse_lev2_K1,mse_L2_K1,mse_per_pois2_K1,mse_ls_per_pois2,
      mse_full_K5,mse_U2_K5, mse_U_F2_K5,mse_lev2_K5,mse_L2_K5,mse_per_pois2_K5,mse_ls_per_pois2,
      mse_full_K10,mse_U2_K10, mse_U_F2_K10,mse_lev2_K10,mse_L2_K10,mse_per_pois2_K10,mse_ls_per_pois2) 
} 

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("mse_full_K1","mse_U2_K1","mse_U_F2_K1","mse_lev2_K1","mse_L2_K1","mse_per_pois2_K1","mse_ls_per_pois2",
                         "mse_full_K5","mse_U2_K5","mse_U_F2_K5","mse_lev2_K5","mse_L2_K5","mse_per_pois2_K5","mse_ls_per_pois2",
                         "mse_full_K10","mse_U2_K10","mse_U_F2_K10","mse_lev2_K10","mse_L2_K10","mse_per_pois2_K10","mse_ls_per_pois2")
result_df
#write.csv(result_df,"E:\\CQR\\s2\\p20-for-K-normal-300.csv")
