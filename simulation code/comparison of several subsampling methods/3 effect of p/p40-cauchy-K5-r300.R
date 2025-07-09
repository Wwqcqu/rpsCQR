
rm(list=ls())
library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)

cores <- 48
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  source("cqr-lp.R")
})
#tau=1/2
tau=c(1/6, 2/6, 3/6, 4/6, 5/6)
#tau=c(1/11,2/11,3/11,4/11,5/11,6/11,7/11,8/11,9/11,10/11)
beta_0 <- c(0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,
            0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2)

n <- N <- 10000
p <- length(beta_0)

H <- 200      ##
r <- 300 

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "mean0-Normal", "cauchy") 
  x <- G[[2]]
  y <- G[[1]]
  
  beta_full <- cqr.fit(X=x,y,tau,method="ip")$beta
  mse_full <- MSE(beta_full, beta_0)
  M <- 10
  
  beta_U <- optsample(y, x, r, 'U',tau=tau)
  mse_U2 <- MSE(beta_0, beta_U)
  
  beta_L <- optsample(y, x, r, 'L',tau=tau)
  mse_L2 <- MSE(beta_0, beta_L)
  
  beta_per_pois <- sumwcqrfit(N,p,y, x, M, r,method='pois',tau=tau)
  mse_per_pois2 <- MSE(beta_0, beta_per_pois)#20秒
  
  beta_U_F <- Unisample(y, x, r,tau)
  mse_U_F2 <- MSE(beta_0, beta_U_F)
  
  beta_lev <- lev_subsample(N,p,y, x, r,tau)###
  mse_lev2 <- MSE(beta_0, beta_lev)
  
  beta_ls_per_pois <- sumwlsfit(N,p,y, x, M, r)
  mse_ls_per_pois2 <- MSE(beta_0, beta_ls_per_pois)
  
  #beta_rr_per_pois <- sumwrfit(N,p,y, x, M, r)
  #mse_rr_per_pois2 <- MSE(beta_0, beta_rr_per_pois)
  
    c(mse_full,mse_U2, mse_U_F2,mse_lev2,mse_L2,mse_per_pois2,mse_ls_per_pois2) 
} 

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("mse_full","mse_U2","mse_U_F2","mse_lev2","mse_L2","mse_per_pois2","mse_ls_per_pois2")
result_df
#write.csv(result_df,"E:\\CQR\\s2\\new-p40-cauchy-K5-r300.csv")


