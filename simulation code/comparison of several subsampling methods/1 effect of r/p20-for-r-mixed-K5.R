
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
r1 <- 100
r2 <- 200
r3 <- 300
r4 <- 400
r5 <- 500

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "mean0-Normal", "Mixed") 
  x <- G[[2]]
  y <- G[[1]]
  
  beta_full <- cqr.fit(X=x,y,tau5,method="ip")$beta
  mse_full <- MSE(beta_full, beta_0)
  M <- 10
  
  beta_U_r100 <- optsample(y, x, r1, 'U',tau=tau5)
  mse_U2_r100 <- MSE(beta_0, beta_U_r100)
  beta_L_r100 <- optsample(y, x, r1, 'L',tau=tau5)
  mse_L2_r100 <- MSE(beta_0, beta_L_r100)
  beta_per_pois_r100 <- sumwcqrfit(N,p,y, x, M, r1,method='pois',tau=tau5)
  mse_per_pois2_r100 <- MSE(beta_0, beta_per_pois_r100)
  beta_U_F_r100 <- Unisample(y, x, r1,tau5)
  mse_U_F2_r100 <- MSE(beta_0, beta_U_F_r100)
  beta_lev_r100 <- lev_subsample(N,p,y, x, r1,tau5)###
  mse_lev2_r100 <- MSE(beta_0, beta_lev_r100)
  beta_ls_per_pois_r100 <- sumwlsfit(N,p,y, x, M, r1)
  mse_ls_per_pois2_r100 <- MSE(beta_0, beta_ls_per_pois_r100)
  
  beta_U_r200 <- optsample(y, x, r2, 'U',tau=tau5)
  mse_U2_r200 <- MSE(beta_0, beta_U_r200)
  beta_L_r200 <- optsample(y, x, r2, 'L',tau=tau5)
  mse_L2_r200 <- MSE(beta_0, beta_L_r200)
  beta_per_pois_r200 <- sumwcqrfit(N,p,y, x, M, r2,method='pois',tau=tau5)
  mse_per_pois2_r200 <- MSE(beta_0, beta_per_pois_r200)
  beta_U_F_r200 <- Unisample(y, x, r2,tau5)
  mse_U_F2_r200 <- MSE(beta_0, beta_U_F_r200)
  beta_lev_r200 <- lev_subsample(N,p,y, x, r2,tau5)###
  mse_lev2_r200 <- MSE(beta_0, beta_lev_r200)
  beta_ls_per_pois_r200 <- sumwlsfit(N,p,y, x, M, r2)
  mse_ls_per_pois2_r200 <- MSE(beta_0, beta_ls_per_pois_r200)
  
  beta_U_r300 <- optsample(y, x, r3, 'U',tau=tau5)
  mse_U2_r300 <- MSE(beta_0, beta_U_r300)
  beta_L_r300 <- optsample(y, x, r3, 'L',tau=tau5)
  mse_L2_r300 <- MSE(beta_0, beta_L_r300)
  beta_per_pois_r300 <- sumwcqrfit(N,p,y, x, M, r3,method='pois',tau=tau5)
  mse_per_pois2_r300 <- MSE(beta_0, beta_per_pois_r300)
  beta_U_F_r300 <- Unisample(y, x, r3,tau5)
  mse_U_F2_r300 <- MSE(beta_0, beta_U_F_r300)
  beta_lev_r300 <- lev_subsample(N,p,y, x, r3,tau5)###
  mse_lev2_r300 <- MSE(beta_0, beta_lev_r300)
  beta_ls_per_pois_r300 <- sumwlsfit(N,p,y, x, M, r3)
  mse_ls_per_pois2_r300 <- MSE(beta_0, beta_ls_per_pois_r300)
  
  beta_U_r400 <- optsample(y, x, r4, 'U',tau=tau5)
  mse_U2_r400 <- MSE(beta_0, beta_U_r400)
  beta_L_r400 <- optsample(y, x, r4, 'L',tau=tau5)
  mse_L2_r400 <- MSE(beta_0, beta_L_r400)
  beta_per_pois_r400 <- sumwcqrfit(N,p,y, x, M, r4,method='pois',tau=tau5)
  mse_per_pois2_r400 <- MSE(beta_0, beta_per_pois_r400)#
  beta_U_F_r400 <- Unisample(y, x, r4,tau5)
  mse_U_F2_r400 <- MSE(beta_0, beta_U_F_r400)
  beta_lev_r400 <- lev_subsample(N,p,y, x, r4,tau5)###
  mse_lev2_r400 <- MSE(beta_0, beta_lev_r400)
  beta_ls_per_pois_r400 <- sumwlsfit(N,p,y, x, M, r4)
  mse_ls_per_pois2_r400 <- MSE(beta_0, beta_ls_per_pois_r400)
  
  beta_U_r500 <- optsample(y, x, r5, 'U',tau=tau5)
  mse_U2_r500 <- MSE(beta_0, beta_U_r500)
  beta_L_r500 <- optsample(y, x, r5, 'L',tau=tau5)
  mse_L2_r500 <- MSE(beta_0, beta_L_r500)
  beta_per_pois_r500 <- sumwcqrfit(N,p,y, x, M, r5,method='pois',tau=tau5)
  mse_per_pois2_r500 <- MSE(beta_0, beta_per_pois_r500)
  beta_U_F_r500 <- Unisample(y, x, r5,tau5)
  mse_U_F2_r500 <- MSE(beta_0, beta_U_F_r500)
  beta_lev_r500 <- lev_subsample(N,p,y, x, r5,tau5)###
  mse_lev2_r500 <- MSE(beta_0, beta_lev_r500)
  beta_ls_per_pois_r500 <- sumwlsfit(N,p,y, x, M, r5)
  mse_ls_per_pois2_r500 <- MSE(beta_0, beta_ls_per_pois_r500)
  
  
    c(mse_full,mse_U2_r100, mse_U_F2_r100,mse_lev2_r100,mse_L2_r100,mse_per_pois2_r100,mse_ls_per_pois2_r100,
      mse_full,mse_U2_r200, mse_U_F2_r200,mse_lev2_r200,mse_L2_r200,mse_per_pois2_r200,mse_ls_per_pois2_r200,
      mse_full,mse_U2_r300, mse_U_F2_r300,mse_lev2_r300,mse_L2_r300,mse_per_pois2_r300,mse_ls_per_pois2_r300,
      mse_full,mse_U2_r400, mse_U_F2_r400,mse_lev2_r400,mse_L2_r400,mse_per_pois2_r400,mse_ls_per_pois2_r400,
      mse_full,mse_U2_r500, mse_U_F2_r500,mse_lev2_r500,mse_L2_r500,mse_per_pois2_r500,mse_ls_per_pois2_r500) 
} 

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("mse_full","mse_U2_r100"," mse_U_F2_r100","mse_lev2_r100","mse_L2_r100","mse_per_pois2_r100","mse_ls_per_pois2_r100",
                         "mse_full","mse_U2_r200"," mse_U_F2_r200","mse_lev2_r200","mse_L2_r200","mse_per_pois2_r200","mse_ls_per_pois2_r200",
                         "mse_full","mse_U2_r300"," mse_U_F2_r300","mse_lev2_r300","mse_L2_r300","mse_per_pois2_r300","mse_ls_per_pois2_r300",
                         "mse_full","mse_U2_r400"," mse_U_F2_r400","mse_lev2_r400","mse_L2_r400","mse_per_pois2_r400","mse_ls_per_pois2_r400",
                         "mse_full","mse_U2_r500"," mse_U_F2_r500","mse_lev2_r500","mse_L2_r500","mse_per_pois2_r500","mse_ls_per_pois2_r500")
result_df
#write.csv(result_df,"E:\\CQR\\s2\\p20-for-r-mixed-K5.csv")


