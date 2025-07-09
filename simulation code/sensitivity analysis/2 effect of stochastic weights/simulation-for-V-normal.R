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


tau=c(1/6, 2/6, 3/6, 4/6, 5/6)
beta_0 <- c(0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2)

N <- 10000
p <- length(beta_0)

H <- 200
r1 <- 100
r2 <- 200
r3 <- 300
r4 <- 400
r5 <- 500

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "mean0-Normal", "Normal") 
  x <- G[[2]]
  y <- G[[1]]
  
  M <- 10
  
  beta_exp_100 <- sumwcqrfit(N,p,y, x, M, r1,method='exp',tau=tau)
  mse_exp_100 <- MSE(beta_0, beta_exp_100)
  beta_geom_100 <- sumwcqrfit(N,p,y, x, M, r1,method='geom',tau=tau)
  mse_geom_100 <- MSE(beta_0, beta_geom_100)
  beta_uni_100 <- sumwcqrfit(N,p,y, x, M, r1,method='uni',tau=tau)
  mse_uni_100 <- MSE(beta_0, beta_uni_100)
  beta_pois_100 <- sumwcqrfit(N,p,y, x, M, r1,method='pois',tau=tau)
  mse_pois_100 <- MSE(beta_0, beta_pois_100)
  beta_neg_bi_100 <- sumwcqrfit(N,p,y, x, M, r1,method='neg_bi',tau=tau)
  mse_neg_bi_100 <- MSE(beta_0, beta_neg_bi_100)
  beta_beta_100 <- sumwcqrfit(N,p,y, x, M, r1,method='beta',tau=tau)
  mse_beta_100 <- MSE(beta_0, beta_beta_100)
  beta_half_norm_100 <- sumwcqrfit(N,p,y, x, M, r1,method='half_norm',tau=tau)
  mse_half_norm_100 <- MSE(beta_0, beta_half_norm_100)
  
  beta_exp_200 <- sumwcqrfit(N,p,y, x, M, r2,method='exp',tau=tau)
  mse_exp_200 <- MSE(beta_0, beta_exp_200)
  beta_geom_200 <- sumwcqrfit(N,p,y, x, M, r2,method='geom',tau=tau)
  mse_geom_200 <- MSE(beta_0, beta_geom_200)
  beta_uni_200 <- sumwcqrfit(N,p,y, x, M, r2,method='uni',tau=tau)
  mse_uni_200 <- MSE(beta_0, beta_uni_200)
  beta_pois_200 <- sumwcqrfit(N,p,y, x, M, r2,method='pois',tau=tau)
  mse_pois_200 <- MSE(beta_0, beta_pois_200)
  beta_neg_bi_200 <- sumwcqrfit(N,p,y, x, M, r2,method='neg_bi',tau=tau)
  mse_neg_bi_200 <- MSE(beta_0, beta_neg_bi_200)
  beta_beta_200 <- sumwcqrfit(N,p,y, x, M, r2,method='beta',tau=tau)
  mse_beta_200 <- MSE(beta_0, beta_beta_200)
  beta_half_norm_200 <- sumwcqrfit(N,p,y, x, M, r2,method='half_norm',tau=tau)
  mse_half_norm_200 <- MSE(beta_0, beta_half_norm_200)
  
  beta_exp_300 <- sumwcqrfit(N,p,y, x, M, r3,method='exp',tau=tau)
  mse_exp_300 <- MSE(beta_0, beta_exp_300)
  beta_geom_300 <- sumwcqrfit(N,p,y, x, M, r3,method='geom',tau=tau)
  mse_geom_300 <- MSE(beta_0, beta_geom_300)
  beta_uni_300 <- sumwcqrfit(N,p,y, x, M, r3,method='uni',tau=tau)
  mse_uni_300 <- MSE(beta_0, beta_uni_300)
  beta_pois_300 <- sumwcqrfit(N,p,y, x, M, r3,method='pois',tau=tau)
  mse_pois_300 <- MSE(beta_0, beta_pois_300)
  beta_neg_bi_300 <- sumwcqrfit(N,p,y, x, M, r3,method='neg_bi',tau=tau)
  mse_neg_bi_300 <- MSE(beta_0, beta_neg_bi_300)
  beta_beta_300 <- sumwcqrfit(N,p,y, x, M, r3,method='beta',tau=tau)
  mse_beta_300 <- MSE(beta_0, beta_beta_300)
  beta_half_norm_300 <- sumwcqrfit(N,p,y, x, M, r3,method='half_norm',tau=tau)
  mse_half_norm_300 <- MSE(beta_0, beta_half_norm_300)
  
  beta_exp_400 <- sumwcqrfit(N,p,y, x, M, r4,method='exp',tau=tau)
  mse_exp_400 <- MSE(beta_0, beta_exp_400)
  beta_geom_400 <- sumwcqrfit(N,p,y, x, M, r4,method='geom',tau=tau)
  mse_geom_400 <- MSE(beta_0, beta_geom_400)
  beta_uni_400 <- sumwcqrfit(N,p,y, x, M, r4,method='uni',tau=tau)
  mse_uni_400 <- MSE(beta_0, beta_uni_400)
  beta_pois_400 <- sumwcqrfit(N,p,y, x, M, r4,method='pois',tau=tau)
  mse_pois_400 <- MSE(beta_0, beta_pois_400)
  beta_neg_bi_400 <- sumwcqrfit(N,p,y, x, M, r4,method='neg_bi',tau=tau)
  mse_neg_bi_400 <- MSE(beta_0, beta_neg_bi_400)
  beta_beta_400 <- sumwcqrfit(N,p,y, x, M, r4,method='beta',tau=tau)
  mse_beta_400 <- MSE(beta_0, beta_beta_400)
  beta_half_norm_400 <- sumwcqrfit(N,p,y, x, M, r4,method='half_norm',tau=tau)
  mse_half_norm_400 <- MSE(beta_0, beta_half_norm_400)
  
  beta_exp_500 <- sumwcqrfit(N,p,y, x, M, r5,method='exp',tau=tau)
  mse_exp_500 <- MSE(beta_0, beta_exp_500)
  beta_geom_500 <- sumwcqrfit(N,p,y, x, M, r5,method='geom',tau=tau)
  mse_geom_500 <- MSE(beta_0, beta_geom_500)
  beta_uni_500 <- sumwcqrfit(N,p,y, x, M, r5,method='uni',tau=tau)
  mse_uni_500 <- MSE(beta_0, beta_uni_500)
  beta_pois_500 <- sumwcqrfit(N,p,y, x, M, r5,method='pois',tau=tau)
  mse_pois_500 <- MSE(beta_0, beta_pois_500)
  beta_neg_bi_500 <- sumwcqrfit(N,p,y, x, M, r5,method='neg_bi',tau=tau)
  mse_neg_bi_500 <- MSE(beta_0, beta_neg_bi_500)
  beta_beta_500 <- sumwcqrfit(N,p,y, x, M, r5,method='beta',tau=tau)
  mse_beta_500 <- MSE(beta_0, beta_beta_500)
  beta_half_norm_500 <- sumwcqrfit(N,p,y, x, M, r5,method='half_norm',tau=tau)
  mse_half_norm_500 <- MSE(beta_0, beta_half_norm_500)
  
  c(mse_exp_100,mse_uni_100,mse_geom_100,mse_neg_bi_100,mse_pois_100,mse_beta_100,mse_half_norm_100,
    mse_exp_200,mse_uni_200,mse_geom_200,mse_neg_bi_200,mse_pois_200,mse_beta_200,mse_half_norm_200,
    mse_exp_300,mse_uni_300,mse_geom_300,mse_neg_bi_300,mse_pois_300,mse_beta_300,mse_half_norm_300,
    mse_exp_400,mse_uni_400,mse_geom_400,mse_neg_bi_400,mse_pois_400,mse_beta_400,mse_half_norm_400,
    mse_exp_500,mse_uni_500,mse_geom_500,mse_neg_bi_500,mse_pois_500,mse_beta_500,mse_half_norm_500)
}

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE_exp_100","MSE_uni_100","MSE_geom_100","mse_neg_bi_100","MSE_pois_100","MSE_beta_100","mse_half_norm_100",
                         "MSE_exp_200","MSE_uni_200","MSE_geom_200","mse_neg_bi_200","MSE_pois_200","MSE_beta_200","mse_half_norm_200",
                         "MSE_exp_300","MSE_uni_300","MSE_geom_300","mse_neg_bi_300","MSE_pois_300","MSE_beta_300","mse_half_norm_300",
                         "MSE_exp_400","MSE_uni_400","MSE_geom_400","mse_neg_bi_400","MSE_pois_400","MSE_beta_400","mse_half_norm_400",
                         "MSE_exp_500","MSE_uni_500","MSE_geom_500","mse_neg_bi_500","MSE_pois_500","MSE_beta_500","mse_half_norm_500")
result_df
#write.csv(result_df,"E:\\CQR\\s1\\p20-simul-for-V_M10\\p20-simul-for-V-normal\\p20-simul-for-V-normal.csv")
