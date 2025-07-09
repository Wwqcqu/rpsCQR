rm(list=ls())
library(mvtnorm)
library(Rfit)
library(doParallel)
library(foreach)

cores <- 48
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {source("main functions.R")})

tau=c(1/6, 2/6, 3/6, 4/6, 5/6)
beta_0 <- c(0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2)
n <-N <- 10000
p <- length(beta_0)

H <- 200
r1 <- 100
r2 <- 200
r3 <- 300
r4 <- 400
r5 <- 500

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  
  G <- generate(beta_0, N, "mean0-Normal", "Mixed") 
  x <- G[[2]]
  y <- G[[1]]
  #beta_full <- wcqr_lp(X=x, y, taus=tau) 
  
  beta_1_100 <- sumwcqrfit(N,p,y, x, M=1, r1,method='pois',tau=tau)
  beta_5_100 <- sumwcqrfit(N,p,y, x, M=5, r1,method='pois',tau=tau)
  beta_10_100 <- sumwcqrfit(N,p,y, x, M=10, r1,method='pois',tau=tau)
  beta_15_100 <- sumwcqrfit(N,p,y, x, M=15, r1,method='pois',tau=tau)
  beta_20_100 <- sumwcqrfit(N,p,y, x, M=20, r1,method='pois',tau=tau)
  beta_25_100 <- sumwcqrfit(N,p,y, x, M=25, r1,method='pois',tau=tau)
  beta_30_100 <- sumwcqrfit(N,p,y, x, M=30, r1,method='pois',tau=tau)
  mse_1_100 <- MSE(beta_0, beta_1_100)
  mse_5_100 <- MSE(beta_0, beta_5_100)
  mse_10_100 <- MSE(beta_0, beta_10_100)
  mse_15_100 <- MSE(beta_0, beta_15_100)
  mse_20_100 <- MSE(beta_0, beta_20_100)
  mse_25_100 <- MSE(beta_0, beta_25_100)
  mse_30_100 <- MSE(beta_0, beta_30_100)
  
  beta_1_200 <- sumwcqrfit(N,p,y, x, M=1, r2,method='pois',tau=tau)
  beta_5_200 <- sumwcqrfit(N,p,y, x, M=5, r2,method='pois',tau=tau)
  beta_10_200 <- sumwcqrfit(N,p,y, x, M=10, r2,method='pois',tau=tau)
  beta_15_200 <- sumwcqrfit(N,p,y, x, M=15, r2,method='pois',tau=tau)
  beta_20_200 <- sumwcqrfit(N,p,y, x, M=20, r2,method='pois',tau=tau)
  beta_25_200 <- sumwcqrfit(N,p,y, x, M=25, r2,method='pois',tau=tau)
  beta_30_200 <- sumwcqrfit(N,p,y, x, M=30, r2,method='pois',tau=tau)
  mse_1_200 <- MSE(beta_0, beta_1_200)
  mse_5_200 <- MSE(beta_0, beta_5_200)
  mse_10_200 <- MSE(beta_0, beta_10_200)
  mse_15_200 <- MSE(beta_0, beta_15_200)
  mse_20_200 <- MSE(beta_0, beta_20_200)
  mse_25_200 <- MSE(beta_0, beta_25_200)
  mse_30_200 <- MSE(beta_0, beta_30_200)
  
  beta_1_300 <- sumwcqrfit(N,p,y, x, M=1, r3,method='pois',tau=tau)
  beta_5_300 <- sumwcqrfit(N,p,y, x, M=5, r3,method='pois',tau=tau)
  beta_10_300 <- sumwcqrfit(N,p,y, x, M=10, r3,method='pois',tau=tau)
  beta_15_300 <- sumwcqrfit(N,p,y, x, M=15, r3,method='pois',tau=tau)
  beta_20_300 <- sumwcqrfit(N,p,y, x, M=20, r3,method='pois',tau=tau)
  beta_25_300 <- sumwcqrfit(N,p,y, x, M=25, r3,method='pois',tau=tau)
  beta_30_300 <- sumwcqrfit(N,p,y, x, M=30, r3,method='pois',tau=tau)
  mse_1_300 <- MSE(beta_0, beta_1_300)
  mse_5_300 <- MSE(beta_0, beta_5_300)
  mse_10_300 <- MSE(beta_0, beta_10_300)
  mse_15_300 <- MSE(beta_0, beta_15_300)
  mse_20_300 <- MSE(beta_0, beta_20_300)
  mse_25_300 <- MSE(beta_0, beta_25_300)
  mse_30_300 <- MSE(beta_0, beta_30_300)
  
  beta_1_400 <- sumwcqrfit(N,p,y, x, M=1, r4,method='pois',tau=tau)
  beta_5_400 <- sumwcqrfit(N,p,y, x, M=5, r4,method='pois',tau=tau)
  beta_10_400 <- sumwcqrfit(N,p,y, x, M=10, r4,method='pois',tau=tau)
  beta_15_400 <- sumwcqrfit(N,p,y, x, M=15, r4,method='pois',tau=tau)
  beta_20_400 <- sumwcqrfit(N,p,y, x, M=20, r4,method='pois',tau=tau)
  beta_25_400 <- sumwcqrfit(N,p,y, x, M=25, r4,method='pois',tau=tau)
  beta_30_400 <- sumwcqrfit(N,p,y, x, M=30, r4,method='pois',tau=tau)
  mse_1_400 <- MSE(beta_0, beta_1_400)
  mse_5_400 <- MSE(beta_0, beta_5_400)
  mse_10_400 <- MSE(beta_0, beta_10_400)
  mse_15_400 <- MSE(beta_0, beta_15_400)
  mse_20_400 <- MSE(beta_0, beta_20_400)
  mse_25_400 <- MSE(beta_0, beta_25_400)
  mse_30_400 <- MSE(beta_0, beta_30_400)
  
  beta_1_500 <- sumwcqrfit(N,p,y, x, M=1, r5,method='pois',tau=tau)
  beta_5_500 <- sumwcqrfit(N,p,y, x, M=5, r5,method='pois',tau=tau)
  beta_10_500 <- sumwcqrfit(N,p,y, x, M=10, r5,method='pois',tau=tau)
  beta_15_500 <- sumwcqrfit(N,p,y, x, M=15, r5,method='pois',tau=tau)
  beta_20_500 <- sumwcqrfit(N,p,y, x, M=20, r5,method='pois',tau=tau)
  beta_25_500 <- sumwcqrfit(N,p,y, x, M=25, r5,method='pois',tau=tau)
  beta_30_500 <- sumwcqrfit(N,p,y, x, M=30, r5,method='pois',tau=tau)
  mse_1_500 <- MSE(beta_0, beta_1_500)
  mse_5_500 <- MSE(beta_0, beta_5_500)
  mse_10_500 <- MSE(beta_0, beta_10_500)
  mse_15_500 <- MSE(beta_0, beta_15_500)
  mse_20_500 <- MSE(beta_0, beta_20_500)
  mse_25_500 <- MSE(beta_0, beta_25_500)
  mse_30_500 <- MSE(beta_0, beta_30_500)
  
  c(mse_1_100, mse_5_100,mse_10_100,mse_15_100,mse_20_100,mse_25_100,mse_30_100,
    mse_1_200, mse_5_200,mse_10_200,mse_15_200,mse_20_200,mse_25_200,mse_30_200,
    mse_1_300, mse_5_300,mse_10_300,mse_15_300,mse_20_300,mse_25_300,mse_30_300,
    mse_1_400, mse_5_400,mse_10_400,mse_15_400,mse_20_400,mse_25_400,mse_30_400,
    mse_1_500, mse_5_500,mse_10_500,mse_15_500,mse_20_500,mse_25_500,mse_30_500)
}
stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("MSE_1_100", "MSE_5_100","MSE_10_100","MSE_15_100","MSE_20_100","MSE_25_100","MSE_30_100",
                         "MSE_1_200", "MSE_5_200","MSE_10_200","MSE_15_200","MSE_20_200","MSE_25_200","MSE_30_200",
                         "MSE_1_300", "MSE_5_300","MSE_10_300","MSE_15_300","MSE_20_300","MSE_25_300","MSE_30_300",
                         "MSE_1_400", "MSE_5_400","MSE_10_400","MSE_15_400","MSE_20_400","MSE_25_400","MSE_30_400",
                         "MSE_1_500", "MSE_5_500","MSE_10_500","MSE_15_500","MSE_20_500","MSE_25_500","MSE_30_500")
result_df
#write.csv(result_df,"E:\\CQR\\s1\\p20-simul-for-M_pois\\p20-simul-for-M-mixed\\p20-simul-for-M-mixed.csv")
