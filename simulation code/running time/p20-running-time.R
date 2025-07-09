
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
tau1=1/2
tau5=c(1/6, 2/6, 3/6, 4/6, 5/6)
tau10=c(1/11,2/11,3/11,4/11,5/11,6/11,7/11,8/11,9/11,10/11)
beta_0 <- c(0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2,0.8,1.2)

n <- N <- 10000
p <- length(beta_0)
r1 <- 100
r2 <- 200
r3 <- 300
H <- 200    

result <- foreach(h = 1:H, .packages = c("mvtnorm", "Rfit")) %dopar% {
  G <- generate(beta_0, N, "mean0-Normal", "t3") 
  x <- G[[2]]
  y <- G[[1]]

  t_K1_M1_r100 <- run_time(N,p,y, x, M=1, r1,method='pois',tau=tau1)
  t_K1_M1_r200 <- run_time(N,p,y, x, M=1, r2,method='pois',tau=tau1)
  t_K1_M1_r300 <- run_time(N,p,y, x, M=1, r3,method='pois',tau=tau1)
  t_K1_M10_r100 <- run_time(N,p,y, x, M=10, r1,method='pois',tau=tau1)
  t_K1_M10_r200 <- run_time(N,p,y, x, M=10, r2,method='pois',tau=tau1)
  t_K1_M10_r300 <- run_time(N,p,y, x, M=10, r3,method='pois',tau=tau1)
  
  t_K5_M1_r100 <- run_time(N,p,y, x, M=1, r1,method='pois',tau=tau5)
  t_K5_M1_r200 <- run_time(N,p,y, x, M=1, r2,method='pois',tau=tau5)
  t_K5_M1_r300 <- run_time(N,p,y, x, M=1, r3,method='pois',tau=tau5)
  t_K5_M10_r100 <- run_time(N,p,y, x, M=10, r1,method='pois',tau=tau5)
  t_K5_M10_r200 <- run_time(N,p,y, x, M=10, r2,method='pois',tau=tau5)
  t_K5_M10_r300 <- run_time(N,p,y, x, M=10, r3,method='pois',tau=tau5)
  
  c(t_K1_M1_r100, t_K1_M1_r200, t_K1_M1_r300, t_K1_M10_r100, t_K1_M10_r200, t_K1_M10_r300,
    t_K5_M1_r100, t_K5_M1_r200, t_K5_M1_r300, t_K5_M10_r100, t_K5_M10_r200, t_K5_M10_r300) 
} 

stopCluster(cl)

resultM <- colMeans(matrix(unlist(result), nrow = H, byrow = TRUE))

result_df <- data.frame(t(resultM))
colnames(result_df) <- c("t_K1_M1_r100", "t_K1_M1_r200", "t_K1_M1_r300", "t_K1_M10_r100", "t_K1_M10_r200", "t_K1_M10_r300",
                         "t_K5_M1_r100", "t_K5_M1_r200", "t_K5_M1_r300", "t_K5_M10_r100", "t_K5_M10_r200", "t_K5_M10_r300")
result_df
#write.csv(result_df,"E:\\CQR\\s2\\p20-running-time.csv")
