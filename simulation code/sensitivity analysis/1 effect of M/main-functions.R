rm(list=ls())
# 加载必要的包
库(quantreg)
库（lpSolve）
库(doParallel)
库(mvtnorm)
库(ggpubr)
库(quantreg)
library(lava)
#install.packages("cqr")
#install.packages("cqrReg")
library(cqrReg)
wcqr_lp <- function (X, y, taus, weights = NULL) 
{
  n <- nrow(X)   
  p <- ncol(X)  
  K <- length(taus)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  A1=kronecker(diag(K), rep(1, n))
  A2=kronecker(diag(K), rep(-1, n))
  A3=kronecker(cbind(rep(1, K), rep(-1, K)), X)
  A4=diag(n*K)
  A5=-diag(n*K)
  A=cbind(A1,A2,A3,A4,A5)
  b=kronecker(rep(1, K),y)
  c1=matrix(0, nrow = 1, ncol = 2*K + 2*p)
  c2=matrix(0, K, n)
  c3=matrix(0, K, n)
  for (i in 1:K) 
  {
    for (j in 1:n)
    {
      c2[i,j]=taus[i]*weights[j]
      c3[i,j]=(1-taus[i])*weights[j]
    }
  }
  c2 <- matrix(c(t(c2)), nrow = 1)
  c3 <- matrix(c(t(c3)), nrow = 1)
  c=cbind(c1,c2,c3)
  const.dir=rep("=", K*n) 
  linprog <- lp("min", c, A, const.dir, b)
  b_tau <- linprog$sol[1:K]-linprog$sol[(K+1):(2*K)]
  beta <- linprog$sol[(2*K+1):(2*K+p)] - linprog$sol[(2*K+p+1):(2*K+p+p)] #
  return(list(intercepts = b_tau, coefficients = beta))
}


#####generate the random sample
####beta is the parameter, N is the size of sample,Xtype is the type of distribution of X,etype is the type of distribution of error 
生成 <- 函数(beta,N,Xtype,etype){
  p <- length(beta)
  mean <- rep(0,p)
  sigma <- matrix(0,p,p)
  对于 (i 在 1:p) {
    对于 (j 在 1:p) {
      sigma[i,j] <- 0.5^(abs(i-j)) 
    }
  }
  如果(Xtype == "mean0-正常"){
    x <- mvtnorm::rmvnorm(N, mean, sigma)
  }
  否则 如果(Xtype == "mean2-正常"){
    x <- mvtnorm::rmvnorm(N,mean=rep(2,p),sigma)
  }
  否则 如果(Xtype == "t3"){
    x <- mvtnorm::rmvt(N,sigma,3)
  }
  else if(Xtype == "t6"){
    x <- mvtnorm::rmvt(N,sigma,6)
  }
  
  if(etype == "Normal"){
    error <- rnorm(N) 
  }
  else if(etype == "Mixed"){
    u <- runif(N)
    error <- rep(0, N)
    for (i in 1:N){
      if(u[i]<0.5){
        error[i] <- rnorm(1,0,1)
      }else{
        error[i] <- rnorm(1,0,0.5^3)
      }
    }
    error=sqrt(6)*error
  }
  else if(etype == "cauchy"){
    error <- rcauchy(N)
  }
  else if(etype=="t3"){#大论文可以再参考2008多考虑一点
    error <- rt(N,3)
  }
  y <- x%*%beta+error
  return(list(y,x,error,sigma))
}

###Calculate the MSE
MSE <- function(x, y) {
  mse <- sum((x - y)^2)
  return(mse)
}

###vector norm
mynorm <- function(x){
  norma <- sqrt(x%*%x)
  return(norma)
}

###Calculate the norm for each row of the matrix
norm_x <- function(x){
  nm <- apply(x, 1, mynorm)
  return(nm)
}

### optimal subsampling for quantile model
optsample <- function(y, x, r, method, tau) {
  tau=tau
  N <- nrow(x)
  p <- ncol(x) 
  K=length(tau)
  
  if (method == "U") {
    prob <- rep(1/N, N)
  } else {
    r0 <- 100   #预抽
    idx0 <- sample(1:N, r0, replace = TRUE, prob = rep(1/N, N))
    x0 <- x[idx0, ]
    y0 <- y[idx0]
    
    result0 <- cqr.fit(x0, y0, tau, method="ip")
    beta_0 <- result0$beta
    b_0 <- result0$intercepts
    
    e <- abs(y - x %*% beta_0)
    
    if (method == "A") {                                 
      C <- crossprod(x) / (N * r)
      prob <- (e * norm_x(t(C %*% t(x)))) / sum(e * norm_x(t(C %*% t(x))))
    } else if (method == "L") {
      ab=0
      for (k in 1:K)
      {
        ab=ab+I((y - x %*% beta_0)<b_0[k])-tau[k]  
      }
      prob <- (abs(ab) * norm_x(x)) / sum(abs(ab ) * norm_x(x))
    }
  }
  
  idx <- sample(1:N, r, replace = TRUE, prob = prob)
  y.opt <- y[idx]
  x.opt <- x[idx, ]
  prob.opt <- prob[idx]
  w=1/prob.opt
  beta <- wcqr_lp(x.opt, y.opt, taus=tau, weights=w)$coefficients
  return(beta)
}

#Uniform subsampling estimator based on sampling   ##without replacement##
Unisample <- function(y, x, r,tau) {
  tau=tau
  N <- nrow(x)
  p <- ncol(x) 
  K=length(tau)
  prob <- rep(1/N, N)
  idx <- sample(1:N, r, replace = FALSE, prob = prob)
  y.opt <- y[idx]
  x.opt <- x[idx, ]
  beta <- cqr.fit(x.opt, y.opt, tau, method="ip")$beta
  返回(beta)
输入：}

####算法 3: 扰动抽样
sumwcqrfit <- 函数(n,p,y, X, M, r,method,tau) {
  tau=tau
  q <- r/n
  K=length(tau)
  betas <- matrix(0, nrow = p, ncol = M)
  
对于 (l 在 1:M) {
    # 子集
    u <- rbinom(n, 1, q)
    
    # 随机加权             
    如果(方法 == 'exp'){  
      v <- rexp(n, q)    
    }else if(method == 'geom'){
      v <- rgeom(n, q)
    }否则 如果(方法 == 'uni'){
      v <- runif(n,0,2/q)
    }否则 如果(方法 == '兴趣点){
      v <- rpois(n, 1/q)
    }否则 如果(方法 == '负二值化){'
      v <- rnbinom(n, size = 1/2/q, prob = 1/2)
    否则如果方法 == '测试版){'
      v <- 3/q*rbeta(n, shape1 = 1, shape2 = 2)
    输入：}
    否则 如果(方法 == '半范数)'
      v <- abs(rnorm(n, 0, sqrt(pi / 2)/q ))
    输入：}
    
    # 权重
    W ← u * v
    
    # 过滤非零权重
    non_zero_idx <- which(W != 0)
    y_per <- y[non_zero_idx]
    X_per <- X[non_zero_idx, ]
    W_per <- W[non_zero_idx]
    
    # 估计
    rq_model <- wcqr_lp(X_per, y_per, taus=tau, weights=W_per)
    betas[, l] <- rq_model$系数
  输入：}
  
  # 组合
  beta_est <- 行平均（betas）
  
  返回（贝塔估计）
输入：}



# 杠杆得分抽样策略
#lev_subsample(n,p,y, X=x, r,tau)
lev_subsample <- 函数(n,p,y, X, r,tau) {函数(n,p,y, X, r,tau) {
  tau=tau
  q <- r/n
  K = 贝塔的长度
  h=rep(0,n)0,n)
  #贝塔值 < - 矩阵(0, nrow = p, ncol = m)
#贝塔值 < - 矩阵(0, nrow = p, ncol = m)
对于(i 在1:n)
  {
    h[i]=X[i,]%*%solve(t(X)%*%X)%*%matrix(X[i,], nrow = p, ncol = 1, byrow = TRUE)1, byrow = TRUE)
  输入：}
  概率 = h / 总和(h)
  idx <- sample(1:n, r, replace = TRUE, prob = prob)
  y.opt <- y[idx]
  x.opt <- X[idx, ]
  prob.opt <- prob[idx]
  #w=1/(r*sqrt(prob.opt))#w=1/(r*sqrt(prob.opt))
  beta <- cqr.fit(x.opt, y.opt, tau, method="ip")$beta"ip")$beta
  返回(beta)
输入：}




#2024姚-伦敦政治经济学院
sumwlsfit <- 函数(n,p,y, X, M, r) {
  q <- r/n
  betas <- matrix(0, nrow = p, ncol = M)0, nrow = p, ncol = M)
  
对于 (l 在 1:M) {
    # 子集
# 子集
    u <- rbinom(n, 1, q)1, q)
    v <- rpois(n, 1/q)   1/q)   
    # 权重
# 权重
    W ← u * v
    
    # 过滤非零权重
# 过滤非零权重
    non_zero_idx <- which(W != 0)0)
    y_per <- y[non_zero_idx]
    X_per <- X[non_zero_idx, ]
    W_per <- W[non_zero_idx]
    
    # 估计# 估计
    model_wls <- lm(y_per ~ 0+X_per, weights = W_per)0+X_per, weights = W_per) 
    betas[, l] <- model_wls$coefficients
  }
  
  # 组合
# 组合
  beta_est <- rowMeans(betas) (测试版)
  
  返回(beta估计)
输入：}
run_time <- 函数(n,p,y, X, M, r,方法,tau) {
  t1 <- proc.time()

  t2 <- proc.time()
  t <- t2 - t1
  return(t[1])1])
输入：}
