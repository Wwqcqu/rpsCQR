rm(list=ls())
# 加载必要的包
library(quantreg)
library(lpSolve)
library(doParallel)
library(mvtnorm)
library(ggpubr)
library(quantreg)
library(lava)
#install.packages("cqr")
#install.packages("cqrReg")
library(cqrReg)
#线性规划求解样本加权复合分位数
wcqr_lp <- function (X, y, taus, weights = NULL) 
{
  n <- nrow(X)   #行数量，即样本量
  p <- ncol(X)   #列数量
  K <- length(taus) #分位数水平个数
  # 如果未提供权重，默认权重为 1
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  A1=kronecker(diag(K), rep(1, n))
  A2=kronecker(diag(K), rep(-1, n))
  A3=kronecker(cbind(rep(1, K), rep(-1, K)), X)
  A4=diag(n*K)
  A5=-diag(n*K)
  A=cbind(A1,A2,A3,A4,A5)#约束系数矩阵
  b=kronecker(rep(1, K),y)#约束右边的向量值
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
  c=cbind(c1,c2,c3)#目标函数系数向量
  const.dir=rep("=", K*n) #N个等式约束
  linprog <- lp("min", c, A, const.dir, b)
  b_tau <- linprog$sol[1:K]-linprog$sol[(K+1):(2*K)]
  beta <- linprog$sol[(2*K+1):(2*K+p)] - linprog$sol[(2*K+p+1):(2*K+p+p)] #
  return(list(intercepts = b_tau, coefficients = beta))
}


#####generate the random sample
####beta is the parameter, N is the size of sample,Xtype is the type of distribution of X,etype is the type of distribution of error 
generate <- function(beta,N,Xtype,etype){
  p <- length(beta)
  mean <- rep(0,p)
  sigma <- matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i,j] <- 0.5^(abs(i-j)) 
    }
  }
  if(Xtype == "mean0-Normal"){
    x <- mvtnorm::rmvnorm(N,mean,sigma)
  }
  else if(Xtype == "mean2-Normal"){
    x <- mvtnorm::rmvnorm(N,mean=rep(2,p),sigma)
  }
  else if(Xtype == "t3"){
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
  return(beta)
}

####Algorithm 3: perturbation subsampling
sumwcqrfit <- function(n,p,y, X, M, r,method,tau) {
  tau=tau
  q <- r/n
  K=length(tau)
  betas <- matrix(0, nrow = p, ncol = M)
  
  for (l in 1:M) {
    # Subset
    u <- rbinom(n, 1, q)
    
    # Stochastic weighting             
    if(method == 'exp'){  
      v <- rexp(n, q)    
    }else if(method == 'geom'){
      v <- rgeom(n, q)
    }else if(method == 'uni'){
      v <- runif(n,0,2/q)
    }else if(method == 'pois'){
      v <- rpois(n, 1/q)
    }else if(method == 'neg_bi'){
      v <- rnbinom(n, size = 1/2/q, prob = 1/2)
    }else if(method == 'beta'){
      v <- 3/q*rbeta(n, shape1 = 1, shape2 = 2)
    }
    else if(method == 'half_norm'){
      v <- abs(rnorm(n, 0, sqrt(pi / 2)/q ))
    }
    
    # Weights
    W <- u * v
    
    # Filter non-zero weights
    non_zero_idx <- which(W != 0)
    y_per <- y[non_zero_idx]
    X_per <- X[non_zero_idx, ]
    W_per <- W[non_zero_idx]
    
    # Estimation
    rq_model <- wcqr_lp(X_per, y_per, taus=tau, weights=W_per)
    betas[, l] <- rq_model$coefficients
  }
  
  # Combination
  beta_est <- rowMeans(betas)
  
  return(beta_est)
}



#leverage score subsampling strategy
#lev_subsample(n,p,y, X=x, r,tau)
lev_subsample <- function(n,p,y, X, r,tau) {
  tau=tau
  q <- r/n
  K=length(tau)
  h=rep(0,n)
  #betas <- matrix(0, nrow = p, ncol = m)
  for(i in 1:n)
  {
    h[i]=X[i,]%*%solve(t(X)%*%X)%*%matrix(X[i,], nrow = p, ncol = 1, byrow = TRUE)
  }
  prob=h/sum(h)
  idx <- sample(1:n, r, replace = TRUE, prob = prob)
  y.opt <- y[idx]
  x.opt <- X[idx, ]
  prob.opt <- prob[idx]
  #w=1/(r*sqrt(prob.opt))
  beta <- cqr.fit(x.opt, y.opt, tau, method="ip")$beta
  return(beta)
}




#2024Yao-LSE
sumwlsfit <- function(n,p,y, X, M, r) {
  q <- r/n
  betas <- matrix(0, nrow = p, ncol = M)
  
  for (l in 1:M) {
    # Subset
    u <- rbinom(n, 1, q)
    v <- rpois(n, 1/q)   
    # Weights
    W <- u * v
    
    # Filter non-zero weights
    non_zero_idx <- which(W != 0)
    y_per <- y[non_zero_idx]
    X_per <- X[non_zero_idx, ]
    W_per <- W[non_zero_idx]
    
    # Estimation
    model_wls <- lm(y_per ~ 0+X_per, weights = W_per) 
    betas[, l] <- model_wls$coefficients
  }
  
  # Combination
  beta_est <- rowMeans(betas)
  
  return(beta_est)
}

run_time <- function(n,p,y, X, M, r,method,tau) {
  t1 <- proc.time()
  beta <- sumwcqrfit(n,p,y, X, M, r,method,tau)
  t2 <- proc.time()
  t <- t2 - t1
  return(t[1])
}
