# ----
# ----  Packages

install.packages("readr")
install.packages("ggplot2")
install.packages("magrittr")
install.packages("FinTS")
install.packages("mvtnorm")
install.packages("fDMA")
install.packages("psych")
install.packages("expm")
install.packages("vrtest")
install.packages("moments")
install.packages("purrr")

# ----
# ---- Input

library(readr)
library(ggplot2)
library(magrittr)
library(FinTS)
library(mvtnorm)
library(fDMA)
library(psych)
library(expm)
library(vrtest)
library(moments)
library(purrr)

# setwd("C:/Users/User/Desktop/Thesis")
setwd("C:/Users/499332ng/OneDrive - Erasmus University Rotterdam/Desktop/Thesis")

# ----
# ---- Extension

# Input
T = 60
n = 12 # or 20
G = 2 
G2= 12
N = 99 #MC runs
S = 1000 #number of simulations

Test = c("Lee King","Engle") #the two univariate ARCH type tests

d1 = 0.8 # or 0.4 or 0.6
d2 = 0.1 # or 0.5 or 0.2
g  = 0.6 # or 0.3 or 0

m  = n/3 #number of equations
alpha = 0.05

# ----
# ---- Functions

sigma_func<-function(i,W){
  w<-W[,i]
  return(sum(w^2)/T)
}

test_procedure <- function(Test, std_res, G, sigma) {
  subsum = sum1 = sum2 = sum3 = sum4 = 0 # for LK
  if (Test == "Engle") {
    p <- archtest(std_res, lag = G)$p.value
  }
  else if (Test == "Lee King") {
    sum1   = sum(sapply((G+1):T,function(t){((std_res[t]^2/sigma)-1)*sum(sapply(1:G,function(g){std_res[t-g]^2}))}))
    sum2   = sum(sapply((G+1):T,function(t){((std_res[t]^2/sigma)-1)^2}))
    sum3   = sum(sapply((G+1):T,function(t){sum(sapply(1:G,function(g){std_res[t-g]^2}))^2}))
    sum4   = sum(sapply((G+1):T,function(t){sum(sapply(1:G,function(g){std_res[t-g]^2}))}))
    
    LK = ((T-G)*sum1/sqrt(sum2))/sqrt((T-G)*sum3-sum4^2)
    p <- pnorm(LK, mean = 0, sd = 1, lower.tail = FALSE)
  }
}

hosking <- function(std_res, G) {
  HM = 0
  Cu <- function(Z,g) {
    CuG = matrix(0,n,n)
    for (t in (g+1):T) {
      CuG = CuG + Z[t,] %*% t(Z[t-g,]) *(1/T)
    }
    return(CuG)
  }
  for (g in 1:G) {
    Cu0 <- Cu(std_res,0)
    Cug <- Cu(std_res,g)
    HM = HM + 1/(T-g) * tr(solve(Cu0) %*% Cug %*% solve(Cu0) %*% t(Cug)) 
  }
  p <- pchisq(T^2*HM, df=n^2*G, lower.tail=FALSE)
}

MC_func <- function(M,v0) {
  W_hat   = M %*% rmvt(T, sigma = diag(n), df = v0)
  W_tilda = W_hat %*% solve(chol(1/T * t(W_hat) %*% W_hat))
  D       = W_hat %*% solve(t(W_hat)%*%W_hat)%*%t(W_hat)
  return(list(W_tilda,D))
}

# ----
# Bonferroni and MC
count_LK_B = count_E_B = count_H2_B = count_LK_MC = count_E_MC = count_H2_MC = count_LK_B_G2 = count_E_B_G2 = count_H2_B_G2 = count_LK_MC_G2 = count_E_MC_G2 = count_H2_MC_G2 = 0
for (s in 1:S) {
  X = cbind(matrix(1,T,1),as.matrix(rnorm(T,0,1))) # for the #iid normal distribution
  M = diag(T) - X %*% solve(t(X) %*% X) %*% t(X)
  W = rmvnorm(T, mean = rep(0, n), sigma = diag(n))
  U = M %*% W # iid residuals
  H = matrix(1,T,n)
  
  # GJRGARCH #comment in or out when used
  # for (t in 2:T) {
  #   for (i in 1:m) {
  #     H[t,i] = 1 + d1*U[t-1,i]^2 + d2*H[t-1,i] + g*(U[t-1,i] < 0)
  #     U[t,i] = W[t,i]*sqrt(H[t,i])
  #   }
  # }
  
  # APARCH
  # for (t in 2:T) {
  #   for (i in 1:m) {
  #     H[t,i] = 1 + d1*abs(U[t-1,i]) + d2*H[t-1,i] - g*U[t-1,i]
  #     U[t,i] = W[t,i]*H[t,i]
  #   }
  # }
  
  # EGARCH
  # for (t in 2:T) {
  #   for (i in 1:m) {
  #     H[t,i] = 1 + d2*(abs(U[t-1,i]) + g*U[t-1,i])/exp(H[t-1,i]) + d1*H[t-1,i]
  #     U[t,i] = W[t,i]*sqrt(exp(H[t,i]))
  #   }
  # }
  
  # ----
  # ---- Initial value
  p_vec_LK = p_vec_E = p_vec_LK_G2 = p_vec_E_G2 = c(0)
  for (i in 1:n) {
    p_vec_LK[i] <- test_procedure(Test[1], U[,i],G,sigma_func(i,U))
    p_vec_E[i] <- test_procedure(Test[2], U[,i],G,sigma_func(i,U))
    
    p_vec_LK_G2[i] <- test_procedure(Test[1], U[,i],G2,sigma_func(i,U))
    p_vec_E_G2[i] <- test_procedure(Test[2], U[,i],G2,sigma_func(i,U))
    
  }
  p_b_LK = min(p_vec_LK)
  p_b_E  = min(p_vec_E)
  p_b_H2 = hosking(U^2, G) #hosking-ARCH
  
  p_b_LK_G2 = min(p_vec_LK_G2)
  p_b_E_G2  = min(p_vec_E_G2)
  p_b_H2_G2 = hosking(U^2, G2)
  
  count_LK_B = count_LK_B + (p_b_LK < alpha/n) #for bonferroni-type tests divide by n
  count_E_B  = count_E_B + (p_b_E < alpha/n)
  count_H2_B = count_H2_B + (p_b_H2 < alpha)
  
  count_LK_B_G2 = count_LK_B_G2 + (p_b_LK_G2 < alpha/n)
  count_E_B_G2 = count_E_B_G2 + (p_b_E_G2 < alpha/n)
  count_H2_B_G2 = count_H2_B_G2 + (p_b_H2_G2 < alpha)
  
  stat_tilda0_LK = 1 - p_b_LK
  stat_tilda0_E  = 1 - p_b_E
  stat_tilda0_H2 = 1 - p_b_H2
  
  stat_tilda0_LK_G2 = 1 - p_b_LK_G2
  stat_tilda0_E_G2  = 1 - p_b_E_G2
  stat_tilda0_H2_G2 = 1 - p_b_H2_G2
  
  # ----
  # ---- MC
  stat_tilda_LK = stat_tilda_E = stat_tilda_H2 = stat_tilda_LK_G2 = stat_tilda_E_G2 = stat_tilda_H2_G2 = c(0)
  
  for (j in 1:N) {
    # set.seed(s*1234+j)
    W       <- rmvnorm(T, mean = rep(0, n), sigma = diag(n))
    W_hat   <- M %*% W
    W_tilda <- W_hat %*% solve(chol(1/T * t(W_hat) %*% W_hat))
    
    for (i in 1:n) {
      p_vec_LK[i] <- test_procedure(Test[1], W_tilda[,i],G,sigma_func(i,W_tilda))
      p_vec_E[i]  <- test_procedure(Test[2], W_tilda[,i],G,sigma_func(i,W_tilda))
      
      p_vec_LK_G2[i] <- test_procedure(Test[1], W_tilda[,i],G2,sigma_func(i,W_tilda))
      p_vec_E_G2[i] <- test_procedure(Test[2], W_tilda[,i],G2,sigma_func(i,W_tilda))
      
    }
    stat_tilda_LK[j] = 1 - min(p_vec_LK)
    stat_tilda_E[j]  = 1 - min(p_vec_E)
    stat_tilda_H2[j] = hosking(W_tilda^2,G2)
    
    stat_tilda_LK_G2[j] = 1 - min(p_vec_LK_G2)
    stat_tilda_E_G2[j]  = 1 - min(p_vec_E_G2)
    stat_tilda_H2_G2[j] = hosking(W_tilda^2,G2)
  }
  count_LK_MC = count_LK_MC + ((sum(stat_tilda_LK>=stat_tilda0_LK) + 1)/(N+1) < alpha) #count whenever it rejects
  count_E_MC  = count_E_MC + ((sum(stat_tilda_E>=stat_tilda0_E) + 1)/(N+1) < alpha)
  count_H2_MC = count_H2_MC + ((sum(stat_tilda_H2>=stat_tilda0_H2) + 1)/(N+1) < alpha)
  
  count_LK_MC_G2 = count_LK_MC_G2 + ((sum(stat_tilda_LK_G2>=stat_tilda0_LK_G2) + 1)/(N+1) < alpha)
  count_E_MC_G2  = count_E_MC_G2 + ((sum(stat_tilda_E_G2>=stat_tilda0_E_G2) + 1)/(N+1) < alpha)
  count_H2_MC_G2 = count_H2_MC_G2 + ((sum(stat_tilda_H2_G2>=stat_tilda0_H2_G2) + 1)/(N+1) < alpha)
  print(s)
}

p_B  = rbind(cbind(count_LK_B,count_E_B,count_H2_B),cbind(count_LK_B_G2,count_E_B_G2,count_H2_B_G2))/S #number of rejections
p_MC = rbind(cbind(count_LK_MC,count_E_MC,count_H2_MC),cbind(count_LK_MC_G2,count_E_MC_G2,count_H2_MC_G2))/S
