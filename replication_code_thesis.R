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

setwd("C:/Users/User/Desktop/Thesis")
# setwd("C:/Users/499332ng/OneDrive - Erasmus University Rotterdam/Desktop/Thesis")

Data0 <- read_csv("FF_Data.csv")

# subset_no = 1

Data <- Data0[13:492,] #Full
# Data <- Data0[(13+(subset_no-1)*60):(72+(subset_no-1)*60),]

# ----
# ---- Input & Regression
X <- cbind(1,Data$market,Data$SMB, Data$HML)
Y <- as.matrix(Data[2:26])

N0 = N1 = N = 999 #N0 and N1 are used in the CSMMC algo, N for MC
T  = nrow(X) #observations
G  = 12 #lags
n  = 25 #equations (portfolios)
a1 = 0.025 #for CSMMC

multivar_reg <- lm(Y  ~ X[,2:4]) #MLR
U            <- multivar_reg$residuals
Sigma        <- 1/T * t(U) %*% U
W_tilda0     <- U %*% solve(chol(Sigma))
M            <- diag(T) - X %*% solve(t(X) %*% X) %*% t(X)
D            <- U %*% solve(t(U)%*%U)%*%t(U)

Test = "Ljung Box"
Test = "Variance Ratio"
Test = "Lee King"
Test = "Engle"

# ----
# ---- Initial value (Bonferroni)
p_vec <- c(0)
for (i in 1:n) {
  p_vec[i] <- test_procedure(Test, W_tilda0[,i],G,sigma_func(i,W_tilda0))
}
p_b = min(p_vec)

p_b = hosking(W_tilda0)
p_b = hosking(W_tilda0^2)

p_b


# ---- 
# ---- MC 

stat_tilda0 = 1 - p_b

stat_tilda <- c(0)
NGnE0      <- c(0)
p_vec      <- c(0)

for (j in 1:N) {
  W       <- rmvnorm(T, mean = rep(0, n), sigma = diag(n)) #random normal distribution
  W_hat   <- M %*% W
  W_tilda <- W_hat %*% solve(chol(1/T * t(W_hat) %*% W_hat))
  
  for (i in 1:n) {
    p_vec[i] <- test_procedure(Test, W_tilda[,i],G,W_tilda)
  }
  stat_tilda[j] = 1 - min(p_vec)
  # stat_tilda[j] = hosking(W_tilda)
  # stat_tilda[j] = hosking(W_tilda^2)
  
  NGnE0[j] = stat_tilda[j] - stat_tilda0
}

p_G = (sum(NGnE0>=0) + 1)/(N+1)
p_G

# ---- 
# ---- MMC

stat_tilda <- c(0)
NGnE0      <- c(0)
p_MMC      <- c(0)
p_N1 = matrix(0,N1+1,2)

for (k in 2:(T-2-n)) { #in the effective sample size
  # for (k in kappa[-1]) { #CSMMC Stage 2
  for (j in 1:N) {
    MC = MC_func(M,k)
    W_tilda = MC[[1]]
    p_vec = c(0)
    for (i in 1:n) {
      p_vec[i] <- test_procedure(Test, W_tilda[,i],G,sigma_func(i,W_tilda))
    }
    stat_tilda[j] = 1 - min(p_vec)
    # stat_tilda[j] = hosking(W_tilda) #commented out when not used in the particular loop
    # stat_tilda[j] = hosking(W_tilda^2)
    
    NGnE0[j] = stat_tilda[j] - stat_tilda0
  }
  p_MMC[k-1] = (sum(NGnE0>=0) + 1)/(N+1)
}
p_a = max(p_MMC)


# ----
# ---- CSMMC

# Stage 1

kappa = CSK = p_CSK = p_C0 = c(0)
p_C   = E_vec_N1 = E_vec_N = matrix(0,N,2) #since N1 = N
for (k in 2:(T-n-2)) {
  # k=8
  # ---- A1-A3
  ESK_KU_0 <- ESK_KU(N0,M,k)
  E_vec_0  <- ESK_KU_0[[1]]
  E_0      <- c(abs(sum(D^3)/T^2 - ESK_KU_0[[2]][1]),
                abs(sum(diag(D^2))/T - ESK_KU_0[[2]][2]))
  
  # ---- B1 - B4
  ESK_KU_N <- ESK_KU(N,M,k)
  E_vec_N  <- ESK_KU_N[[1]]
  p_C0 <- c((sum(E_vec_N[,1]>=E_0[1]) + 1)/(N+1), 
            (sum(E_vec_N[,2]>=E_0[2]) + 1)/(N+1))
  CSK0  <- 1 - min(p_C0)
  
  # ---- C1 - C5
  ESK_KU_N1 <- ESK_KU(N1,M,k)
  E_vec_N1  <- ESK_KU_N1[[1]]
  for (l in 1:N1) {
    p_C[l,]    <- c((sum(E_vec_N[,1]>=E_vec_N1[l,1]) + 1)/(N+1),
                    (sum(E_vec_N[,2]>=E_vec_N1[l,2]) + 1)/(N+1))
    CSK[l]     <- 1 - min(p_C[l,])
  }
  
  p_CSK[k-1] <- (sum(CSK>=CSK0) + 1)/(N1+1)
  
  if (p_CSK[k-1] > a1) {
    kappa <- append(kappa, k)
  }
}
kappa[-1]

# Stage 2 -> Go to MMC

# ----
# ---- Functions

sigma_func<-function(i,W){ #for the sigma_sqrd values
  w<-W[,i]
  return(sum(w^2)/T)
}

test_procedure <- function(Test, std_res, G, sigma) {
  subsum = sum1 = sum2 = sum3 = sum4 = 0 # for LK
  if (Test == "Ljung Box") {
    p <- Box.test(ts(std_res, frequency=G), lag = G, type = "Ljung-Box")$p.value
  }
  else if (Test == "Variance Ratio") {
    VR = 1 + 2*sum(sapply(1:G, function(g){(1 - g/G)*sum(sapply((g+1):T,function(t){std_res[t]*std_res[t-g]}))/
        sum(sapply(1:T,function(t){std_res[t]^2}))}))
    p  <- 2*pnorm(VR, mean = 1, sd = 2*(2*G-1)*(G-1)/(3*G), lower.tail = FALSE) #two sided test, so multiply by two
  }
  else if (Test == "Engle") {
    p <- archtest(std_res, lag = G)$p.value
  }
  else if (Test == "Lee King") { #complicated, so i split it into four sums
    sum1   = sum(sapply((G+1):T,function(t){((std_res[t]^2/sigma)-1)*sum(sapply(1:G,function(g){std_res[t-g]^2}))}))
    sum2   = sum(sapply((G+1):T,function(t){((std_res[t]^2/sigma)-1)^2}))
    sum3   = sum(sapply((G+1):T,function(t){sum(sapply(1:G,function(g){std_res[t-g]^2}))^2}))
    sum4   = sum(sapply((G+1):T,function(t){sum(sapply(1:G,function(g){std_res[t-g]^2}))}))
    
    LK = ((T-G)*sum1/sqrt(sum2))/sqrt((T-G)*sum3-sum4^2) #merge the sums here
    print(LK)
    p <- pnorm(LK, mean = 0, sd = 1, lower.tail = FALSE)
  }
}

hosking <- function(std_res) {
  HM = 0
  Cu <- function(Z,g) { #autocovariance function
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

MC_func <- function(M,v0) { #generator used in MMC
  W_hat   = M %*% rmvt(T, sigma = diag(n), df = v0)
  W_tilda = W_hat %*% solve(chol(1/T * t(W_hat) %*% W_hat))
  D       = W_hat %*% solve(t(W_hat)%*%W_hat)%*%t(W_hat)
  return(list(W_tilda,D))
}

ESK_KU <- function(N,M,k) { #for CSMMC
  MC_all = vector(mode = "list", length = N)
  SK = KU = matrix(0,N,1)
  E_vec = matrix(0,N,2)
  
  for (i in 1:N) {
    MC_all <- MC_func(M,k)
    D_i    <- MC_all[[2]]
    SK[i] <- sum(D_i^3)/T^2
    KU[i] <- sum(diag(D_i^2))/T
  }
  SK_KU_bar = c(mean(SK),mean(KU))
  E_vec[,1] <- abs(SK - SK_KU_bar[1])
  E_vec[,2] <- abs(KU - SK_KU_bar[2])
  return(list(E_vec,SK_KU_bar))
}