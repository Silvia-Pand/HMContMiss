rm(list=ls())
if(!("LMest"%in%installed.packages())) install.packages("LMest")
if(!("mvtnorm"%in%installed.packages())) install.packages("mvtnorm")
require(LMest)
require(mvtnorm)
source('lmbasic.cont.MISS.R')
source('lk_comp_cont_MISS.R')
source("draw.lm.basic.cont.MISS.R")


n <- 1000   #number of individuals
TT <- 5     #number of time occasions
k <- 3      #number of latent states
r <- 3      #number of response variables
nrep <- 9
#initial probabilities (including the k+1 absorbing state)
piv <- rep(1/k,k)
piv <- c(piv,0)

#transition probabilities (including the k+1 absorbing state)
Pi <- matrix(c(0.65,0.09,0.01,0.25,0.08,0.59,0.08,0.25,0.01,0.09,0.65,0.25,0,0,0,1), k+1, k+1,byrow=TRUE)
Pi <- array(Pi, c(k+1, k+1, TT))
Pi[,,1] <- 0

#Mean vectors for each state
Mu <- matrix(c(-2,-2,0,0,0,0,0,2,2), r, k)

#Variance-covariance matrix
Si <- matrix(0.5,r,r)
diag(Si) <- 1

#frequency of intermittent missing data 
pmiss <- 0.25

#draw a sample from the HM model with dropout and intermittent missing data  
sim <- draw.lm.basic.cont.MISS(piv, Pi, Mu, Si, n, pmiss)
YY <- sim$Y
Utrue <- sim$U

est <- lmbasic.cont.MISS(YY,k=k,modBasic=1)
lktrace <- est$lk

#random starting values
if(nrep > 0) for(h in 1:nrep){
  cat("random init", h)
  esth <- lmbasic.cont.MISS(YY,k=k,modBasic=1,start=1)
  lktrace <- c(lktrace ,esth$lk)
  if(esth$lk > est$lk) 	est = esth
}

print(summary(est))