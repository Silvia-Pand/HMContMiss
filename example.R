rm(list=ls())
require(LMest)
require(mvtnorm)
require(MultiLCIRT)
source('lmbasic.cont.MISS.R')
source('lmcovlatent.cont.MISS.R')
source('lk_comp_cont_MISS.R')
source("lk_obs_latent_cont_MISS.R")
source("lk_comp_latent_cont_MISS.R")
source("prob_post_cov_cont.R")
source("est_multilogit.R")
source("prob_multilogit.R")
source("bootstrap.MISS.R")

#### load simulated data ####
load("example_data.RData")

##YYs: array containing the continuous outcomes 
#with intermittent missing (coded as NA) and dropout (coded as 999)

##XXs: array containing individual covariates

dim(YYs)
n <- dim(YYs)[1]                   # number of individuals
TT <- dim(YYs)[2]                  # number of time occasions
r <- dim(YYs)[3]                #number of continuous outcomes            
head(YYs[,1,])
dim(XXs)
head(XXs[,1,])

#### Basic HM model with a given number of states ####
est = lmbasic.cont.MISS(YYs,k=3,modBasic=1)

#individual covariates affecting the initial probabilities
X1 <- XXs[,1,]
#individual covariates affecting the transition probabilities
X2 <- XXs[,2:TT,]

#### Estimate  HM models with covariates ####
# searching for the optimal number of states
# the initialization of the EM algorithm is based on
# both a deterministic and a random rule

estcov <- lmcovlatent.cont.MISS(Y=YYs,X1=X1,X2=X2,k=2,
                                start=0, 
                                output=TRUE,
                                out_se=TRUE)

nrep <-2
Kmax<-4
modv <- vector("list",1)
bicv <- numeric(Kmax)
for(k in 1:Kmax){
  print(k)
  modv[[k]] <- lmcovlatent.cont.MISS(Y=YYs,X1=X1,X2=X2,k=k,
                                    start=0, 
                                    output=TRUE,
                                    out_se=TRUE)
  lktrace <- modv[[k]]$lk
  if(k>1){
    for(k1 in 1:(nrep*(k-1))){
      print(c(k,k1))
      tmp <- lmcovlatent.cont.MISS(Y=YYs,k=k,
                                  X1=X1,X2=X2,
                                  start=1, 
                                  output=TRUE,
                                  out_se=TRUE)
      lktrace <-c(lktrace, tmp$lk)
      if(tmp$lk>modv[[k]]$lk){
        modv[[k]] <- tmp
        print("change")}
    }
  }
  bicv[k] <- modv[[k]]$bic
}

k <- which.min(bicv)

#### Apply non parametric bootstrap with B bootrstrap samples ####
outb <- bootstrap.MISS(Y = YYs,X1,X2,
                  Mu = modv[[k]]$Mu,
                  Si=modv[[k]]$Si,
                  Be = modv[[k]]$Be,
                  Ga=modv[[k]]$Ga,
                  B = 2)
