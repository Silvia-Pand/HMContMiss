lmcovlatent.cont.MISS <- function(Y,X1=NULL,X2=NULL,
                             yv = rep(1,nrow(Y)),k,start=0,
                             tol=10^-8,maxit=1000,
                             Mu=NULL,Si=NULL,Be=NULL,Ga=NULL,
                             output=FALSE, out_se=FALSE, miss = FALSE){

  # Fit the LM model for continuous outcomes with individual covariates in the distribution of the latent process
  # dealing with intermittent missingness and dropout
  # INPUT:
  # Y = array of available continuous outcome (n x TT x r)
  # X1 = matrix of covariates affecting the initial probabilities
  # X2 = array of covariates affecting the transition probabilities
  # yv = vector of frequencies
  # k = number of latent states
  # start = initialization (0 = deterministic, 1 = random, 2 = initial values in input)
  # maxit = maximum number of iterations
  # param = type of parametrization for the transition probabilities:
  #         multilogit = standard multinomial logit for every row of the transition matrix
  #         difflogit  = multinomial logit based on the difference between two sets of parameters
  # Mu = conditional means of the response variables (if start=2)
  # Si = var-cov matrix common to all states (if start=2)
  # Be = parameters on the initial probabilities (if start=2)
  # Ga = parameters on the transition probabilities (if start=2)
  # output = to return additional output

  # Preliminaries
  param="multilogit"
  check_der = FALSE # to check score and info
  fort=FALSE
  sY = dim(Y)
  n = sY[1]
  TT = sY[2]
  if(length(sY)==2) r = 1
  else r = sY[3]

 
  ## Check and inpute for missing data
  miss = any(is.na(Y))
  drop = any(Y==999,na.rm=TRUE)
  R= NULL
  if(miss){
    R = (!is.na(Y))
    if(fort) RR = array(as.integer(1*R),c(n,TT,r))
    Y[is.na(Y)] = 0
    cat("Missing data in the dataset.\n")
  }

  Yv = matrix(Y,n*TT,r)
  
  if(drop){
    Yv2 = Yv[Yv[,1]!=999,] 
  }
    # Covariate structure and related matrices: initial probabilities
  if(k == 2){
    GBe = as.matrix(c(0,1))
  }else{
    GBe = diag(k); GBe = GBe[,-1]
  }
  if(is.null(X1)){
    nc1=0
    Xlab = rep(1,n)
    nameBe = NULL
  }else{
    if(is.vector(X1)) X1 = matrix(X1,n,1)
    nc1 = dim(X1)[2] # number of covariates on the initial probabilities
    if(n!= dim(X1)[1]) stop("dimension mismatch between S and X1")
    nameBe = colnames(X1)
    out = aggr_data(X1)
    Xdis = out$data_dis
    if(nc1==1) Xdis = matrix(Xdis,length(Xdis),1)
    Xlab = out$label
  }
  Xndis = max(Xlab)
  XXdis = array(0,c(k,(k-1)*(nc1+1),Xndis))
  for(i in 1:Xndis){
    if(nc1==0) xdis = 1 else xdis = c(1,Xdis[i,])
    XXdis[,,i] = GBe%*%(diag(k-1)%x%t(xdis))
  }


  # for the transition probabilities
  if(is.null(X2)){
    nc2 = 0
    Zlab = rep(1,n*(TT-1))
    nameGa = NULL
    Zndis = max(Zlab)
  }else{
    nc2 = dim(X2)[3] # number of covariates on the transition probabilities
    if(n!= dim(X2)[1]) stop("dimension mismatch between S and X2")
    nameGa = colnames(aperm(X2,c(1,3,2)))
    Z = NULL
    for(t in 1:(TT-1)) Z = rbind(Z,X2[,t,])
    if(nc2==1) Z = as.vector(X2)
    out = aggr_data(Z); Zdis = out$data_dis; Zlab = out$label; Zndis = max(Zlab)
    if(nc2==1) Zdis=matrix(Zdis,length(Zdis),1)
  }
    if(drop) ZZdis = array(0,c(k+1,((k+1)-1)*(nc2+1),Zndis,k))
    else ZZdis = array(0,c(k,(k-1)*(nc2+1),Zndis,k))
    for(h in 1:k){
      if(drop){
        if(k+1==2){
          if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
        }else{
          GGa = diag(k+1); GGa = GGa[,-h]
        }
      }
      else{
        if(k==2){
          if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
        }else{
          GGa = diag(k); GGa = GGa[,-h]
        }
      }
      for(i in 1:Zndis){
        if(nc2==0) zdis = 1 else zdis = c(1,Zdis[i,])
        if(drop) ZZdis[,,i,h] = GGa%*%(diag((k+1)-1)%x%t(zdis))
        else ZZdis[,,i,h] = GGa%*%(diag(k-1)%x%t(zdis))
      }
    }

  # When there is just 1 latent class
  if(k == 1){
    Piv = rep(1,n); Pi = 1
    yvv = rep(yv,TT)
    Mu = colSums(yvv*Yv,na.rm=TRUE)/sum(yvv)
    Di = Yv-rep(1,n*TT)%o%Mu
    Si = t(Di)%*%(yvv*Di)/sum(yvv)
    lk = 0
    if(r==1){
      for (i in 1:n) for(t in 1:TT){
        indo = R[i,t,]
        if(sum(R[i,t,])==1) lk = lk +dnorm(Y[i,t,][indo],Mu[indo],Si[indo,indo],log=TRUE)        }
    }else{
      for (i in 1:n) for(t in 1:TT){
        indo = R[i,t,]
        if(sum(R[i,t,])==1)
        {
          lk = lk +dnorm(Y[i,t,][indo],Mu[indo],Si[indo,indo],log=TRUE)
        }else if(sum(R[i,t,])>1)
        {
          lk = lk + dmvnorm(Y[i,t,][indo],Mu[indo],Si[indo,indo],log=TRUE)
        }
      }
    }
    np = k*r+r*(r+1)/2
    aic = -2*lk+np*2
    bic = -2*lk+np*log(n)
    Mu = matrix(Mu,r,k)
    nameY <- dimnames(Y)[[3]]
    dimnames(Mu) <- list(nameY,state=1:k)
    out = list(lk=lk,Piv=Piv,Pi=Pi,Mu=Mu,Si=Si,np=np,k = k, aic=aic,bic=bic,lkv=NULL,V=NULL, n = n, TT = TT, paramLatent=param )
    class(out)="LMlatentcont"
    return(out)
  }
  # Starting values: deterministic initialization
  if(start == 0){
    if(drop){
      mu = colMeans(Yv2,na.rm=TRUE)
      Si = cov(Yv2,use = "complete.obs"); std = sqrt(diag(Si))
    }else{
      mu = colMeans(Yv,na.rm=TRUE)
      Si = cov(Yv,use = "complete.obs"); std = sqrt(diag(Si))
    }
    std = sqrt(diag(Si))
    qt = qnorm((1:k)/(k+1))
    Mu = matrix(0,r,k)
    for(u in 1:k) Mu[,u] = qt[u]*std+mu
    
    # parameters on initial probabilities
    be = array(0,(nc1+1)*(k-1))
    out = prob_multilogit(XXdis,be,Xlab)
    Piv = out$P; Pivdis = out$Pdis
    
    # parameters on transition probabilities
    if(drop){
        Ga = matrix(0,(nc2+1)*((k+1)-1),k)
        Ga[1+(0:((k+1)-2))*(nc2+1),] = -log(10)
        PIdis = array(0,c(Zndis,k+1,k+1)); PI = array(0,c(k+1,k+1,n,TT))
      }else{
        Ga = matrix(0,(nc2+1)*(k-1),k)
        Ga[1+(0:(k-2))*(nc2+1),] = -log(10)
        PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))  
      }
      
      for(h in 1:k){
        tmp = ZZdis[,,,h]
        if(nc2==0) tmp = array(tmp,c(k,(k-1),Zndis))
        out = prob_multilogit(tmp,Ga[,h],Zlab)
        PIdis[,,h] = out$Pdis
        if(drop) PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k+1,n,TT-1))
        else PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
      }
  }
  # random initialization
  if(start==1){
    Mu = matrix(0,r,k)
    if(drop){
      mu = colMeans(Yv2,na.rm=TRUE)
      Si = cov(Yv2,use = "complete.obs")
    }else{
      mu = colMeans(Yv,na.rm=TRUE)
      Si = cov(Yv,use = "complete.obs")
    }
    for(u in 1:k) Mu[,u] = rmvnorm(1,mu,Si)
    
    # parameters on initial probabilities
    be = c(rnorm(1),rep(0,nc1))
    if(k>2) for(h in 2:(k-1)) be = c(be,rnorm(1),rep(0,nc1))
    out = prob_multilogit(XXdis,be,Xlab)
    Piv = out$P; Pivdis = out$Pdis
    # parameters on transition probabilities
    if(drop){
        Ga = matrix(0,(nc2+1)*((k+1)-1),k)
        Ga[1+(0:((k+1)-2))*(nc2+1),] = -abs(rnorm(((k+1)-1)))     
        PIdis = array(0,c(Zndis,k+1,k+1)); PI = array(0,c(k+1,k+1,n,TT))
      }else{
        Ga = matrix(0,(nc2+1)*(k-1),k)
        Ga[1+(0:(k-2))*(nc2+1),] = -abs(rnorm((k-1)))
        PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
      }
      for(h in 1:k){
        tmp = ZZdis[,,,h]
        if(nc2==0) tmp = array(tmp,c(k,(k-1),Zndis))
        out = prob_multilogit(tmp,Ga[,h],Zlab)
        PIdis[,,h] = out$Pdis
        if(drop) PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k+1,n,TT-1))
        else PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
      }
  }
  # initialization as input
  if(start==2){
    # parameters on initial probabilities
    be = as.vector(Be)
    out = prob_multilogit(XXdis,be,Xlab)
    Piv = out$P; Pivdis = out$Pdis
    # parameters on transition probabilities
    if(is.list(Ga)) stop("invalid mode (list) for Ga")
      if(drop){
        Ga = matrix(Ga,(nc2+1)*((k+1)-1),k)
        PIdis = array(0,c(Zndis,k+1,k+1)); PI = array(0,c(k+1,k+1,n,TT))
      }else{
        Ga = matrix(Ga,(nc2+1)*(k-1),k)
        PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
      }  
      for(h in 1:k){
        out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab)
        PIdis[,,h] = out$Pdis 
        if(drop) PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k+1,n,TT-1))
        else PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
      }
  }
  if(drop){
    Piv = cbind(Piv,0)
    PI[k+1,k+1,,2:TT] = 1
  }
  ###### standard EM #####
  out = lk_comp_latent_cont_miss(Y,R,yv,Piv,PI,Mu,Si,k,drop=drop)
  lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  it = 0; lko = lk-10^10; lkv = NULL
  par = c(as.vector(Piv),as.vector(PI),as.vector(Mu),as.vector(Si))
  if(any(is.na(par))) par = par[-which(is.na(par))]
  paro = par
  # Iterate until convergence
  # display output
  cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
  cat("      k     |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
  cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
  cat(sprintf("%11g",c(k,start,0,lk)),"\n",sep=" | ")
  while((lk-lko)/abs(lk)>tol & it<maxit){
    Mu0 = Mu; Si0 = Si; Piv0 = Piv; PI0 = PI
    it = it+1
    # ---- E-step ----
    # Compute V and U
    out = prob_post_cov_cont(Y,yv,Mu,Si,Piv,PI,Phi,L,pv)
    U = out$U; V = out$V
    # If required store parameters
    # ---- M-step ----
    # Update Mu
    if(drop) Vv = matrix(aperm(V,c(1,3,2)),n*TT,k+1) else Vv = matrix(aperm(V,c(1,3,2)),n*TT,k)
    
    if(miss){
        if(drop) Y1 = array(Y,c(n,TT,r,k+1)) else Y1 = array(Y,c(n,TT,r,k))
        Var = array(0,c(n,TT,r,r))
        for(i in 1:n) for(t in 1:TT){
            nr = sum(R[i,t,])
            if(nr==0){
              Y1[i,t,,1:k] = Mu
              Var[i,t,,] = Si
            }else if(nr<r){
              indo = R[i,t,]; indm = !R[i,t,]
              Tmp = Si[indm,indo]%*%solve(Si[indo,indo])
              Var[i,t,indm,indm] = Si[indm, indm]-Tmp%*%Si[indo,indm]
              for(u in 1:k) Y1[i,t,indm,u] = Mu[indm,u]+Tmp%*%(Y[i,t,indo]-Mu[indo,u])
            }
          }
        Mu = matrix(0,r,k)
        for(u in 1:k){
          Yv1 = matrix(Y1[,,,u],n*TT)
          Mu[,u] = (t(Yv1)%*%Vv[,u])/sum(Vv[,u])
        }
      
        Sitmp = matrix(0,r,r)
        for(u in 1:k){
          Yv1 = matrix(Y1[,,,u],n*TT)
          Var1 = array(Var,c(n*TT,r,r))
          Tmp = Yv1-rep(1,n*TT)%*%t(Mu[,u])
          Sitmp = Sitmp+t(Tmp)%*%(Vv[,u]*Tmp)+apply(Vv[,u]*Var1,c(2,3),sum)
        }
        if(drop) Si  = Sitmp/dim(Yv2)[1] else Si = Sitmp/(n*TT) 
      }else{
      for(u in 1:k) Mu[,u] = (t(Yv)%*%Vv[,u])/sum(Vv[,u]) 
      # Update Si
      Si = matrix(0,r,r)
      for(u in 1:k) Si= Si+ t(Yv-rep(1,n*TT)%o%Mu[,u])%*%diag(Vv[,u])%*%
        as.matrix(Yv-rep(1,n*TT)%o%Mu[,u])
      #Si = Si/(sum(yv)*TT)
      if(drop) Si  = Si/dim(Yv2)[1] else Si = Si/(n*TT)
    }

    # Update piv
    out = est_multilogit(V[,1:k,1],XXdis,Xlab,be,Pivdis)
    be = out$be; Pivdis = out$Pdi; Piv = out$P
    if(drop) Piv = cbind(Piv,0)
    # Update Pi
    for(h in 1:k){
        UU = NULL
        for(t in 2:TT) UU = rbind(UU,t(U[h,,,t]))
        tmp = ZZdis[,,,h]
        if(nc2==0) tmp = array(tmp,c(k,(k-1),Zndis))
        tmp2 = PIdis[,,h]
        if(Zndis==1) tmp2 = matrix(tmp2,1,k)
        out = est_multilogit(UU,tmp,Zlab,Ga[,h],tmp2)
        PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = t(out$P); Ga[,h] = out$be
        }
       # Compute log-likelihood
    paro = par; par = c(as.vector(Piv),as.vector(PI),as.vector(Mu),as.vector(Si))
    if(any(is.na(par))) par = par[-which(is.na(par))]
    lko = lk
    out = lk_comp_latent_cont_miss(Y,R,yv,Piv,PI,Mu,Si,k,drop=drop)
    lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    # Display output
    if(it/10 == floor(it/10)){
      cat(sprintf("%11g",c(k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
    }
    lkv = c(lkv,lk)
  }
  if(it/10 > floor(it/10))  cat(sprintf("%11g",c(k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
  if(out_se){
    th = NULL
    th = c(th,as.vector(Mu))
    th = c(th,Si[upper.tri(Si,TRUE)])
    th = c(th, be)
    for(h in 1:k) th = c(th, Ga[,h])
   
    out = lk_obs_latent_cont_MISS(th,Y,R,yv,XXdis,Xlab,ZZdis,Zlab,param,drop=drop)
    lk0 = out$lk; sc0 = out$sc

    lth = length(th)
    scn = rep(0,lth); Fn = matrix(0,lth,lth)

    for(h in 1:lth){
      thh = th; thh[h] = thh[h]+10^-5
      outh = lk_obs_latent_cont_MISS(thh,Y,R,yv,XXdis,Xlab,ZZdis,Zlab,param,drop=drop)
      scn[h] = (outh$lk-lk0)/10^-5
      Fn[,h] = (outh$sc-sc0)/10^-5
    }

    J = -(Fn+t(Fn))/2
    iJ = ginv(J)
    se = sqrt(diag(iJ))

    nMu = r*k
    nSi = r*(r+1)/2
    nbe = (1+nc1)*(k-1)
      if(drop) nga=(1+nc2)*((k+1)-1)*k else nga=(1+nc2)*(k-1)*k
    seMu = se[1:nMu]
    seSi = se[nMu+(1:nSi)]
    sebe = se[nMu+nSi+(1:nbe)]
    sega = se[nMu+nSi+nbe+(1:nga)]
  }

  # Compute number of parameters
  np = k*r+r*(r+1)/2 #Mu e Si
  np = np+(k-1)*(nc1+1) #Be
  if(drop) np = np+((k+1)-1)*(nc2+1)*k 
  else np = np+(k-1)*(nc2+1)*k 
  
  aic = -2*lk+np*2
  bic = -2*lk+np*log(n)
  # local decoding
  Ul = matrix(0,n,TT)
  for(i in 1:n) for(t in 1:TT){
    Ul[i,t] = which.max(V[i,,t])
  }

  Be = matrix(be,nc1+1,k-1)
  if (is.null(nameBe)){
    if(nc1==0) nameBe = c("Intercept") else nameBe = c("intercept",paste("X1",1:nc1,sep=""))
  }else{
    nameBe = c("intercept",nameBe)
  }

  dimnames(Be) = list(nameBe,logit=2:k)
  if(out_se) {seBe = matrix(sebe,nc1+1,k-1); dimnames(seBe) = list(nameBe,logit=2:k)}

    if(is.null(nameGa)){
      if(nc2==0) nameGa = c("Intercept") else nameGa = c("intercept", paste("X2",1:nc2,sep=""))
    }else{
      nameGa = c("intercept",nameGa)
    }
    if(k>2) {
      if(drop){
        Ga = array(as.vector(Ga),c(nc2+1,(k+1)-1,k))
        dimnames(Ga) = list(nameGa,logit=2:(k+1),logit=1:k)
      } else{
        Ga = array(as.vector(Ga),c(nc2+1,k-1,k))
        dimnames(Ga) = list(nameGa,logit=2:k,logit=1:k)
      } 
    }else if(k==2){
      if(drop){
        Ga = array(as.vector(Ga),c(nc2+1,(k+1)-1,k))
        dimnames(Ga) = list(nameGa,logit=2:(k+1),logit=1:k)
      } else{
      dimnames(Ga) = 	list(nameGa,logit=1:k)
      }
    }

    if(out_se){
      if(k==2){
        if(drop){
          seGa = array(as.vector(sega),c(nc2+1,(k+1)-1,k))
          dimnames(seGa) = list(nameGa,logit=2:(k+1),logit=1:k)
        } else{
          seGa = matrix(sega,nc2+1,2)
          dimnames(seGa) = list(nameGa,logit=1:k)
        }
      }else if(k>2){
        if(drop){
          seGa = array(as.vector(sega),c(nc2+1,(k+1)-1,k))
          dimnames(Ga) = list(nameGa,logit=2:(k+1),logit=1:k)
        } else{
          seGa = array(as.vector(sega),c(nc2+1,k-1,k))
          dimnames(seGa) = list(nameGa,logit=2:k,logit=1:k)
        }
      }
    }

  # adjust output
  lk = as.vector(lk)
  if(output){
    if(drop) {
      dimnames(Piv)=list(subject=1:n,state=1:(k+1))
      dimnames(PI)=list(state=1:(k+1),state=1:(k+1),subject=1:n,time=1:TT)
    }else{
      dimnames(Piv)=list(subject=1:n,state=1:k)
      dimnames(PI)=list(state=1:k,state=1:k,subject=1:n,time=1:TT)
    }
  }
  nameY <- dimnames(Y)[[3]]
  dimnames(Mu) <- list(nameY,state=1:k)
  if(r==1) colnames(Si) <- nameY else dimnames(Si) <- list(nameY,nameY)
  
  out = list(lk=lk,Be=Be,Ga=Ga,Mu=Mu,Si=Si,np=np,k = k,aic=aic,bic=bic,lkv=lkv, n = n, TT = TT,paramLatent=param )
  if(out_se){
    seMu = matrix(seMu,r,k)
    seSi2 = matrix(0,r,r)
    seSi2[upper.tri(seSi2,TRUE)]=seSi
    if(r>1) seSi2 = seSi2+t(seSi2-diag(diag(seSi2)))
    seSi = seSi2
    dimnames(seMu)=list(nameY,state=1:k)
    if(r==1) colnames(seSi) <- nameY else dimnames(seSi) <- list(nameY,nameY)
    
    out$seMu = seMu
    out$seSi = seSi
    out$seBe = seBe
    out$seGa = seGa
  }
  # final output
  if(miss) out$Y = Y
  if(output){
    out$PI = PI
    out$Piv = Piv
    out$Ul = Ul
  }
  cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
  class(out)="LMlatentcont"
  return(out)
}
