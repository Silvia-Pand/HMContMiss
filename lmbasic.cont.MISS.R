lmbasic.cont.MISS  <- function(Y,k,start=0,modBasic=0,tol=10^-8,maxit=1000,out_se=FALSE,piv=NULL,Pi=NULL,Mu=NULL,Si=NULL,indb=NULL){


    # Preliminaries
    check_der = TRUE  # to check derivatives
    sY = dim(Y)
    n = as.integer(sY[1])
    k = as.integer(k)
    TT = as.integer(sY[2])
    mod <- modBasic
    if(is.data.frame(Y)){
      warning("Data frame not allowed for Y")
    }
    if(length(sY)==2){
      r = 1
      if(is.matrix(Y)) Y = array(Y,c(dim(Y),1))
    }else r = sY[3]
    r = as.integer(r)

    ## Check and inpute for missing data
    miss = any(is.na(Y))
    drop = any(Y==999,na.rm=TRUE)
    R= NULL
    if(miss){
      R = (!is.na(Y))
      Y[is.na(Y)] = 0
      cat("Missing data in the dataset.\n")
    }

    Yv = matrix(Y,n*TT,r)
    if(drop) Yv2 = Yv[Yv[,1]!=999,]
    th = NULL; sc = NULL; J = NULL
    if(out_se){
      B = cbind(-rep(1,k-1),diag(k-1))
      Bm = rbind(rep(0,k-1),diag(k-1))
      C = array(0,c(k-1,k,k))
      Cm = array(0,c(k,k-1,k))
      for(u in 1:k){
        C[,,u] = rbind(cbind(diag(u-1),-rep(1,u-1),matrix(0,u-1,k-u)),
                       cbind(matrix(0,k-u,u-1),-rep(1,k-u),diag(k-u)))
        Cm[,,u] = rbind(cbind(diag(u-1),matrix(0,u-1,k-u)),
                        rep(0,k-1),
                        cbind(matrix(0,k-u,u-1),diag(k-u)))
      }
    }

    # When there is just 1 latent class
    if(k == 1){
      piv = 1; Pi = 1
      Mu = as.matrix(colMeans(Yv,na.rm=TRUE))
      Si = as.matrix(cov(Yv,use = "complete.obs"))
      rownames(Mu) <- dimnames(Y)[[3]]
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
      ent = NULL
      clc = NULL
      icl.bic = NULL
      rownames(Mu) <- dimnames(Y)[[3]]
      out = list(lk=lk,piv=piv,Pi=Pi,Mu=Mu,Si=Si,np=np, k = k, ent = ent, clc = clc, 
                 icl.bic = icl.bic,
                 aic=aic,bic=bic,lkv=NULL,V=NULL,n = n, TT = TT, modBasic = mod )
      class(out)="LMbasiccont"
      return(out)
    }
    # Starting values
    if(start == 0){
       if(drop){
         mu = colMeans(Yv2,na.rm=TRUE)
         Si = cov(Yv2,use = "complete.obs"); std = sqrt(diag(Si))
       }else{
        mu = colMeans(Yv,na.rm=TRUE)
        Si = cov(Yv,use = "complete.obs"); std = sqrt(diag(Si))
      }
      qt = qnorm((1:k)/(k+1))
      Mu = matrix(0,r,k)
      for(u in 1:k) Mu[,u] = qt[u]*std+mu
      if(length(indb)>0) Mu[indb,] = colMeans(Yv)[indb]
      piv = rep(1,k)/k
      if(drop){
        Pi = matrix(1,k+1,k+1)+9*diag(k+1); Pi = diag(1/rowSums(Pi))%*%Pi;
        Pi = array(Pi,c(k+1,k+1,TT)); Pi[,,1] = 0
      }else{
        Pi = matrix(1,k,k)+9*diag(k); Pi = diag(1/rowSums(Pi))%*%Pi;
        Pi = array(Pi,c(k,k,TT)); Pi[,,1] = 0
      }
    }
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
      if(drop) Pi = array(runif((k+1)^2*TT),c(k+1,k+1,TT)) else Pi = array(runif(k^2*TT),c(k,k,TT))
      for(t in 2:TT) Pi[,,t] = diag(1/rowSums(Pi[,,t]))%*%Pi[,,t]
      Pi[,,1] = 0
      piv = runif(k); piv = piv/sum(piv)
    }
    if(start==2){
      if(is.null(piv)) stop("initial value of the initial probabilities (piv) must be given in input")
      if(is.null(Pi)) stop("initial value of the transition probabilities (Pi) must be given in input")
      if(is.null(Mu)) stop("initial value of the conditional means of the response variables (Mu) must be given in input")
      if(is.null(Si)) stop("initial value of the var-cov matrix common to all states (Si) must be given in input")
      piv = piv
      Pi = Pi
      Mu = Mu
      Si = Si
    }
    if(drop){
      piv = c(piv,0)
      Pi[k+1,1:k,]=0 
      Pi[k+1,k+1,2:TT] = 1
    }
        # Compute log-likelihood
    out = lk_comp_cont_MISS(Y,R,piv,Pi,Mu,Si,k,drop=drop)
    lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("     mod    |      k      |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat(sprintf("%11g",c(mod,k,start,0,lk)),"\n",sep = " | ")
    it = 0; lko = lk-10^10; lkv = NULL
    par = c(piv,as.vector(Pi),as.vector(Mu),as.vector(Si))
    if(any(is.na(par))) par = par[-which(is.na(par))]
    paro = par
    # Iterate until convergence
    while((lk-lko)/abs(lk)>tol & it<maxit){
      # t0 = proc.time()
      Mu0 = Mu; Si0 = Si; piv0 = piv; Pi0 = Pi
      it = it+1;
      # ---- E-step ----
      # Compute V and U
      if(drop){
        V = array(0,c(n,k+1,TT)); U = array(0,c(k+1,k+1,TT))
        M = matrix(1,n,k+1)
      }else{
        V = array(0,c(n,k,TT)); U = array(0,c(k,k,TT))
        M = matrix(1,n,k)
      }  
      if(n==1) V[,,TT] = L[,,TT]/sum(L[1,,TT])
      else V[,,TT] = L[,,TT]/rowSums(L[,,TT])
      for(i in 1:n){
          Tmp = (L[i,,TT-1]%o%Phi[i,,TT])*Pi[,,TT]
          U[,,TT] = U[,,TT]+Tmp/sum(Tmp)
        }
      if(TT>2){
        for(t in seq(TT-1,2,-1)){
          M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1])
          M = M/rowSums(M)
          V[,,t] = L[,,t]*M
          if(n==1) V[,,t] = V[,,t]/sum(V[1,,t])
          else V[,,t] = V[,,t]/rowSums(V[,,t])
            for(i in 1:n){
              Tmp = (L[i,,t-1]%o%(Phi[i,,t]*M[i,]))*Pi[,,t]
              U[,,t] = U[,,t]+Tmp/sum(Tmp)
            }
        }
      }
      M = (Phi[,,2]*M)%*%t(Pi[,,2])
      M = M/rowSums(M)
      V[,,1] = L[,,1]*M
      if(n==1) V[,,1] = V[,,1]/sum(V[1,,1])
      else V[,,1] = V[,,1]/rowSums(V[,,1])

           # If required store parameters
      # ---- M-step ----
      # Update Mu
     if(drop) Vv = matrix(aperm(V,c(1,3,2)),n*TT,k+1) else Vv = matrix(aperm(V,c(1,3,2)),n*TT,k)
      if(miss){
        Mu00 = Mu; itc = 0  
        while((max(abs(Mu00-Mu))>10^-10 || itc==0) & itc<10){
          Mu00 = Mu; itc = itc+1
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
          
          Mub = matrix(0,r,k)
          for(u in 1:k){
            Yv1 = matrix(Y1[,,,u],n*TT)
            Mub[,u] = (t(Yv1)%*%Vv[,u])/sum(Vv[,u])
          }
          if(length(indb)==0){
            Mu = Mub
          }else{
            mu0 = Mu[indb,]%*%colSums(Vv)/sum(Vv)
            iSi = solve(Si)
            Tmp = solve(iSi[-indb,-indb])%*%iSi[-indb,indb]
            for(u in 1:k) Mu[-indb,u] = Mub[-indb,u]+Tmp%*%(Mub[indb,u]-mu0)
            Mu[indb,] = mu0
          }
          Sitmp = matrix(0,r,r)
          for(u in 1:k){
            Yv1 = matrix(Y1[,,,u],n*TT)
            Var1 = array(Var,c(n*TT,r,r))
            Tmp = Yv1-rep(1,n*TT)%*%t(Mu[,u])
            Sitmp = Sitmp+t(Tmp)%*%(Vv[,u]*Tmp)+apply(Vv[,u]*Var1,c(2,3),sum)
          }
          if(drop) Si  = Sitmp/dim(Yv2)[1] else Si = Sitmp/(n*TT) #SP: controllare se questo è necessario
          
          if(length(indb)==0) Mu00 = Mu
        }
      }else{
        Mu00 = Mu; itc = 0  
        Mub = matrix(0,r,k)
        for(u in 1:k) Mub[,u] = (t(Yv)%*%Vv[,u])/sum(Vv[,u])
        while((max(abs(Mu00-Mu))>10^-10 || itc==0) & itc<10){
          Mu00 = Mu; itc = itc+1
          if(length(indb)==0){
            Mu = Mub
          }else{
            mu0 = colMeans(Yv)[indb]
            iSi = solve(Si)
            Tmp = solve(iSi[-indb,-indb])%*%iSi[-indb,indb]
            for(u in 1:k) Mu[-indb,u] = Mub[-indb,u]+Tmp%*%(Mub[indb,u]-mu0)
            Mu[indb,] = mu0
          }
          Si = matrix(0,r,r)
          for(u in 1:k) Si= Si+ t(Yv-rep(1,n*TT)%*%t(Mu[,u]))%*%(Vv[,u]*as.matrix(Yv-rep(1,n*TT)%*%t(Mu[,u]))) #FB: velocizzato togliendo diag
          #Si = Si/(n*TT)
         if(drop) Si  = Si/dim(Yv2)[1] else Si = Si/(n*TT) #SP: controllare se questo è necessario
          if(length(indb)==0) Mu00 = Mu
        }
      }
  
      # Update piv and Pi
      piv = colSums(V[,,1])/n
      U = pmax(U,10^-300)
      if(mod==0) for(t in 2:TT) Pi[,,t] = diag(1/rowSums(U[,,t]))%*%U[,,t]
      if(mod==1){
        Ut = apply(U[,,2:TT],c(1,2),sum)
        if(drop) Pi[,,2:TT] = array(diag(1/rowSums(Ut))%*%Ut,c(k+1,k+1,TT-1))
        else Pi[,,2:TT] = array(diag(1/rowSums(Ut))%*%Ut,c(k,k,TT-1))
      }
      if(mod>1){
        Ut1 = U[,,2:mod]
        if(length(dim(Ut1))>2) Ut1 = apply(Ut1,c(1,2),sum)
        Ut2 = U[,,(mod+1):TT]
        if(length(dim(Ut2))>2) Ut2 = apply(Ut2,c(1,2),sum)
        if(drop){
          Pi[,,2:mod] = array(diag(1/rowSums(Ut1,2))%*%Ut1,c(k+1,k+1,mod-1))
          Pi[,,(mod+1):TT] = array(diag(1/rowSums(Ut2,2))%*%Ut2,c(k+1,k+1,TT-mod))
        }else{
          Pi[,,2:mod] = array(diag(1/rowSums(Ut1,2))%*%Ut1,c(k,k,mod-1))
          Pi[,,(mod+1):TT] = array(diag(1/rowSums(Ut2,2))%*%Ut2,c(k,k,TT-mod))
        }
      }
      if(drop){
        Pi[k+1,1:k,]=0
        Pi[k+1,k+1,2:TT] = 1
      }
      # Compute log-likelihood
      paro = par; par = c(piv,as.vector(Pi),as.vector(Mu),as.vector(Si))
      if(any(is.na(par))) par = par[-which(is.na(par))]
      lko = lk
      # print(c(4.5,proc.time()-t0))
      out = lk_comp_cont_MISS(Y,R,piv,Pi,Mu,Si,k,drop=drop)
      # print(c(4.7,proc.time()-t0))
      lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
      if(it/1 == round(it/1)) cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
      lkv = c(lkv,lk)
      # print(c(5,proc.time()-t0))
    }
    V2 = aperm(V,c(1,3,2))
    if(drop) V2 = aperm(array(V2,c(n,TT,k+1,r)),c(1,2,4,3)) else V2 = aperm(array(V2,c(n,TT,k,r)),c(1,2,4,3))
    if(miss) Yimp = apply(Y1*V2,c(1,2,3),sum) 
    # Compute information matrix if required
    if(out_se){
      th = NULL
      th = c(th,as.vector(Mu))
      th = c(th,Si[upper.tri(Si,TRUE)])
      th = c(th,B%*%log(piv))
      if(mod==0) for(t in 2:TT) for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,t]))
      if(mod==1) for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,2]))

      th0 = th-10^-5/2
      out = lk_obs_cont_miss(th0,Bm,Cm,k,Y,R,TT,r,mod)
      lk0 = out$lk; sc0 = out$sc
      lth = length(th)
      scn = rep(0,lth)
      J = matrix(0,lth,lth)
      for(j in 1:lth){
        thj = th0; thj[j] = thj[j]+10^-5
        out = lk_obs_cont_miss(thj,Bm,Cm,k,Y,R,TT,r,mod)
        scn[j] = (out$lk-lk0)/10^-5
        J[,j] = (out$sc-sc0)/10^-5
      }
      J = -(J+t(J))/2
      if(check_der){
        print(c(lk,lk0))
        print(round(cbind(scn,sc0,round(scn-sc0,4)),5))
      }
      #  se = sqrt(diag(ginv(J)))
      Va = ginv(J)
      nMu = r*k
      nSi = r*(r+1)/2
      Va2 = Va[1:(nMu+nSi),1:(nMu+nSi)]
      se2 = sqrt(diag(Va2))

      Va = Va[-(1:(nMu+nSi)),-(1:(nMu+nSi))]
      Om = diag(piv)-tcrossprod(piv,piv)
      M = Om%*%Bm
      if(mod==0){
        for(t in 2:TT) for(u in 1:k){
          Om = diag(Pi[u,,t])-Pi[u,,t]%o%Pi[u,,t]
          M = blkdiag(M,Om%*%Cm[,,u])
        }
      }
      if(mod==1){
        for(u in 1:k){
          Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
          M = blkdiag(M,Om%*%Cm[,,u])
        }
      }
      if(mod>1){
        for(u in 1:k){
          Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
          M = blkdiag(M,Om%*%Cm[,,u])
        }
        for(u in 1:k){
          Om = diag(Pi[u,,mod+1])-Pi[u,,mod+1]%o%Pi[u,,mod+1]
          M = blkdiag(M,Om%*%Cm[,,u])
        }
      }
      M = as.matrix(M)
      Va = M%*%Va%*%t(M)
      dVa = diag(Va)
      if(any(dVa<0)) warning("Negative elements in the estimated variance-covariance matrix for the parameters estimates")
      se = sqrt(abs(dVa))
      # Divide parameters
      se = c(se2,se)
      seMu = se[1:nMu]
      seSi = se[nMu+(1:nSi)]
      sepiv = se[nMu+nSi+(1:k)]

      if(mod==0) sePi = se[nMu+nSi+k+(1:(k*k*(TT-1)))]
      if(mod==1) sePi = se[nMu+nSi+k+(1:(k*k))]
      if(mod>1) sePi = se[nMu+nSi+k+(1:(k*k*2))]
    }
    # Compute number of parameters
    rb = length(indb)
    np = (k-1)+k*(r-rb)+rb+r*(r+1)/2
    if(drop) {
      if(mod==0) np = np+(TT-1)*k*k
      if(mod==1) np = np+k*k
      if(mod>1) np = np+2*k*k
    }
    else{
    if(mod==0) np = np+(TT-1)*k*(k-1)
    if(mod==1) np = np+k*(k-1)
    if(mod>1) np = np+2*k*(k-1)
    }
    aic = -2*lk+np*2
    bic = -2*lk+np*log(n)

    cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
    # adjust output
    #	if(any(yv!=1)) V = V/yv

    lk = as.vector(lk)
    if(drop) dimnames(Pi)=list(state=1:(k+1),state=1:(k+1),time=1:TT)
    else dimnames(Pi)=list(state=1:k,state=1:k,time=1:TT)
    
    #	dimnames(Mu) = list(dimnames(Y)[[3]],state=1:k)
    #	dimnames(Si) = list(dimnames(Y)[[3]],dimnames(Y)[[3]])
    
    if(r==1) dimnames(Mu) = list(item=1,state=1:k) else dimnames(Mu)=list(item=1:r,state=1:k)
    dimnames(Si)=list(item=1:r,item=1:r)

    rownames(Mu) <- dimnames(Y)[[3]]
    ent = -sum(V*log(pmax(V,10^-300))) # entropy
    clc = -2*lk+2*ent
    icl.bic = bic + 2*ent
    out = list(lk=lk,piv=piv,Pi=Pi,Mu=Mu,Si=Si,np=np,k = k,
               aic=aic,bic=bic, ent = ent, clc = clc, icl.bic= icl.bic,
               lkv=lkv,V=V, n = n, TT = TT, modBasic = mod)
    if(miss){
      out$Y = Y
      out$Yimp = Yimp  
    }
    if(out_se){
      seMu = matrix(seMu,r,k)
      seSi2 = matrix(0,r,r)
      seSi2[upper.tri(seSi2,TRUE)]=seSi
      seSi2 = seSi2+t(seSi2-diag(diag(seSi2)))
      seSi = seSi2
      sePi0 = sePi
      sePi = array(0,c(k,k,TT))
      if(mod>1){
        sePi0 = array(sePi0,c(k,k,2))
        sePi0 = aperm(sePi0,c(2,1,3))
        sePi[,,2:mod] = sePi0[,,1]
        sePi[,,(mod+1):TT] = sePi0[,,2]
      } else {
        sePi[,,2:TT] = sePi0
        sePi = aperm(sePi,c(2,1,3))
      }
      dimnames(sePi) = list(state=1:k,state=1:k,time=1:TT)
      if(r==1) dimnames(seMu) = list(item=1,state=1:k) else dimnames(seMu)=list(item=1:r,state=1:k)

      out$sepiv = sepiv
      out$sePi = sePi
      out$seMu = seMu
      out$seSi = seSi
    }
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    class(out)="LMbasiccont"
    return(out)
  }
