bootstrap.MISS <- function(Y,X1,X2,Mu,Si,Be,Ga,B=100){

# preliminaries
	mMu = mSi = mBe = 0
	m2Mu = m2Si = m2Be = 0
	mGa = 0
	m2Ga = 0

  if(is.vector(Mu)){
    	r =1
    	k = length(Mu)
  }else{
    	r = nrow(Mu)
    	k = ncol(Mu)
  }

  for (b in 1:B) {
      cat("non-parametric boostrap sample n. ",b,"\n")
      n <- dim(Y)[1]
      ind = sample(n,n,replace=T)
      Yb = Y[ind,,]
      X1b = X1[ind,]
      X2b = X2[ind,,]
    
    out <- lmcovlatent.cont.MISS(Y=Yb,X1=X1b,X2=X2b,k=k,output=TRUE)
    
    mMu = mMu + out$Mu/B
    mSi = mSi + out$Si/B
    mBe = mBe + out$Be/B
   	mGa = mGa + out$Ga/B
    m2Ga = m2Ga + out$Ga^2/B   	
    
    m2Mu = m2Mu + out$Mu^2/B
    m2Si = m2Si + out$Si^2/B
    m2Be = m2Be + out$Be^2/B
 
  }
  seMu = sqrt(m2Mu - mMu^2)
  seSi = sqrt(m2Si - mSi^2)
  seBe = sqrt(m2Be - mBe^2)
  seGa = sqrt(m2Ga - mGa^2)
  out = list(mMu = mMu, mSi = mSi, mBe = mBe, mGa = mGa,
                 seMu = seMu, seSi = seSi, seBe = seBe, seGa = seGa)
  
}