draw.lm.basic.cont.MISS <-
function(piv,Pi,Mu,Si,n,pmiss){

# Draw a sample of size n from a Basic Latent Markov model 
# for continuous data with dropout (coded as 999) with parameter piv, Pi, Mu and Si


# Preliminaries
	if(is.vector(Mu)){
    	r = 1
    	k = length(Mu)
    	Mu = matrix(Mu,r,k)
    }else{
    	r = nrow(Mu)
    	k = ncol(Mu)
    }
	TT = dim(Pi)[3]
	if(r==1) Si = matrix(Si,r,r)
# For each subject
    Y = array(0,c(n,TT,r))
    cat("------------|\n")
    cat(" sample unit|\n")
    cat("------------|\n")
    U = matrix(0,n,TT)
    for(i in 1:n){
    		if(i/1000==floor(i/1000)) cat(sprintf("%11g",i),"\n",sep=" | ")
        piv<-piv[1:k]
    		U[i,1] = k+1-sum(runif(1)<cumsum(piv))
		    Y[i,1,] = rmvnorm(1,Mu[,U[i,1]],Si)
   		  for(t in 2:TT){
   		    if(U[i,t-1]==(k+1)){
   		      Y[i,t,] = 999
   		      U[i,t] = (k+1)
   		    }else{
   		      U[i,t] = which(rmultinom(1,1,Pi[U[i,t-1],,t])==1)
   		      if(U[i,t]==(k+1)){
   		        Y[i,t,] = 999
   		      } else{
   		        Y[i,t,] = rmvnorm(1,Mu[,U[i,t]],Si)
   		      }
    			} 
   		}
    }
    #generate intermittent missingness
    nmiss <- pmiss*(n*TT)
    Ti <- rep(0,n)
    for(i in 1:n) Ti[i] <- sum(U[i,]<(k+1))
    
    #generate intermittent missing data 
    indmissn <- sample(1:n,size=nmiss,replace=TRUE)
    indmissT <- rep(0,nmiss)
    c <-0
    for(j in indmissn){
      c <- c+1
      indmissT[c] <- sample(1:Ti[j],size=1)
    }
    
    for(i in 1:nmiss){
      if(indmissT[i]==1){
        indmissr <- sample(1:r,size=1)
        Y[indmissn[i],indmissT[i],indmissr] = NA
      }else Y[indmissn[i],indmissT[i],] = NA
    } 
    
    cat("------------|\n")
    	out = list(Y=Y,U=U)
}
