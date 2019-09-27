## Evaluating the forecast pseudo-residuals

allprobs <- matrix(1, nrow=n, ncol=N)
for(j in 1:N) 
  allprobs[which(!is.na(x)), j] <- dnorm(x, mean=mu[j], sd=stdev[j])

HMM.lalpha <- function(allprobs, gamma, delta, n,N, mu){
  lalpha<-matrix(NA,N,n)
  
  lscale <- 0
  foo <- delta*allprobs[1,]
  lscale <- 0
  lalpha[,1] <- log(foo)+lscale
  sumfoo  <- sum(foo)
  for (i in 2:n){
    foo <- foo%*%gamma*allprobs[i,]
    sumfoo <- sum(foo); lscale <- lscale+log(sumfoo); foo <- foo/sumfoo # scaling
    lalpha[,i] <- log(foo)+lscale
  }
  lalpha
}



HMM.psres <- function(x, allprobs, gamma, n,N, mu){

  
  delta <- solve(t(diag(N)-gamma + 1), rep(1, N))
  
  la<-moveHMM.lalpha(allprobs, gamma, delta, n,N, mu)
  
  pstepmat<-matrix(NA,n,N)
  fres <-rep(NA,n)
  ind.step <- which(!is.na(x))
  
  for(j in 1:length(ind.step)) {
    pstepmat[ind.step[j],1] <- pnorm(x[ind.step[j]], mean=mu[1], sd=stdev[1])
    pstepmat[ind.step[j],2] <- pnorm(x[ind.step[j]], mean=mu[2], sd=stdev[2])
  }
  
  
  if (!is.na(x[1])) fres[1]<-qnorm(rbind(c(1,0))%*%pstepmat[1,])
  for (i in 2:n){
    
    c<-max(la[,i-1])     
    a<-exp(la[,i-1]-c)
    if (!is.na(x[i])) fres[i]<-qnorm(t(a)%*%(gamma/sum(a))%*%pstepmat[i,])
  }
  return(list(fres=fres))
}

