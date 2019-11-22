## forward backward algorithm for steps and turning angles

log_sum_exp <- function(x1, x2){
  
  msf = max(x1 + x2)
  exps = sum(exp(x1+x2 - msf))
             
  lse = msf + log(exps)
  
  return(lse)
               
}

forwardbackward <- function(N, TT, log_gamma_tr, log_gamma, log_delta, steps, angles, shape, rate, loc, kappa){

stateProbs <- matrix(NA, nrow=TT, ncol=N)
logalpha <- matrix(NA, nrow=TT, ncol=N)
logbeta <- matrix(NA, nrow=TT, ncol=N)


## forward variables
for(n in 1:N){
  logalpha[1,n] = log_delta[n] 
  if(steps[1]>0)
      logalpha[1,n] = logalpha[1,n] + log(dgamma(x=steps[1], shape = shape[n], rate=rate[n])) 
  if(angles[1]>-pi)
    logalpha[1,n] = logalpha[1,n] + log(dvm(theta=angles[1], mu = loc[n], kappa = kappa[n]))
}

for(t in 2:TT){
  for(n in 1:N){
    logalpha[t,n] = log_sum_exp(log_gamma_tr[n,], logalpha[t-1,]) 
    if(steps[t]>0)
      logalpha[t,n] = logalpha[t,n] + log(dgamma(x=steps[t], shape = shape[n], rate=rate[n])) 
    if(angles[t]>-pi)
      logalpha[t,n] = logalpha[t,n] + log(dvm(theta=angles[t], mu = loc[n], kappa = kappa[n]))
  }
}


## backward variables
logbeta[TT,] = 0

for(t0 in (TT-1):1){
  
  temp <- logbeta[t0 + 1,]
  
  for(n in 1:N){
    if(steps[t0 + 1]>0)
      temp[n] = temp[n] + log(dgamma(x=steps[t0 + 1], shape = shape[n], rate=rate[n])) 
    if(angles[t0 + 1]>-pi)
       temp[n]= temp[n] + log(dvm(theta=angles[t0 + 1], mu = loc[n], kappa = kappa[n]))
  }
  
  for(n in 1:N){
    logbeta[t0,n] = log_sum_exp(log_gamma[n,], temp) 
  }
}

## log likelihood
llk = log_sum_exp(logalpha[TT,], rep(0, N))

## state probabilties
for(t in 1:TT){
  stateProbs[t,] = exp(logalpha[t,] + logbeta[t,] - llk)
}

return(stateProbs)

}
