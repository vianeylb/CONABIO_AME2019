// forward-backward algorithm in Stan


generated quantities {

  real stateProbs[T,N];
  vector[N] lp;
  vector[N] lp_p1;
  
  // forward-backward algorithm (state probabilities)
  {
    real logalpha[T,N];
    real logbeta[T,N];
    real llk;
    
    // log alpha probabilities
    for(t in 1:T) {
      if(t==1 || ID[t]!=ID[t-1]) {
        for(n in 1:N)
          lp[n] = -log(N);
      }
      
      for (n in 1:N) {
        lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t,n]) + lp);
        if(steps[t]>=0)
          lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t] | shape[n], rate[n]);
        if(angles[t]>=(-pi())) {
          lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
        }
        logalpha[t,n] = lp_p1[n];
      }
      lp = lp_p1;
    }
    
    // log beta probabilities
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      
      if(t==T || ID[t+1]!=ID[t]) {
        for(n in 1:N)
          lp_p1[n] = 0;
      } else {
        for(n in 1:N) {
          if(steps[t+1]>=0)
            lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t+1] | shape[n], rate[n]);
          if(angles[t+1]>=(-pi()))
            lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t+1] | loc[n], kappa[n]);
        }
        
        for(n in 1:N){
         lp_p1[n] = log_sum_exp(to_vector(log_theta[t+1,n]) + lp_p1);
        }
      }
      lp = lp_p1;
      for(n in 1:N)
        logbeta[t,n] = lp[n];
    }

    // state probabilities
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      if(t==T || ID[t+1]!=ID[t])
        llk = log_sum_exp(logalpha[t]);
      for(n in 1:N)
        stateProbs[t,n] = exp(logalpha[t,n] + logbeta[t,n] - llk);
        
      //stateProbs[t] = stateProbs[t]/sum(stateProbs[t]);
    }
  }
}

