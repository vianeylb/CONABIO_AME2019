// Markov-Switching Regression Model

data {
  
  //observations
  int<lower=1> TT;
  vector[TT] y;
  
  //states
  int<lower=1> N;

  //covariates
  int<lower=0> ncovs;
  vector[ncovs + 1] xmat[TT];
  
}

parameters {
  
  //observation process
  //real mu[N];
  vector[ncovs+1] beta[N];
  real<lower=0> sigma[N];
  
  //state process
  simplex[N] tpm[N];
  simplex[N] init;
}

model {
  
  vector[N] lp;
  vector[N] lp_temp;
  vector[N] log_tpm_tr[N];
  
  
  sigma ~ student_t(3, 0, 1);
  for(i in 1:2){
    beta[i] ~ normal(0, 3);
  }
  
  for(n1 in 1:N){
    for(n2 in 1:N){
      log_tpm_tr[n1, n2] = log(tpm[n2, n1]);
    }
  }
  
  //forward algorithm  
  for(n in 1:N){
    lp[n] = log(init[n]) + normal_lpdf(y[1]| beta[n].*xmat[1], sigma[n]);
  }
  
  for(t in 2:TT){
    for(n in 1:N){
      lp_temp[n] = log_sum_exp(log_tpm_tr[n] + lp) + 
                   normal_lpdf(y[t] | sum(beta[n].*xmat[t]), sigma[n]);
    }
    
    lp = lp_temp;
  }

  target += log_sum_exp(lp);

}

