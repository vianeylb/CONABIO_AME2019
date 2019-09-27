// Markov-Switching Splines Model

data {
  
  //observations
  int<lower=1> TT;
  vector[TT] y;
  
  //states
  int<lower=1> N;

  //covariates
  vector[TT] xmat;
  int<lower=1> num_basis; 
  matrix[num_basis, TT] B; 
  
}

parameters {
  
  //observation process
  row_vector[num_basis] a_raw[N]; 
  real a0[N]; 
  real<lower=0> tau[N]; 
  real<lower=0> sigma[N];
  
  //state process
  simplex[N] tpm[N];
  simplex[N] init;
}

transformed parameters {
  row_vector[num_basis] a[N];
  vector[TT] y_hat[N];
  
  for(n in 1:N){
    a[n,1] = a_raw[n,1];
  }
  
  for (i in 2:num_basis){
    for(n in 1:N)
      a[n,i] = a[n,i-1] + a_raw[n,i]*tau[n];
  }
  
  for(n in 1:N){
    y_hat[n] = a0[n]*to_vector(xmat) + to_vector(a[n]*B);
  }
}

model {
  
  vector[N] lp;
  vector[N] lp_temp;
  vector[N] log_tpm_tr[N];
  
  
  // Priors
  for(n in 1:N)
   a_raw[n] ~ normal(0, 1);
   
  a0 ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  for(n1 in 1:N){
    for(n2 in 1:N){
      log_tpm_tr[n1, n2] = log(tpm[n2, n1]);
    }
  }
  
  //forward algorithm  
  for(n in 1:N){
    lp[n] = log(init[n]) + normal_lpdf(y[1]| y_hat[n,1], sigma[n]);
  }
  
  for(t in 2:TT){
    for(n in 1:N){
      lp_temp[n] = log_sum_exp(log_tpm_tr[n] + lp) + 
                   normal_lpdf(y[t] | y_hat[n,t], sigma[n]);
    }
    
    lp = lp_temp;
  }

  target += log_sum_exp(lp);

}
