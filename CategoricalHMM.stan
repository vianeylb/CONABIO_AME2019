// Categorical HMM

data {

  int<lower=0> N; // number of states
  int<lower=1> T; // length of data set
  int<lower=1> y[T]; // observations
  int<lower=1> cats;

  vector[cats] as[N];

}


parameters {

  simplex[cats] theta[N];

  simplex[N] init;
  simplex[N] tpm[N]; // N x N tpm
  ordered[N] mu; // state-dependent parameters

}  

model {

  vector[N] log_tpm_tr[N];
  vector[N] lp;
  vector[N] lp_p1;
  
  // prior for theta
  
  theta[1] ~ dirichlet(as[1]);
  theta[2] ~ dirichlet(as[2]);
  

  // transpose the tpm and take natural log of entries
  for (n1 in 1:N)
    for (n2 in 1:N)
    log_tpm_tr[n2, n1] = log(tpm[n1, n2]);

  
  // forward algorithm implementation
  for(n in 1:N) // first observation
    lp[n] = log(init[n]) + categorical_lpmf(y[1] | theta[n]);

  for (t in 2:T) { // looping over observations
    for (n in 1:N) // looping over states
      lp_p1[n] = log_sum_exp(log_tpm_tr[n] + lp) + 
        categorical_lpmf(y[t] | theta[n]); 
      
    lp = lp_p1;
  }

  target += log_sum_exp(lp);
}

