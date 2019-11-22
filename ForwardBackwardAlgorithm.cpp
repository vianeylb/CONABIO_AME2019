/* Forward/backward algorithm using the log_sum_exp function to ensure
 * numerical stability
 * 
 * N = # of states in HMM
 * 
 * n = length of time series
 * 
 * log_gamma = transition probability matrix for a homogeneous state process
 *             log transformed 
 *             :: NxN matrix
 * 
 * log_tr_gamma = transition probability matrix for a homogeneous state process
 *                transposed and log transformed 
 *                :: NxN matrix
 * 
 * log_allprobs = logarithm of the state-dependent densities evaluated 
 *                at each observation 
 *                :: n x N matrix
 *                
 * log_delta  = log initial distribution 
 *              :: vector of length N
 *              
 */


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_sum_exp(NumericVector foo1, NumericVector foo2){

  double msf;
  NumericVector sumfoos = foo1 + foo2;
  double exps;
  double lse;
  int n = foo1.size();

  msf = max(sumfoos);
  exps = 0.0; 
  for(int i=0; i< n; i++){
    exps+=exp(sumfoos[i] -msf);
  }

  lse = msf + log(exps);
  
  return(lse);
}

// [[Rcpp::export]]
NumericMatrix forbackalg(int n, int N, NumericVector log_foo, NumericMatrix log_tr_gamma, NumericMatrix log_gamma, NumericMatrix log_allprobs) {

  double lscale;
  NumericVector zero(N); /* defaults to a vector of zeros */

  NumericMatrix logalpha(n, N); /* a numeric matrix of zeros */
  NumericMatrix logbeta(n, N); /* a numeric matrix of zeros */
  NumericMatrix stateprobs(n, N); /* a numeric matrix of zeros */
  
  
  /* forward algorithm */
  
  for(int j=0; j<N; j++){ 
    logalpha(0,j)= log_foo[j] + log_allprobs(0,j);
  }

  for (int i=1; i < n; i++){
    for(int j=0; j<N; j++){ 
      logalpha(i,j) = log_sum_exp(logalpha.row(i-1), log_tr_gamma.row(j)) + log_allprobs(i,j);
    }
  }
  
  
  /* likelihood */
  
  lscale = log_sum_exp(logalpha.row(n-1), zero);
  
  /*backward algorithm */
  
   for(int j=0; j<N; j++){ 
     logbeta(n-1,j) = 0;
  }
  
  for (int i = n-2 ; i > -1; i--){
    for(int j=0; j<N; j++){ 
      logbeta(i,j) = log_sum_exp(logbeta.row(i+1), log_gamma.row(j)) + log_allprobs(i+1,j);
    }
    
  }
  
  /* evaluating state probabilities */
  
  for(int i=0; i < n; i++){
    for(int j=0; j< N; j++){
      stateprobs(i,j) = exp(logalpha(i,j) + logbeta(i,j) - lscale);
    }
  }
  
  
  return(stateprobs);
  
}


