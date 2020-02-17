data {
  int N; // number of observations
  int P; // number of parameters
  int num_A_levels; // number of dose dummies
  matrix[N, P] L; // matrix of confounders (including constant)
  matrix[N, num_A_levels] A; //matrix of dose dummies
  vector[N] Y; // outcome
}

parameters {
  vector[P] beta; // coefficients of confounder vector
  vector[num_A_levels] theta; // coefficients of dose 
  real<lower=0> phi; // conditional outcome distribution variance
}

model {
  
  // set priors as described in paper tau_k = 1 for all k
  theta[1] ~ normal( 0, 10 );
  theta[2] ~ normal( 2*theta[1], 1);
  
  for(j in 3:num_A_levels){
    theta[j] ~ normal( 2*theta[j-1] - theta[j-2], 1 );
  }
  
  beta ~ normal(0, 10); 
  phi ~ cauchy(0,10);
  
  // specify likelihood 
  Y ~ normal(L*beta + A*theta , phi);
}


// transform thetas to form Psi's (the dose-response curve)
generated quantities {
  vector[num_A_levels] Psi;
  
  Psi[1] = theta[1];
  for(k in 2:num_A_levels){
    Psi[k] = theta[k] - theta[k-1];
  }
}
