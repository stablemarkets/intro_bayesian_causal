data {
  int N; // number of observations
  vector[N] A; // treatment
  vector[N] L; // confounder
  vector[N] Y; // outcome
}

parameters {
  real theta; // treatment differential estimate.
  real beta; // confounder coefficient
  real<lower=0> phi; // conditional outcome distribution variance
  real<lower=0> delta; // sensitivity parameter
  real delta1; // sensitivity parameter
}

model {
  
  // set priors 
  theta ~ normal(0, 3);
  beta ~ normal(0, 3);
  delta ~ gamma(1, 1);
  delta1 ~ normal(0,1);
  
  // specify likelihood 
  Y ~ normal(A*theta + beta*L, phi);
}


// transform thetas to form Psi's (the dose-response curve)
generated quantities {
  real psi1;
  real psi2;
  real psi3;
  
  psi1 = theta + 1*delta;
  psi2 = theta - 1*delta;
  psi3 = theta + delta1;
}
