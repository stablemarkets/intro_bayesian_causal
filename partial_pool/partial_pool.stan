data {
  int N; // number of observations
  int Pv;
  int Pw;
  matrix[N, Pv] V; // matrix of strata dummies
  matrix[N, Pw] W; // matrix of confounders (including constant in 1st entry)
  vector[N] A; // treatment indicator
  int Y[N]; // outcome
}

parameters {
  vector[Pw] beta_w; // confounder effects and intercept
  vector[Pv] beta_v; // strata main effects
  vector[Pv] theta_1q; // interaction effects
  real theta; // treatment main effect
  real mu;
}

model {
  
  // specify priors
  beta_w ~ normal(0, 1); 
  beta_v ~ normal(0, 1); 
  theta ~ normal( mu, 3 );
  
  theta_1q ~ normal(mu - theta, .5);
  
  // specify likelihood 
  for(i in 1:N){
    Y[i] ~ bernoulli_logit( W[i]*beta_w + V[i]*beta_v + (theta + V[i]*theta_1q)*A[i] );
  }
}


// transform thetas to form Psi's (the dose-response curve)
generated quantities {
  // take draw of bayesian bootstrap weights
  vector[N] bb_weights = dirichlet_rng( rep_vector( inv(N), N) ) ;
  
  vector[N] cond_mean_y1;
  vector[N] cond_mean_y0;
  real marg_mean_y1;
  real marg_mean_y0;
  vector[Pv+1] odds_ratio;
  
  odds_ratio[1] = exp(theta);
  
  // cycle through strata of interest and compute causal Odds Ratio for each.
  for( v in 1:(Pv) ){
    
    // compute mean under strata v, under both A=1 and A=0
    for(i in 1:N){
      cond_mean_y1[i] = inv_logit( W[i]*beta_w + beta_v[v] + theta + theta_1q[v] );
      cond_mean_y0[i] = inv_logit( W[i]*beta_w + beta_v[v] );
    }
    
    // taking average over bayesian bootstrap weights
    marg_mean_y1 = bb_weights' * cond_mean_y1 ;
    marg_mean_y0 = bb_weights' * cond_mean_y0 ;
    
    // compute odds ratio 
    odds_ratio[v+1] = (marg_mean_y1/(1 - marg_mean_y1)) / (marg_mean_y0/(1 - marg_mean_y0));
  }
  
}
