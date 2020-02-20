data {
  int N; // number of observations
  int Nt; // number of time points per subject
  
  // counters for parameter tracking
  int n_coeffs;
  int st[Nt-1];
  int ed[Nt-1];

  
  matrix[N, Nt] A; // matrix of strata dummies
  matrix[N, Nt] L; // matrix of confounders (including constant in 1st entry)
  vector[N] Y; // outcome
  
}

parameters {
  
  // confounder model parameters
  vector[Nt] beta_int; // vector of confounder model intercepts
  vector[n_coeffs] beta_L; // vector of confounder history effects
  vector[n_coeffs] theta; // vector of treatment history effects
  
  // outcome model parameters
  real int_y;
  vector[Nt] beta_Ly; // confounder history
  vector[Nt] theta_y; // treatment history
  
  // variance parameters
  real<lower=0> phi_L0; 
  real<lower=0> phi_Lt; 
  real<lower=0> phi_y; 
}

model {
  // specify ridge prior shrinking aggressivly backward in time.
  for(t in 2:Nt){
    for(b in st[t-1]:ed[t-1] ){
     beta_L[b] ~ normal(0, ( 1.5^(b - st[t-1]) ) * 1 );
     theta[b] ~  normal(0, ( 1.5^(b - st[t-1]) ) * 1 );
    }
  }
  
  // specify standard priors on remaining parameters (moderately informative)
  int_y ~ normal(0,.5);
  beta_Ly ~ normal(0,.5);
  phi_y  ~ cauchy(0,2);
  theta_y ~ normal(0,.5);
  
  beta_int ~ normal(0, 1);
  phi_Lt ~ cauchy(0,2);
  
  phi_L0 ~ cauchy(0, 2);

  
  // model for confounder at first time point
  L[,1] ~ normal( beta_int[1], phi_L0);
  
  // sequence of confounder models, each conditioning on complete past history
  for(t in 2:Nt){
    L[ ,t] ~ normal( beta_int[t] + L[, 1:(t-1) ]*beta_L[st[t-1]:ed[t-1]] + A[, 1:(t-1) ]*theta[st[t-1]:ed[t-1]] ,  phi_Lt ) ; 
  }
  
  // outcome model conditional on treatment and confounder histories
  Y ~ normal( int_y + L*beta_Ly + A*theta_y, phi_y );
  
}


generated quantities {
  row_vector[Nt] L_pred1;
  row_vector[Nt] L_pred0;
  
  real mu1;
  real mu0;

  // draw sequentially from confounder distribution under A=1 and A=0, 
  // and compute the mean under each, mu1 and mu0 - as described in paper
  L_pred1[1] = normal_rng(  beta_int[1], phi_L0 );
  for(t in 2:Nt){
    L_pred1[t] = normal_rng( beta_int[t] + L_pred1[1:(t-1) ]*beta_L[st[t-1]:ed[t-1]] + sum(theta[st[t-1]:ed[t-1]]) ,  phi_Lt ); 
  }
  mu1 =  L_pred1*beta_Ly + sum(theta_y);

  L_pred0[1] = normal_rng(  beta_int[1], phi_L0 );
  for(t in 2:Nt){
    L_pred0[t] = normal_rng( beta_int[t] + L_pred0[1:(t-1) ]*beta_L[st[t-1]:ed[t-1]],  phi_Lt ); 
  }
  mu0 =  L_pred0*beta_Ly;

}
