functions {
  
  // generate draw of function from GP posterior
  vector gp_pred_rng(matrix x_test, vector y_train, matrix x_train, real alpha,
                     real rho, real sigma, real delta) {
    int N_train = rows(x_train);
    int N_test = rows(x_test);
    
    vector[N_test] f_test;
    
    {
      matrix[N_train, N_train] L_K;
      vector[N_train] K_div_y1;
      matrix[N_train, N_test] k_x1_x2;
      matrix[N_train, N_test] v_pred;
      vector[N_test] f2_mu;
      matrix[N_test, N_test] cov_f2;
      matrix[N_test, N_test] diag_delta;
      matrix[N_train, N_train] K;
      
      for (n_i in 1:N_train){
        for(n_j in 1:N_train){
          K[n_i, n_j] = cov_exp_quad(x_train[n_i], x_train[n_j], alpha, rho);  
        }
      }
      
      for (n in 1:N_train)
        K[n, n] = K[n,n] + square(sigma);
        
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';

      for (n_i in 1:N_train){
        for(n_j in 1:N_test){
          k_x1_x2[n_i, n_j] = cov_exp_quad(x_train[n_i], x_test[n_j], alpha, rho);  
        }
      }
      
      f_test_mu = (k_x1_x2' * K_div_y1);
      
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      
      cov_f_test = cov_exp_quad(x_test, alpha, rho) - v_pred' * v_pred;
      
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f_test = multi_normal_rng(f_test_mu, cov_f_test + diag_delta);
    }
    return f_test;
  }
}

data {
  int<lower=1> N_train; // training data sample size
  int<lower=1> P; // number of columns in model matrix
  
  // training data, could set these to be row vectors for several variables
  matrix[N_train, P] x_train[N1]; 
  vector[N_train] y_train; // training outcome

  int<lower=1> N_test; // test data sample size
  matrix[N_test, P] x_test; // test data 
}

transformed data {
  vector[N_train] mu = rep_vector(0, N_train);
  real delta = 1e-9;
}

parameters {
  // inference for hyperparameters
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  matrix[N_train, N_train] L_K;

  {
    matrix[N_train, N_train] K;
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n_i in 1:N_train){
      for(n_j in 1:N_train){
        K[n_i, n_j] = cov_exp_quad(x_train[n_i], x_train[n_j], alpha, rho);  
      }
    }
    
    for(n in 1:N_train) K[n, n] = K[n,n] + sq_sigma;
    
    L_K = cholesky_decompose(K);
  }

  rho ~ inv_gamma(1, 1);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 3 );

  y_train ~ multi_normal_cholesky(mu, L_K);
}
generated quantities {
  vector[N_train] f_train;
  vector[N_test] f_test;

  f_train = gp_pred_rng(x_train, y1, x_train, alpha, rho, sigma, delta);
  f_test = gp_pred_rng(x_test, y1, x_train, alpha, rho, sigma, delta);

}
