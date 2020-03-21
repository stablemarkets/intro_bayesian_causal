functions {
  
  // generate draw of function from GP posterior
  vector gp_pred_rng(vector[] x2, vector y1, vector[] x1, real alpha,
                     real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    
    {
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      
      K = cov_exp_quad(x1, alpha, rho);
      
      for (n in 1:N1)
        K[n, n] = K[n,n] + square(sigma);
        
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}

data {
  int<lower=1> N1; // training data sample size
  int<lower=1> D;
  // training data, could set these to be row vectors for several variables
  vector[D] x1[N1];
  vector[N1] y1; // training outcome

  int<lower=1> N2; // test data sample size
  vector[D] x2[N2]; // test data 
}

transformed data {
  vector[N1] mu = rep_vector(0, N1);
  real delta = 1e-9;
}

parameters {
  // inference for hyperparameters
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  matrix[N1, N1] L_K;

  {
    matrix[N1, N1] K = cov_exp_quad(x1, alpha, rho);
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n1 in 1:N1)
      K[n1, n1] = K[n1, n1] + sq_sigma;

    L_K = cholesky_decompose(K);
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();

  y1 ~ multi_normal_cholesky(mu, L_K);
}
generated quantities {
  vector[N2] f_test;
  vector[N1] f_train;

  f_test = gp_pred_rng(x2, y1, x1, alpha, rho, sigma, delta);
  f_train = gp_pred_rng(x1, y1, x1, alpha, rho, sigma, delta);

}
