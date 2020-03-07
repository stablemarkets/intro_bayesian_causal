devtools::install_github("StableMarkets/ChiRP",ref = 'fast-reg')
library(ChiRP)

library(BayesTree)
library(rstan)
library(LaplacesDemon)
library(gtools)
library(latex2exp)

## function that takes a set N x M matrix of M posterior predictions for each of 
## the N subjects under intervention 1 (mu_a1) and intervention 0 (mu_a0)
## Performs Bayesian bootstrap and returns vector containing M posterior draws
## of Psi.
bayes_boot = function(mu_a1, mu_a0){
  n = nrow(mu_a1)
  M = ncol(mu_a1)
  bb_weights = rdirichlet(1, rep(1/n, n))
  
  psi_post = numeric(M)
  for(m in 1:M){
    psi_post[m] = sum(bb_weights*( mu_a1[, m] - mu_a0[ , m] ))
  }
  
  return(psi_post)
}

#### ------------------------ Simulate Data --------------------------------####
set.seed(1)

N = 500
L = rnorm(N)
A = rbinom(N, 1, invlogit( 0 + 1*L ) )
Y = rnorm(N, 10 + (0 + L + L^2 )*A , 1)
  
plot(L, Y, col=ifelse(A==1, 'red','black'))

d_train = data.frame(Y=scale(Y), A=A, L=scale(L))

d_a1 = data.frame(A=1, L=d_train$L)
d_a0 = data.frame(A=0, L=d_train$L)

d_test = rbind(d_a1, d_a0)

#### ------------------------ Dirichlet Process ----------------------------####

set.seed(2)
res=fDPMix(d_train = d_train, formula = Y ~ L + A, d_test = d_test, 
           iter=1000, burnin=500,init_k = 10)

psi_dp = bayes_boot( mu_a1 = res$test[1:N,], mu_a0 = res$test[(N+1):(2*N), ])

#### ----------------------- BART Model         ----------------------------####

set.seed(3)
bart_res = bart(x.train = d_train[, c('L','A') ], y.train = d_train$Y, 
                x.test = d_test, ndpost = 500)

bart_res_pred = t(bart_res$yhat.test)

mu_a1 = bart_res_pred[1:N, ]
mu_a0 = bart_res_pred[(N+1):(2*N), ]

psi_bart = bayes_boot( mu_a1 = mu_a1, mu_a0 = mu_a0)

plot(d_train$L, d_train$Y)
points(d_train$L, rowMeans(mu_a1), col='red')
points(d_train$L, rowMeans(mu_a0), col='blue')

#### ----------------------- Gaussian Process   ----------------------------####

x_train = model.matrix(~ -1 + A + L, data = d_train )
x_test = model.matrix(~ -1 + A + L, data = d_test )

stan_data <- list(N2=nrow(d_test), 
                  D = ncol(x_test),
                  N1 = nrow(d_train),
                  x2 = x_test, 
                  x1 = x_train, y1 = d_train$Y)

## https://mc-stan.org/docs/2_22/stan-users-guide/gaussian-processes-chapter.html
## GP code from Stan site that implements Guassian Process regression with 
## inference for hyperparamters.
## see above link for guidance on hyperprior calibration and other great info.
stan_mod = stan_model("gaussian_process_with_HPs_multi.stan")

## In practice, should run for many more iterations and 
## check convergence diagnostics
stan_res = sampling(stan_mod, data=stan_data, 
                    chains=1, warmup=500, iter=1000, seed=1000)

gp_res_test=extract(stan_res, pars=c('f_test'))[[1]]

gp_res_pred = t(gp_res_test)

mu_a1 = gp_res_pred[1:N, ]
mu_a0 = gp_res_pred[(N+1):(2*N), ]

psi_gp = bayes_boot( mu_a1 = mu_a1, mu_a0 = mu_a0)
hist(psi_gp)


#### ----------------------- Linear Regression    --------------------------####

## compute bootstrap distribution for psi_lm
set.seed(123)
B=500
psi_lm_boot = numeric(B)
for(b in 1:B){
  id = sample(1:N, size = N, replace = T)
  
  d_train_b = d_train[id,]
  
  d_a1_b = data.frame(A=1, L=d_train_b$L)
  d_a0_b = data.frame(A=0, L=d_train_b$L)
  
  d_test_b = rbind(d_a1_b, d_a0_b)

  lm_b = lm(data = d_train_b, Y ~ A + L)
  pred_lm_b = predict(lm_b, d_test_b)
  psi_lm_boot[b] = mean(pred_lm_b[1:N] - pred_lm_b[(N+1):(2*N)])
  
}


#### ----------------------- Plot Results       ----------------------------####

png("bnp_ate.png")
all_psi = cbind(psi_bart, psi_dp, psi_gp, psi_lm_boot)

post_mean = colMeans(all_psi)


plot(post_mean, pch=20, col='blue', ylim= c(0, 2),
     ylab=TeX("$\\Psi$"), xlab=TeX("Model"), axes=F )

axis_labs = c('BART', 'Dirichlet Process', 'Gaussian Process', 'Linear Model')

axis(1, at = 1:4, labels = axis_labs, padj = .5 )
axis(2, at = seq(0,2, .5), labels = seq(0,2, .5) )


### Plot posterior credible Intervals 
colfunc <- colorRampPalette(c("white", "lightblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

for(i in ci_perc){
  pci = apply(all_psi, 2, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  segments(1:4,pci[1,], 1:4, pci[2,], col=colvec[as.character(i)], lwd=10 )
}
###

points(post_mean, col='steelblue', pch=20, lwd=8)
abline(h=1, col='red', lty=2)

legend('topleft',
       legend = c('Posterior Mean/Credible Band', 'True Value'),
       col = c('lightblue', 'red'), pch=c(20,NA), lty=c(1, 2), lwd=c(8,1), bty='n')

legend('topleft',
       legend = c('', ''),
       col = c('steelblue', 'red'), pch=c(20,NA), lty=c(NA, 2), lwd=c(8,1), bty='n')
dev.off()




