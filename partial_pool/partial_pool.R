###--------------------------------------------------------------------------###
### Author: Arman Oganisian
### Conduct partial pooling analysis
###--------------------------------------------------------------------------###

## Load Packages
library(rstan)
library(LaplacesDemon)
library(latex2exp)
set.seed(66)

####------------------------ Simulate Data ---------------------------------####
N = 500 # sample size
warmup = 1000
iter = 2000
n_draws = iter - warmup

#simulate standard normal confounder (L), dose (A), and outcome (Y)
W = rnorm(n = N) 
Wmat = cbind(1, W)

V = sample(x = 1:5, size = N, replace = T, prob = c(3/10,3/10,2/10,1/10,1/10))
Vmat = model.matrix(~ -1 + as.factor(V) )[,-1]

A = rbern(N, prob = invlogit( 0 + 1*W + Vmat %*% c( -.5,.5,.5,-.5 ) ) )

Y = rbern(N, prob = invlogit( -1 + 1*W  + (1 + Vmat %*% c( -.5,0,.5, .6 ) )*A ) ) 

stan_data = list(Y=Y, A=A, V=Vmat, W = Wmat,
                 Pw = ncol(Wmat), Pv = ncol(Vmat), N=N)

####------------------------ Sample Posterior    ---------------------------####
partial_pool_model = stan_model(file = "partial_pool.stan")

stan_res = sampling(partial_pool_model, data = stan_data, 
                    pars = c("odds_ratio", 'mu', "overall_odds_ratio" ),
                    warmup = warmup, iter = iter, chains=1, seed=1)

Psi_draws = extract(stan_res, pars='odds_ratio')[[1]]
overall_mean = extract(stan_res, pars='overall_odds_ratio')[[1]]

####------------------- Compute Frequentist Estimates   --------------------####

vfac = factor(V)

## estimate outcome model, which we'll then integrate over p(W)
freq_reg = glm(Y ~  W + vfac + A + vfac*A, family=binomial('logit') )

## loop through strata of interest
Psi_freq = numeric(5) ## shell to contain causal odds ratios
for(vf in 1:5){
  vval= as.factor(vf)
  
  ## standardize over empirical distribution p(W)
  marg_mean_y1 = mean( predict(freq_reg, data.frame(W, vfac=vval, A=1 ), type='response') )
  marg_mean_y0 = mean( predict(freq_reg, data.frame(W, vfac=vval, A=0 ), type='response') )
  
  ## compute causal odds ratio
  Psi_freq[vf] = (marg_mean_y1/(1-marg_mean_y1)) / (marg_mean_y0/(1-marg_mean_y0))
}


####-------------------         Plot Results            --------------------####
v_strata = 1:(ncol(Vmat)+1)
post_mean = colMeans(Psi_draws)

png("ppooling_plot.png")

plot(post_mean, pch=20, col='blue', ylim=c(0,5),
     ylab=TeX("$\\Psi(v)$"), xlab=TeX("$V$"), axes=F )

axis_labs = paste0(v_strata, "\n(n=",table(V),")")

axis(1, at = 1:5, labels = axis_labs, padj = .5 )
axis(2, at = 0:5, labels = 0:5 )


### Plot posterior credible Intervals 
colfunc <- colorRampPalette(c("white", "lightblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

for(i in ci_perc){
  pci = apply(Psi_draws, 2, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  segments(1:5,pci[1,], 1:5, pci[2,], col=colvec[as.character(i)], lwd=10 )
}
###

points(post_mean, col='steelblue', pch=20, lwd=8)
points(Psi_freq, col='black', pch=20, lwd=5)
abline(h= mean(overall_mean), col='steelblue', lty=2)

legend('topleft', 
       legend = c('Posterior Mean/Credible Band', 'MLE', "Overall Effect" ),
       col = c('steelblue', 'black', 'steelblue' ), pch=c(20,20,NA),
       lty = c(NA,NA,2), bty='n')

dev.off()