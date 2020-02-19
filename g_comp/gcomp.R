###--------------------------------------------------------------------------###
### Author: Arman Oganisian
### Run g-computation using Ridge prior
###--------------------------------------------------------------------------###

## Load Packages
library(rstan)
library(LaplacesDemon)
library(latex2exp)
set.seed(1)


####------------------------ Simulate Data ---------------------------------####
Nt = 10 # number of time points
N = 100 # number of subjects 

L = matrix(NA, ncol=Nt, nrow=N)
A = matrix(NA, ncol=Nt, nrow=N)

L[,1] = rnorm(N)
A[,1] = rbern(N, prob = invlogit(0 + 1*L[,1] ) )

for(t in 2:Nt){
  hist_idx = 1:(t-1)
  
  beta_L = 2*(1 - 1/hist_idx)
  theta =  -.5*(1 - 1/hist_idx)
  
  L[,t] = rnorm(N, mean = L[, hist_idx, drop=F] %*% beta_L + A[,hist_idx, drop=F] %*% theta, sd = 1)  
  A[,t] = rbern(N, prob = invlogit(0 + A[,hist_idx, drop=F] %*% theta + L[,1:t, drop=F] %*% ( 1-1/(1:t) )  ) )
  
}

Y = rnorm(N, A %*% rep(1, Nt) + L %*% rep(-1, Nt), 1 )

## create position trackers for stan coefficients (this is purely computational)
## as Stan doesn't allow ragged arrays.
beta_len  = function(x) sum( x - 1:(x-1) ) ## get lenght of parameter vector for model at time X
endidx = sapply(2:Nt, beta_len   )
stidx = sapply(2:Nt, function(x) beta_len(x-1) + 1 )
stidx[1] = 1

## create list of data to feed into Stan
stan_data = list(Y=Y, A=A, L=L, st=stidx, ed=endidx, 
                 N=N, Nt=Nt, n_coeffs = max(endidx) )

####------------------------ Run Stan Model --------------------------------####


gcomp_model = stan_model(file = "gcomp.stan")

stan_res = sampling(gcomp_model, data = stan_data, 
                    control = list(max_treedepth = 10),
                    warmup = 500, iter = 1000, chains=1, seed=1)

beta_draws = extract(stan_res, pars='beta_L')[[1]]
mu1 = extract(stan_res, pars='mu1')[[1]]
mu0 = extract(stan_res, pars='mu0')[[1]]

## confounder history coefficients for L_10 model are entries 37-45
stidx[9]
endidx[9]
## order: beta_1, beta_2, ..., beta_9.
## reverse to go from most
beta_draws_last = beta_draws[, 37:45]
post_means = rev(colMeans(beta_draws_last))


png("gcomp_sparse.png")
plot(0:8, post_means, type='l', ylim=c(-.2,3), 
     ylab='Coefficient Value', axes=F, 
     xlab=TeX("Coefficients $\\beta_t$: effect of $L_t$ on $L_9$, for $t\\in 0,...8$ "))
axis(side = 1, 
     at = 0:8, 
     labels = c(TeX("$\\beta_8$"),TeX("$\\beta_7$"), TeX("$\\beta_6$"),
                TeX("$\\beta_5$"), TeX("$\\beta_4$"), TeX("$\\beta_3$"),
                TeX("$\\beta_2$"), TeX("$\\beta_1$"), TeX("$\\beta_0$"))  )
axis(side = 2, at = seq(-1,3,.5), labels = seq(-1,3,.5))

### Plot posterior credible Band 
colfunc <- colorRampPalette(c("white", "skyblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

for(i in ci_perc){
  pci = apply(beta_draws_last, 2, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  polygon(c(0:8,rev(0:8)),c( rev(pci[1,]),rev( rev(pci[2,]) )),
          col = colvec[as.character(i)], border = FALSE)
}
###
lines(0:8, post_means, type='l', col='steelblue')
abline(h=0, lty=2)

legend("topright", 
       legend = c('Posterior Mean/Credible Band', 'No Effect'),
       col = c('steelblue', 'black'), lty=c(1,2), bty='n' )
dev.off()

png('gcomp_psi.png')
hist( (mu1 - mu0 )/ sd(mu1 - mu0), freq=F,
      xlab=TeX('Posterior Draws of $\\Psi$'), col='skyblue',
      main='', ylim=c(0,.401))
lines(density((mu1 - mu0 )/ sd(mu1 - mu0)), col='steelblue')
dev.off()

