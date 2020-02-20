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



beta_L_all = rev(c(-1,1,-.5,.5,-.25,.25, 0, 0, 1,0))

for(t in 2:Nt){
  hist_idx = 1:(t-1)
  
  beta_L = beta_L_all[hist_idx]
  theta =  -.5*(1 - 1/hist_idx)
  
  L[,t] = rnorm(N, mean = L[, hist_idx, drop=F] %*% beta_L + A[,hist_idx, drop=F] %*% theta, sd = 10)  
  A[,t] = rbern(N, prob = invlogit(0 + A[,hist_idx, drop=F] %*% theta + L[,1:t, drop=F] %*% ( 1-1/(1:t) )  ) )
  
}

Y = rnorm(N, A %*% rep(1, Nt) - L %*% rep(-1, Nt), 1 )

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
                    warmup = 1000, iter = 3000, chains=1, seed=1)

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
plot(0:8, post_means, type='l', ylim=c(-2,2), 
     ylab='Coefficient Value', axes=F, 
     xlab=TeX("Coefficients $\\beta_t$: effect of $L_t$ on $L_9$, for $t\\in 0,...8$ "))
axis(side = 1, 
     at = 0:8, 
     labels = c(TeX("$\\beta_8$"),TeX("$\\beta_7$"), TeX("$\\beta_6$"),
                TeX("$\\beta_5$"), TeX("$\\beta_4$"), TeX("$\\beta_3$"),
                TeX("$\\beta_2$"), TeX("$\\beta_1$"), TeX("$\\beta_0$"))  )
axis(side = 2, at = seq(-2,2,.5), labels = seq(-2,2,.5))

### Plot posterior credible Band 
colfunc <- colorRampPalette(c("white", "skyblue"))
ci_perc = seq(.95,.05,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

for(i in ci_perc){
  pci = apply(beta_draws_last, 2, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  polygon(c(0:8,rev(0:8)),c( rev(pci[1,]),rev( rev(pci[2,]) )),
          col = colvec[as.character(i)], border = FALSE)
}
###
lines(0:8, post_means, type='o', col='steelblue', pch=20)
abline(h=0, lty=2)

tt=lm(L[,Nt]~ L[,1:(Nt-1)] + A[,1:(Nt-1)] )
points(0:8, rev(tt$coefficients[2:10]), pch=20 )
points(0:8, rev(beta_L_all[1:9]), col='red', pch=20)

legend("topleft", 
       legend = c('Posterior Mean/Credible Band', 'MLE', 'Truth'),
       col = c('steelblue', 'black', 'red'), lty=c(1,NA, NA), pch=c(20,20, 20), bty='n' )
dev.off()

png('gcomp_psi.png')
hist( (mu1 - mu0 ), freq=F,
      xlab=TeX('Posterior Draws of $\\Psi$'), col='skyblue',
      main='', breaks=30)
lines(density((mu1 - mu0 )), col='steelblue')
dev.off()


