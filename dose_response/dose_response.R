###--------------------------------------------------------------------------###
### Author: Arman Oganisian
### Conduct dose-response analysis using autoregressive prior.
###--------------------------------------------------------------------------###

## Load Packages
library(rstan)
library(LaplacesDemon)
library(latex2exp)
set.seed(1)

####------------------------ Simulate Data ---------------------------------####
K = 10 # number of dose levels
N = 100 # sample size
warmup = 500
iter = 1000
n_draws = iter - warmup


#simulate standard normal confounder (L), dose (A), and outcome (Y)
L = rnorm(n = N) 

A = numeric(length = N)

## higher doses are less likely to be assigned seq(1,-1,length.out = K)
## L influences assignment
D = 0:(K-1)
for(i in 1:N){
  A[i] = sample(D, 1, prob = invlogit( 1 - (2/9)*(D) - .5*(D)*L[i] + L[i]  ) )
}

## higher doses are more sparsely populated
table(A)/N

## Gaussian outcome
# pnorm(A-5) implies strictly increasing dose effect that initially is steep, 
## has an inflection point at A=5, and begins having diminishing effects on Y
Y = rnorm(n = N, mean =  5*pnorm( A - 5) - 5*L, sd = 2 )

# P=2 dimensions of model matrix
A_mat = model.matrix(Y ~ -1 + as.factor(A) )
X = model.matrix(Y ~ L )

stan_data = list(Y=Y, L=X, A = A_mat[,-1],
                 N=N, num_A_levels = K-1, P = ncol(X) )

####------------------------ Sample Posterior    ---------------------------####
DR_model = stan_model(file = "DR_model.stan")

stan_res = sampling(DR_model, data = stan_data, 
                    warmup = warmup, iter = iter, chains=1, seed=1)

Psi_draws = extract(stan_res, pars='Psi')[[1]]


####------------------- Compute Frequentist Estimates   --------------------####

freq_reg = lm(Y ~  L + as.factor(A))
summary(freq_reg)

Psi_freq = numeric(length = K-1)

Psi_freq[1] = freq_reg$coefficients[3]
Psi_freq[2:(K-1)] = freq_reg$coefficients[4:(K-1+2)] - freq_reg$coefficients[3:(K-2+2)]

####-------------------         Plot Results            --------------------####
dose = 1:(K-1)
true_Psi = 5*pnorm(dose-5) - 5*pnorm((dose-1)-5)

png("dose_response_curve.png",width = 600, height = 500)
plot( dose, colMeans(Psi_draws), ylim=c(-4,4), 
      col='blue', pch=20, ylab=TeX("$\\Psi(k)$"), type='o',axes=F)

axis_labs = paste0(dose, " \n (n=",table(A[A!=0]),")")
axis(side = 1, seq(1:(K-1)), labels = axis_labs,tick = T, padj = .5)
axis(side = 2, seq(-4,4,1), labels = seq(-4,4,1), tick = T)

### Plot posterior credible Band 
colfunc <- colorRampPalette(c("white", "lightblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

for(i in ci_perc){
  pci = apply(Psi_draws, 2, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  polygon(c(dose,rev(dose)),c(pci[1,],rev(pci[2,])),
          col = colvec[as.character(i)], border = FALSE)
}
###

lines(dose, colMeans(Psi_draws), col='steelblue', type='o', pch=20)

points(dose, Psi_freq, col='black', pch=20, type='o')
points(dose, true_Psi, col='red', pch=20, type='o' )
legend('bottomleft', 
       legend = c('Posterior Mean/Credible Band', 'MLE', 'True Curve'),
       col = c('steelblue', 'black', 'red'), pch=c(20,20,20), lty=c(1,1,1),
       bty='n')
dev.off()