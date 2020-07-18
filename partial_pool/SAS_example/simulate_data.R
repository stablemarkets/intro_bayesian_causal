###--------------------------------------------------------------------------###
### Author: Arman Oganisian
### Simulate partial pooling analysis data to import to SAS
###--------------------------------------------------------------------------###

## Load Packages
library(LaplacesDemon)
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

write.csv( data.frame(Y=Y, A=A, W=W, V=V), file = 'simulated_data.csv')
