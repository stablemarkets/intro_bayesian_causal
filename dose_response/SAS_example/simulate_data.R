###--------------------------------------------------------------------------###
### Author: Arman Oganisian
### Simulate dose-response analysis data for import to SAS
###--------------------------------------------------------------------------###

## Load Packages
library(LaplacesDemon)
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

d = data.frame(Y=Y, A=A, L=L )

write.csv(d, "simulatd_dose_data.csv")
