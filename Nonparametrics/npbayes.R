devtools::install_github("StableMarkets/ChiRP",ref = 'fast-reg')
library(ChiRP)

library(BayesTree)
library(rstan)
library(gtools)

#### ------------------------ Simulate Data --------------------------------####
set.seed(1)

N = 500
X = seq(0, 2*pi, length.out = N)

X_01 = ((X - min(X)) / ( max(X) - min(X) ))

oscil1 = exp(-1*X)*cos(2*pi*X)
oscil2 = rev((exp(-1*X)*cos(2*pi*X)))

Y = (X_01)*oscil1 + (rev(X_01))*oscil2  + rnorm(N,0, .02)

test = ifelse(X>2.5 & X<3.75, 1, 0)


plot(X, Y, col=ifelse(test==1, 'red', 'blue'), pch=20)

d = data.frame(X=scale(X), Y=scale(Y))

#### ------------------------ Dirichlet Process ----------------------------####

d_train = d[test==0, ]
d_test = d[test==1, ]

set.seed(2)
res=fDPMix(d_train = d_train, formula = Y ~ X, d_test = d_test, 
           iter=1000, burnin=500,init_k = 10)


png('dp_oscil.png')
plot(d$X, d$Y, pch=20, 
     xlim=c(-2,2), ylim=c(-3,3), 
     col=ifelse(test==1,'pink','lightblue'),
     axes=F, xlab='X', ylab='Y')

### Plot posterior credible Band 
colfunc <- colorRampPalette(c("white", "skyblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

X_draw = c(d_train$X, d_test$X)
ord = order(X_draw)

X_draw = X_draw[ord]
res_stack = rbind(res$train, res$test)
res_stack = res_stack[ord,]

for(i in ci_perc){
  pci = apply(res_stack, 1, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  # polygon(c(X_draw,rev(X_draw)),c(pci[1,],rev(pci[2,])),
  #         col = colvec[as.character(i)], border = F)
  segments(X_draw,pci[1,], X_draw, pci[2,], col=colvec[as.character(i)], lwd=4 )
}
###

axis(1, -2:2, -2:2)
axis(2, -3:3, -3:3)

points(d$X, d$Y, pch=20, col=ifelse(test==1,'pink','gray'))
lines(d$X, rowMeans(res_stack), col='steelblue', pch=20)
dev.off()

#### ----------------------- BART Model         ----------------------------####

set.seed(3)
bart_res = bart(x.train = d_train$X, y.train = d_train$Y, x.test = d_test$X)

png('bart_oscil.png')
plot(d$X, d$Y, pch=20, 
     xlim=c(-2,2), ylim=c(-3,3), 
     col=ifelse(test==1,'pink','lightblue'),
     axes=F, xlab='X', ylab='Y')

### Plot posterior credible Band 
colfunc <- colorRampPalette(c("white", "skyblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

X_draw = c(d_train$X, d_test$X)
ord = order(X_draw)

X_draw = X_draw[ord]
res_stack = rbind(t(bart_res$yhat.train), t(bart_res$yhat.test))
res_stack = res_stack[ord,]

for(i in ci_perc){
  pci = apply(res_stack, 1, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  polygon(c(X_draw,rev(X_draw)),c(pci[1,],rev(pci[2,])),
          col = colvec[as.character(i)], border = F)
}
###

axis(1, -2:2, -2:2)
axis(2, -3:3, -3:3)

points(d$X, d$Y, pch=20, col=ifelse(test==1,'pink','gray'))
points(d$X, rowMeans(res_stack), col='steelblue', pch=20, type='l', lwd=1)
dev.off()

#### ----------------------- Gaussian Process   ----------------------------####
stan_data <- list(N1=nrow(d_train), 
                  x1 = d_train$X, y1=d_train$Y,
                  x2 = d_test$X, N2=length(d_test$X))

## https://mc-stan.org/docs/2_22/stan-users-guide/gaussian-processes-chapter.html
## GP code from Stan site that implements Guassian Process regression with 
## inference for hyperparamters.
## see above link for guidance on hyperprior calibration.
stan_mod = stan_model("gaussian_process_with_HPs.stan")

## In practice, should run for many more iterations and check convergence diagnostics
stan_res = sampling(stan_mod, data=stan_data, 
                    chains=1, iter=100, seed=100)

post_f_test=extract(stan_res, pars=c('f_test'))[[1]]
post_f_train=extract(stan_res, pars=c('f_train'))[[1]]

png('gp_oscil.png')

plot(d$X, d$Y, pch=20, 
     xlim=c(-2,2), ylim=c(-3,3), 
     col=ifelse(test==1,'pink','lightblue'),
     axes=F, xlab='X', ylab='Y')

### Plot posterior credible Band 
colfunc <- colorRampPalette(c("white", "skyblue"))
ci_perc = seq(.99,.01,-.01)
colvec = colfunc(length(ci_perc))
names(colvec) = ci_perc

X_draw = c(d_train$X, d_test$X)
ord = order(X_draw)

X_draw = X_draw[ord]
res_stack = rbind(t(post_f_train), t(post_f_test))
res_stack = res_stack[ord,]

for(i in ci_perc){
  pci = apply(res_stack, 1, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
  polygon(c(X_draw,rev(X_draw)),c(pci[1,],rev(pci[2,])),
          col = colvec[as.character(i)], border = F)
}
###

axis(1, -2:2, -2:2)
axis(2, -3:3, -3:3)

points(d$X, d$Y, pch=20, col=ifelse(test==1,'pink','gray'))
points(d$X, rowMeans(res_stack), col='steelblue', pch=20, type='l', lwd=1)
dev.off()
