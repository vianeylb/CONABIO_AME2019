## Markov-Switching Splines Stan Setup
## code partially taken from the splines case study found at: 
## https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

library("splines")
library("rstan")
library(tidyverse)

set.seed(17)

## Simulate data from the two components
xmat <- seq(from=-5, to=5, by=.1) # generating inputs
B <- t(bs(xmat, knots=seq(-5,5,1), degree=3, intercept = TRUE)) # creating the B-splines
TT <- length(xmat); num_basis <- nrow(B)
a0 <- c(0.2, 0.9) # intercept
a <- rbind(rnorm(num_basis, 0, 1), 
           rnorm(num_basis, -2, 0.5)) # coefficients of B-splines
Y_true <- rbind(as.vector(a0[1]*xmat + a[1,]%*%B),
                as.vector(a0[2]*xmat + a[2,]%*%B)) # generating the output
Y <- Y_true + rnorm(length(xmat),0,0.5) # adding noise

## scatter plot of observations
ggplot(data=tibble(x=rep(xmat, 2), y=c(Y[1,], Y[2,]), group=factor(rep(c(1,2), each=101))), 
       aes(x, y, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")

## scatter plot of observations with true curves
ggplot(data=tibble(x=rep(xmat, 2), 
                   y=c(Y[1,], Y[2,]), 
                       ytrue=c(Y_true[1,], Y_true[2,]), 
                   group=factor(rep(c(1,2), each=101))), 
       aes(x, y, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  geom_line(aes(x, ytrue, group=group), size=1) + 
  theme_minimal() + scale_color_discrete(name="State")


## Simulating from a Markov-switching spline model

tpm <- matrix(c(0.95, 0.05, 0.05, 0.95), 
              nrow=2, byrow=T)
init <- c(0.5, 0.5)

len <- 3000
states <- numeric(len)
obs <- numeric(len)

states[1] <- sample(size = 1, x=1:2, prob = init)
for(j in 2:len){
  states[j] = sample(size=1, x=1:2, prob=tpm[states[j-1],])
}

plot(ts(states))

xobs <- sample(size=len, x=xmat, replace = T, prob=rep(1/101, 101))
Bobs <- t(bs(xobs, knots=seq(-5,5,1), degree=3, intercept = TRUE))

for(j in 1:len){
  obs[j] = Y_true[states[j],which(xobs[j]==xmat)] + rnorm(1, 0, 0.5)
}

ggplot(data=tibble(xobs, obs, group=factor(states)), 
       aes(xobs, obs, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")

ggplot(data=tibble(x=1:len, obs, group=factor(states)), 
       aes(x, obs, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")



## Fitting the model 

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan.data.splines <- list(TT = len, y = obs, N=2, xmat=xobs, num_basis = num_basis, B=Bobs)
fit.splines <- stan("MarkovSwitchingSplines.stan", data=stan.data.splines, chains = 1, iter=1000)






