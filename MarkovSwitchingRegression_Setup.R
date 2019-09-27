## Markov-Switching Regression Stan Setup

library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(17)

## Simulate data from the two components
xmat <- seq(from=-5, to=5, by=.1) # generating inputs
betas <- matrix(c(-1, 4, 0.5, -1), byrow=T, nrow=2)
means <- betas%*%t(cbind(1, xmat))
sigmas <- c(1, 3)
  
Y <- rbind(means + rnorm(101, 0, sigmas[1]), 
           means + rnorm(101, 0, sigmas[2])) 

## scatter plot of observations
ggplot(data=tibble(x=rep(xmat, 2), y=c(Y[1,], Y[2,]), group=factor(rep(c(1,2), each=101))), 
       aes(x, y, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")

## scatter plot of observations with true curves
ggplot(data=tibble(x=rep(xmat, 2), 
                   y=c(Y[1,], Y[2,]), 
                   means=c(means[1,], means[2,]), 
                   group=factor(rep(c(1,2), each=101))), 
       aes(x, y, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  geom_line(aes(x, means, group=group), size=1) + 
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

for(j in 1:len){
  obs[j] = means[states[j],which(xobs[j]==xmat)] + rnorm(1, 0, sigmas[states[j]])
}

ggplot(data=tibble(xobs, obs, group=factor(states)), 
       aes(xobs, obs, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")

ggplot(data=tibble(x=1:len, obs, group=factor(states)), 
       aes(x, obs, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")



## Fitting the model 
stan.data.regression <- list(TT = len, y = obs, N=2, ncovs=1, xmat = cbind(1, xobs))
fit.regression <- stan("MarkovSwitchingRegression.stan", data=stan.data.regression, chains = 3)



