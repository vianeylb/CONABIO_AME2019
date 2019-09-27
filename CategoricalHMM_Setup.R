## Categorical HMM Setup

## Setup for a basic HMM
library(tidyverse)
library(rstan)
library(MCMCpack)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

N <- 2
init <- rep(1/N, N)

tpm <- matrix(c(0.99, 0.01, 0.03, 0.97), ## row corresponds to state 3
              byrow=T, nrow=2)

len <- 3000
theta <- matrix(c(0.1, 0.2, 0.7, 
               0.95, 0.04, 0.01), 
               byrow=T, nrow=2)

states <- numeric(len)
obs <- numeric(len)

## sample function 
states[1] <- sample(size=1, x = 1:N, prob = init)

for(j in 2:len){
  states[j] <- sample(size=1, x = 1:N, prob = tpm[states[j-1],])
}

for(j in 1:len) obs[j] <- sample(x = 1:3, size=1, prob = theta[states[j],])

## scatter plot of observations
ggplot(data=tibble(x=1:len, y=obs, group=factor(states)), 
       aes(x, y, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")

ggplot(data=tibble(x=1:len, y=obs, group=factor(states)), 
       aes(x, y)) + geom_point(aes(color=group)) + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="States") + 
  geom_path(aes(x=1:len, y=obs), alpha=0.3)


state1dur <- dgeom(x = 0:100, prob = 0.01)
state2dur <- dgeom(x=  0:100, prob = 0.03)

stateDur <- tibble(x=rep(1:101, 2),
                   dur = c(state1dur, state2dur), 
                   group=factor(rep(1:2, each=101)))


ggplot(stateDur, aes(x, dur, color=group)) + 
  geom_point(size=3) + theme_minimal() + 
  xlim(0, 25) + geom_path(size=2, alpha=0.4) +
  scale_color_discrete(name="State") + 
  ylab("Duration") + 
  ggtitle("State Duration Distributions")

## Fitting the model

as <- matrix(c(1, 1, 5, 
               100, 1, 1), 
             nrow=2, byrow=T)

stan.data.cathmm <- list(N = N, y=obs, T=len, cats=3, as=as)
fit.cathmm <- stan(file="CategoricalHMM.stan", data=stan.data.cathmm, chains=3)

