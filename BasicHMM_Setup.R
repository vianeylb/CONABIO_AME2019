## Setup for a basic HMM
library(tidyverse)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

N <- 2
init <- rep(1/N, N)

tpm <- matrix(c(0.8, 0.2, 0.3, 0.7), ## row corresponds to state 3
              byrow=T, nrow=2)

len <- 3000
mu <- c(5, 10)
sd <- c(2, 1)

states <- numeric(len)
obs <- numeric(len)

## sample function 
states[1] <- sample(size=1, x = 1:N, prob = init)

for(j in 2:len){
  states[j] <- sample(size=1, x = 1:N, prob = tpm[states[j-1],])
}

obs <- rnorm(len, mean=mu[states], sd = sd[states])

## scatter plot of observations
ggplot(data=tibble(x=1:len, y=obs, group=factor(states)), 
       aes(x, y, color=group)) + geom_point() + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="State")

ggplot(data=tibble(x=1:len, y=obs, group=factor(states)), 
       aes(x, y)) + geom_point(aes(color=group)) + xlab("X") + ylab("Y") + 
  theme_minimal() + scale_color_discrete(name="States") + 
  geom_path(aes(x=1:len, y=obs), alpha=0.3)


state1dur <- dgeom(x = 0:100, prob = 0.2)
state2dur <- dgeom(x=  0:100, prob = 0.3)

stateDur <- tibble(x=rep(1:101, 2),
                   dur = c(state1dur, state2dur), 
                   group=factor(rep(1:2, each=101)))

ggplot(stateDur, aes(x, dur, color=group)) + 
  geom_point(size=2) + theme_minimal()

ggplot(stateDur, aes(x, dur, color=group)) + 
  geom_point(size=3) + theme_minimal() + 
  xlim(0, 25) + geom_path(size=2, alpha=0.4) + 
  ggtitle("State Duration Distributions")

## Fitting the model

stan.data.bhmm <- list(N = N, y=obs, T=len)
fit.bhmm <- stan(file="BasicHMM.stan", data=stan.data.bhmm, chains=3, 
                 iter=2000, warmup = 1000)

plot(fit.bhmm, pars=c("mu"), plotfun="hist")
plot(fit.bhmm, pars=c("sigma"), plotfun="hist")
plot(fit.bhmm, pars=c("tpm"), plotfun="hist")

#install.packages("shinystan")
#library(shinystan)
#launch_shinystan(fit.bhmm)



