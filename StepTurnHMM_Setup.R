library(rstan)
library(moveHMM)
library(tidyverse)
library(readr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

rawhaggis <- read_csv("haggis.csv")
# derive step lengths and turning angles from locations
data <- prepData(rawhaggis, type="UTM")
# data <- prepData(data.frame(rawhaggis), type="UTM")
# plot(data)
#data <- as_tibble(prepData(data.frame(rawhaggis), type="UTM"))
#plot(data)

ggplot(data, aes(step)) + geom_histogram(fill="white", col="black") + 
  theme_minimal() + xlab("Step Length")

ggplot(data, aes(angle)) + geom_histogram(fill="white", col="black") + 
  theme_minimal() + xlab("Turning Angle")


# set NAs to out-of-range values
data$step[is.na(data$step)] <- -10
data$angle[is.na(data$angle)] <- -10
data$ID <- as.numeric(data$ID)

# full data set - no covariates
stan.data.nocovs <- list(T=nrow(data), ID=data$ID, steps=data$step, angles=data$angle, N=2, nCovs=0)

inits <- list(list(mu=c(1,5), sigma=c(1,5), xangle=c(-1,3), 
                   yangle=c(0,0), beta=matrix(c(-2,-2,0,0,0,0,0,0),nrow=2)),
              list(mu=c(1,5), sigma=c(1,5), xangle=c(-1,3), 
                   yangle=c(0,0), beta=matrix(c(-2,-2,0,0,0,0,0,0),nrow=2)))

fit <- stan(file="StepTurnHMM.stan", data=stan.data, init=inits,
            control=list(adapt_delta=0.9), chains=2)


# full data set - takes a bit to run
stan.data <- list(T=nrow(data), ID=data$ID, steps=data$step, angles=data$angle, N=2, nCovs=3,
                  covs=cbind(1, data$temp-10, data$slope-20, (data$slope-20)^2))

inits <- list(list(mu=c(1,5), sigma=c(1,5), xangle=c(-1,3), 
                   yangle=c(0,0), beta=matrix(c(-2,-2,0,0,0,0,0,0),nrow=2)),
              list(mu=c(1,5), sigma=c(1,5), xangle=c(-1,3), 
                   yangle=c(0,0), beta=matrix(c(-2,-2,0,0,0,0,0,0),nrow=2)))

fit <- stan(file="StepTurnHMM.stan", data=stan.data, init=inits,
            control=list(adapt_delta=0.9), chains=2)
