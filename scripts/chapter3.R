## ---------------------------
##
## Script name: chapter3.R
##
## Purpose of script:
##
## Exercises of chapter 3 from the book
##
## How to use:
##
## Just run it
##
## Author: Raphael Scherrer
##
## Date Created: 2021-05-11
##
## This script comes with no guarantee whatsoever.
##
## Copyright (c) Raphael Scherrer, 2021
##
## Find me on GitHub at https://github.com/rscherrer
##
## Email:
## r.scherrer@rug.nl
## raphael.scherrer@evobio.eu
## raph.rjfs@hotmail.fr
##
## ---------------------------

rm(list = ls())

library(tidyverse)
library(rethinking)

theme_set(theme_classic())

# Set-up copied from the book
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

#### Exercise 3E1 ####

# Proportion of the posterior sample below 0.2
length(which(samples < 0.2)) / length(samples)

#### Exercise 3E2 ####

# Proportion of the posterior sample above 0.8
length(which(samples > 0.8)) / length(samples)

#### Exercise 3E3 ####

# Proportion of the posterior sample between 0.2 and 0.8
length(which(samples >= 0.2 & samples <= 0.8)) / length(samples)

#### Exercise 3E4 ####

# 20% of the distribution is below...
quantile(samples, 0.2)

#### Exercise 3E5 ####

# 20% of the distribution is above...
quantile(samples, 0.8)

#### Exercise 3E6 ####

# Highest posterior density interval (narrowest interval containing 66%)
HPDI(samples, prob = 0.66)

#### Exercise 3E7 ####

# Percentile interval (centered interval containing 66%)
PI(samples, prob = 0.66)

#### Exercise 3M1 ####

# Construct a posterior for a new dataset

# Observed data
nwaters <- 8
ntosses <- 15

# Possible values of the parameter
pwaters <- seq(0, 1, 0.001)

# Flat prior
prior <- 1

# Likelihood of each under a binomial
likelihoods <- dbinom(nwaters, ntosses, prob = pwaters)

# Compute the posterior densities
posteriors <- likelihoods * prior
posteriors <- posteriors / sum(posteriors)

#### Exercise 3M2 ####

# Generate a posterior sample
samples <- sample(pwaters, size = 10000, prob = posteriors, replace = TRUE)

# 90% HPDI
HPDI(samples, 0.9)

#### Exercise 3M3 ####

# Construct a posterior predictive check for this model and data. This means
# simulate the distribution of samples, averaging over the posterior
# uncertainty in p. What is the probability of observing 8 water in 15 tosses?

# Simulate draws from a weighted range of probabilities
simulated <- rbinom(10000, size = ntosses, prob = samples)

# Probability to draw 8 in 15
length(which(simulated == 8)) / length(simulated)

#### Exercise 3M4 ####

# Using the posterior distribution constructed from the new (8/15) data,
# now calculate the probability of observing 6 water in 9 tosses.

simulated <- rbinom(10000, size = 9, prob = samples)
length(which(simulated == 6)) / length(simulated)

#### Exercise 3M5 ####

# New, biased prior
priors <- rep(1, length(pwaters))
priors[pwaters < 0.5] <- 0

# Compute the posterior density
posteriors <- likelihoods * priors
posteriors <- posteriors / sum(posteriors)

# Generate a posterior sample
samples <- sample(pwaters, size = 10000, prob = posteriors, replace = TRUE)

# 90% HPDI
HPDI(samples, 0.9)

# Simulate draws from a weighted range of probabilities
simulated <- rbinom(10000, size = ntosses, prob = samples)

# Probability to draw 8 in 15
length(which(simulated == 8)) / length(simulated)

# Probability to draw 6 in 9
simulated <- rbinom(10000, size = 9, prob = samples)
length(which(simulated == 6)) / length(simulated)

# The new prior makes both 8 in 15 and 6 in 9 more likely than before

# Compare with the true value for 8 in 15...
simulated <- rbinom(10000, size = 15, prob = 0.7)
length(which(simulated == 8)) / length(simulated)

# ... and for 6 in 9
simulated <- rbinom(10000, size = 9, prob = 0.7)
length(which(simulated == 6)) / length(simulated)

# Both priors give an overestimated number of 8 in 15, but the better one even
# more so.

#### Exercise 3M6 ####

# Function to get the mean span of the 99% PI of a posterior sample
get_mean_span <- function(
  ntosses, pwater = 0.7, ndraws = 10, nsamples = 10000
) {

  # ntosses is the number of tosses of the globe
  # pwater is the assumed true proportion of water
  # ndraws is the number of draws when sampling some data
  # nsamples is the number of draws when sampling from the posterior

  # Several draws
  nwaters <- rbinom(ndraws, ntosses, prob = pwater)

  # For each draw...
  spans <- purrr::map_dbl(nwaters, function(nwaters) {

    # Likelihood of each under a binomial
    likelihoods <- dbinom(nwaters, ntosses, prob = pwaters)

    # Compute the posterior density
    posteriors <- likelihoods * priors
    posteriors <- posteriors / sum(posteriors)

    # Generate a posterior sample
    samples <- sample(pwaters, size = nsamples, prob = posteriors, replace = TRUE)

    # 99% PI
    interval <- PI(samples, 0.99)

    # What is the span of that interval?
    diff(interval)

  })

  # Return the mean span
  return(mean(spans))

}

# The answer is somewhere between these two
get_mean_span(1000)
get_mean_span(3000)

# Refine our inference
tibble(ntosses = seq(1000, 3000, 10)) %>%
  mutate(span = purrr::map_dbl(ntosses, get_mean_span)) %>%
  ggplot(aes(x = ntosses, y = span)) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  xlab("Number of tosses") +
  ylab("99% PI mean span")

ggsave("results/exercise_3M6.png", width = 4, height = 2, dpi = 300)

# This was a computationally intensive approach because we have to average
# over lots of draws for each number of tosses. Plus, averaging over draws
# ignores the variability in spans for any given number of tosses. Our approach
# also assumes a specific true value, and I do not know if the estimate will
# change depending on that. Nevertheless, knowing that the Earth is ~ 70% water
# we can say that we need ~ 2200 tosses of the globe to get on average a
# posterior whose 99% PI is 0.05 in span

#### Exercise 3H1 ####

# Data from the book
data(homeworkch3)

birth1
birth2

# Observed data
nboys <- sum(birth1) + sum(birth2)
nbirths <- length(birth1) + length(birth2)

# Possible values for the prob. of a birth to be a boy
pboys <- seq(0, 1, 0.001)

# Flat prior
prior <- 1

# Likelihood of the observed data for each probability
likelihoods <- dbinom(nboys, nbirths, prob = pboys)

# Compute the posterior probabilities
posteriors <- likelihoods * prior
posteriors <- posteriors / sum(posteriors)

# What is the probability of having a boy that maximizes the posterior?
pboys[which.max(posteriors)]

#### Exercise 3H2 ####

# Sample from the posterior
samples <- sample(pboys, size = 10000, replace = TRUE, prob = posteriors)

# Various highest posterior density intervals
HPDI(samples, prob = c(50, 89, 97))

#### Exercise 3H3 ####

# Simulate many datasets from the posterior samples
simulated <- rbinom(10000, size = nbirths, prob = samples)

# Simulated and observed numbers of boys
tibble(nboys = simulated) %>%
  ggplot(aes(x = nboys)) +
  geom_density() +
  geom_vline(xintercept = nboys, linetype = 2) +
  xlab("Number of boys") +
  ylab("Density")

ggsave("results/exercise_3H3.png", width = 3, height = 2, dpi = 300)

#### Exercise 3H4 ####

# Draw many sets of 100 first births
simulated <- rbinom(10000, size = 100, prob = samples)

# Simulated and observed numbers of boys in the first birth
tibble(nboys = simulated) %>%
  ggplot(aes(x = nboys)) +
  geom_density() +
  geom_vline(xintercept = sum(birth1), linetype = 2) +
  xlab("Number of boys") +
  ylab("Density")

ggsave("results/exercise_3H4.png", width = 3, height = 2, dpi = 300)

# The model overestimates the number of first born boys

#### Exercise 3H5 ####

# Number of first-born girls
n_girl_first <- sum(!birth1)

# Simulate many sets of as many births and see how many are boys
simulated <- rbinom(10000, size = n_girl_first, prob = samples)

# Number of boys following girls
n_boy_after_girl <- sum(birth2[birth1 == 0])

# Simulated and observed numbers of boys following girls
tibble(nboys = simulated) %>%
  ggplot(aes(x = nboys)) +
  geom_density() +
  geom_vline(xintercept = n_boy_after_girl, linetype = 2) +
  xlab("Number of boys") +
  ylab("Density")

ggsave("results/exercise_3H5.png", width = 3, height = 2, dpi = 300)
