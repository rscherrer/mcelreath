## ---------------------------
##
## Script name: chapter2.R
##
## Purpose of script:
##
## Exercises of chapter 2 from the book
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

# Exercises from chapter 2

rm(list = ls())

library(tidyverse)

theme_set(theme_classic())

#### Exercise 2M1 ####

# Function to get the posterior density from an observed sample
get_posterior <- function(nwater, n) {

  # nwater: number of observed W
  # n: sample size

  # Proportions of water cover
  pwaters <- seq(0, 1, 0.001)

  # Uniform prior
  prior <- 1

  # Compute the likelihood given binomial process
  likelihoods <- dbinom(nwater, n, prob = pwaters)

  # Compute the posterior
  posteriors <- prior * likelihoods

  # Standardize the posterior
  posteriors <- posteriors / sum(posteriors)

  return(tibble(pwater = pwaters, posterior = posteriors))

}

# Plot the posterior distributions for a few observed samples
tibble(
  nwater = c(3, 3, 5),
  n = c(3, 4, 7)
) %>%
  mutate(
    case = paste("W =", nwater, "out of", n),
    posterior = map2(nwater, n, get_posterior)
  ) %>%
  unnest(cols = c(posterior)) %>%
  ggplot(aes(x = pwater, y = posterior)) +
  geom_line() +
  facet_wrap(. ~ case) +
  xlab("Proportion of water") +
  ylab("Posterior density")

ggsave("results/exercise_2M1.png", width = 6, height = 2, dpi = 300)

#### Exercise 2M2 ####

# Now with a truncated prior distribution
get_posterior2 <- function(nwater, n) {

  # nwater: number of observed W
  # n: sample size

  # Proportions of water cover
  pwaters <- seq(0, 1, 0.001)

  # Uniform prior
  priors <- rep(1, length(pwaters))
  priors[pwaters < 0.5] <- 0

  # Compute the likelihood given binomial process
  likelihoods <- dbinom(nwater, n, prob = pwaters)

  # Compute the posterior
  posteriors <- priors * likelihoods

  # Standardize the posterior
  posteriors <- posteriors / sum(posteriors)

  return(tibble(pwater = pwaters, posterior = posteriors))

}

# Plot the posterior distributions for a few observed samples
tibble(
  nwater = c(3, 3, 5),
  n = c(3, 4, 7)
) %>%
  mutate(
    case = paste("W =", nwater, "out of", n),
    posterior = map2(nwater, n, get_posterior2)
  ) %>%
  unnest(cols = c(posterior)) %>%
  ggplot(aes(x = pwater, y = posterior)) +
  geom_line() +
  facet_wrap(. ~ case) +
  xlab("Proportion of water") +
  ylab("Posterior density")

ggsave("results/exercise_2M2.png", width = 6, height = 2, dpi = 300)

#### Exercise 2M3 ####

# Proportion of water on Earth and Mars
pwaters <- c(earth = 0.7, mars = 0)

# Probability of samping "land" in one draw under the two models
likelihoods <- dbinom(0, 1, prob = pwaters)

# Equal priors
priors <- c(earth = 0.5, mars = 0.5)

# Compute the posterior
posteriors <- likelihoods * priors
posteriors / sum(posteriors)

# Indeed 0.23 for Earth



#### Exercise 2M4 ####

# Priors once we know a black side has been drawn
priors <- c(bb = 1/2, bw = 1/2, ww = 0)

# Probabilities of black sides for each card
pblacks <- c(bb = 1, bw = 0.5, ww = 0)

# Probabilities of drawing one black side for each card (= pblacks cause 1 draw)
likelihoods <- dbinom(1, 1, prob = pblacks)

# Probabilities to draw a second black given we already have one
posteriors <- likelihoods * priors
posteriors / sum(posteriors)

# The probability that the other side of the black-sided card on the table
# is also black is 2/3

#### Exercise 2M5 ####

# Same but with one extra bb card in the deck

# Priors once we know a black side has been drawn
priors <- c(bb = 2/3, bw = 1/3, ww = 0)

# Probabilities of black sides for each card
pblacks <- c(bb = 1, bw = 0.5, ww = 0)

# Probabilities of drawing one black side for each card (= pblacks cause 1 draw)
likelihoods <- dbinom(1, 1, prob = pblacks)

# Probabilities to draw a second black given we already have one
posteriors <- likelihoods * priors
posteriors / sum(posteriors)

# Now the probability is 80%

#### Exercise 2H1 ####

# T = having one twin
# T2 = having another pair of twins
# Species A or B
#
# P(T2|T) prob that the female has another pair of twins given already one
# = P(T|A) P(A|T1) + P(T|B) P(B|T1)
# in other words, prob that female has twins given she is A, times the prob that
# she is A given she had already one twin, plus the same thing for B
#
# where P(T|A) = 0.1 and P(T|B) = 0.2
#
# and P(A|T1) and P(B|T1) are posterior probabilities, so:
# P(A|T1) ~ P(T|A) * P(A) with equal priors = 0.1 * 1/2
# P(B|T1) ~ P(T|B) * P(B) = 0.2 * 1/2
# After standardization: 1/3 and 2/3
#
# Plugging these posteriors in the equation above gives 1 chance out of 6 that
# the female has another pair of twins

### Exercise 2H2 ####

# Prob. that panda is A given it had twins already is P(A|T1) = 1/3 (see above)

#### Exercise 2H3 ####

# T: number of twin births
# N: number of births
# P(A|T = 1, N = 2): prob. A given one twin event out of two births
# prop. to P(T = 1, N = 2|A) * P(A)
# where prior P(A) = P(A|T = 1, N = 1), i.e. posterior from after observing
# the first birth = 1/3
# and where P(T = 1, N = 2|A) = 0.1 * (1 - 0.1) = 0.09
# To normalize we need P(T = 1, N = 2|B) * P(B) = 0.2 * (1 - 0.2) * 2/3 = 0.107
# So normalizing constant is 0.09 * 1/3 + 0.107 = 0.137
# And P(A|T = 1, N = 2) = 0.09 * 1/3 / 0.137 = 0.219

#### Exercise 2H4 ####

# Ignoring births, compute P(A|test says A)

p_test_A_given_is_A <- 0.8
p_test_A_given_is_B <- 1 - 0.65
p_is_A <- p_is_B <- 1/2
p_is_A_given_test_A <- p_test_A_given_is_A * p_is_A
p_is_B_given_test_A <- p_test_A_given_is_B * p_is_B
norm <- p_is_A_given_test_A + p_is_B_given_test_A
p_is_A_given_test_A <- p_is_A_given_test_A / norm
p_is_A_given_test_A

# Now taking into account that one birth event out of two was twins

# Data:
# Test says A
# One twin event out of two births

# These probs we know
p_test_A_given_is_A <- 0.8
p_test_A_given_is_B <- 0.35

# Priors from what we know from the previous exercise based on births
p_is_A <- 0.219
p_is_B <- 1 - p_is_A

# Posterior derived by adding info from the test
p_is_A_given_data <- p_test_A_given_is_A * p_is_A
p_is_B_given_data <- p_test_A_given_is_B * p_is_B

# Normalizing factor
norm <- p_is_A_given_data + p_is_B_given_data

# Probability is...
p_is_A_given_data / norm
