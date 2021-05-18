## ---------------------------
##
## Script name: chapter4.R
##
## Purpose of script:
##
## Exercises of chapter 4 from the book
##
## How to use:
##
## Just run it
##
## Author: Raphael Scherrer
##
## Date Created: 2021-05-17
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
library(splines)
library(patchwork)

theme_set(theme_classic())

#### Exercise 4E1 ####

# The first line is the likelihood

#### Exercise 4E2 ####

# The posterior has two parameters of the model: mu and sigma

#### Exercise 4E3 ####

# P(mu, sigma|y_i) = standardize( Product_i Normal(y_i|mu, sigma) *
# Normal(mu|0,10) * Exponential(sigma|1) )

#### Exercise 4E4 ####

# The second line is the linear model

#### Exercise 4E5 ####

# There are three parameters in the posterior: alpha, beta and sigma

#### Exercise 4M1 ####

# Simulate data from the prior
rnorm(1e4, mean = rnorm(1e3, mean = 0, sd = 10), sd = rexp(1e3, rate = 1))

#### Exercise 4M2 ####

# Quap formula:

# y ~ dnorm(mu , sigma),
# mu ~ dnorm(0, 0, 10),
# sigma ~ dexp(1)

#### Exercise 4M3 ####

# Mathematical model definition:

# y_i ~ Normal(mu_i, sigma)
# mu_i = alpha + beta * x_i
# alpha ~ Normal(0, 10)
# beta ~ Uniform(0, 1)
# sigma ~ Exponential(1)

#### Exercise 4M4 ####

# Regression of height onto year

# height_i ~ Normal(mu_i, sigma)
# mu_i = alpha + beta * year_i
# alpha ~ Normal(170, 10) (ppl are probably around 170cm to start with)
# beta ~ Normal(0, 1) (change through time could be + or - but probably small)
# sigma ~ Exponential(1) (couple of cm of deviation on avg)

#### Exercise 4M5 ####

# If I know now that students get taller through time I would change the prior
# on beta by:
# beta ~ Exponential(1)
# symbolizing that the students may get taller by a few cm on avg but may not
# get smaller

#### Exercise 4M6 ####

# If the variance is never more than 64cm then we could revise by sth like:
# sigma ~ Uniform(0, 8) (8 is the standard deviation, 64 = 8^2 is the variance)

#### Exercise 4M7 ####

# Refit model m4.3 from the chapter, but omit the mean weight xbar this time. Compare the
# new model’s posterior to that of the original model. In particular, look at the covariance among the
# parameters. What is different? Then compare the posterior predictions of both models.

# Load data
data(Howell1)
d <- Howell1 %>% filter(age >= 18)

# Define the average weight, x-bar
xbar <- mean(d$weight)

# Fit model
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma) ,
    mu <- a + b * (weight - xbar),
    a ~ dnorm(178 , 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data = d
)

m4.3
vcov(m4.3)

# Fit model without mean weight
m4.3.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma) ,
    mu <- a + b * weight,
    a ~ dnorm(178 , 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data = d
)

# Compare posteriors
m4.3
m4.3.1

# a seems underestimated

# Compare covariance matrices
vcov(m4.3)
vcov(m4.3.1)

# There is a large variance in a
sim4.3 <- t(sim(m4.3, d, n = 100))
sim4.3.1 <- t(sim(m4.3.1, d, n = 100))
colnames(sim4.3) <- colnames(sim4.3.1) <- paste0("sim", 1:100)

# Assemble the two simulated datasets into one and plot
map2(list(sim4.3, sim4.3.1), c("m4.3", "m4.3.1"), function(simdata, model) {

  as_tibble(cbind(d, simdata)) %>%
    mutate(model = model) %>%
    pivot_longer(sim1:sim100)

}) %>%
  do.call("rbind", .) %>%
  ggplot(aes(x = weight, y = value, color = model)) +
  geom_point(alpha = 0.05)

# Predictions of both models look similar though. This is because omitting xbar
# changes the origin, so the meaning of the intercept a. It makes sense that
# there is a larger variance in its estimation when weight is not measured
# relative to the mean weight, because then the intercept is located well
# outside of the domain of validity of the model, i.e. where most points are
# (no-one has weight zero).

#### Exercise 4M8 ####

# Load data
data(cherry_blossoms)
d <- cherry_blossoms
precis(d)
d <- d[complete.cases(d$doy),] # complete cases on doy

# Function to make a spline plot
blossom_plot <- function(num_knots, quap_list = NULL) {

  knot_list <- quantile(d$year, probs = seq(0, 1, length.out = num_knots))

  # Bases of the spline
  B <- bs(
    d$year,
    knots = knot_list[-c(1, num_knots)],
    degree = 3,
    intercept = TRUE
  )

  # Fit model
  if (is.null(quap_list)) {

    m4.7 <- quap(
      alist(
        D ~ dnorm(mu, sigma),
        mu <- a + B %*% w,
        a ~ dnorm(100, 10),
        w ~ dnorm(0, 10),
        sigma ~ dexp(1)
      ),
      data = list(D = d$doy, B = B),
      start = list(w = rep(0, ncol(B)))
    )

  } else {

    m4.7 <- quap(
      quap_list,
      data = list(D = d$doy, B = B),
      start = list(w = rep(0, ncol(B)))
    )

  }

  # Posterior samples
  post <- extract.samples(m4.7)
  w <- apply(post$w, 2, mean)

  # Posterior means for each year
  mu <- link(m4.7)

  # Interval of those means
  mu_PI <- apply(mu, 2, PI, 0.97)

  # Plot
  mu_PI <- as_tibble(t(mu_PI))
  colnames(mu_PI) <- c("q02", "q98")

  d %>%
    cbind(mu_PI) %>%
    as_tibble() %>%
    ggplot(aes(x = year, y = doy)) +
    geom_point(alpha = 0.3) +
    geom_ribbon(mapping = aes(xmax = year, ymin = q02, ymax = q98), alpha = 0.5) +
    xlab("Year") +
    ylab("Day of the year")

}

# Make the different plots
p1 <- blossom_plot(15) + ggtitle("15 knots")
p2 <- blossom_plot(20) + ggtitle("20 knots")
p3 <- blossom_plot(25) + ggtitle("25 knots")

# Show the plots
p1 / p2 / p3

ggsave("results/exercise_4M8-1.png", width = 4, height = 6, dpi = 300)

# The spline gets a bit wigglier

# Now with different standard deviations for the prior of w
p1 <- blossom_plot(
  15,
  quap_list = alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 1),
    sigma ~ dexp(1)
  )
) + ggtitle(parse(text = "sigma[w]==1"))

p2 <- blossom_plot(
  15,
  quap_list = alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  )
) + ggtitle(parse(text = "sigma[w]==10"))

p3 <- blossom_plot(
  15,
  quap_list = alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 20),
    sigma ~ dexp(1)
  )
) + ggtitle(parse(text = "sigma[w]==20"))

# Show the plots
p1 / p2 / p3
ggsave("results/exercise_4M8-2.png", width = 4, height = 6, dpi = 300)

# Does not change much
