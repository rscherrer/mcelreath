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

data(Howell1)
data(cherry_blossoms)

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


#### Exercise 4H1 ####

# Load the Kalahari data
data(Howell1)
data <- Howell1

# Useful values
mean_weight <- mean(data$weight)
sd_weight <- sd(data$weight)

# Add standardized and higher order predictors
data <- data %>%
  mutate(
    weight_s = (weight - mean_weight) / sd_weight,
    weight_s2 = weight_s^2,
    weight_s3 = weight_s^3
  )

# Estimate the posterior of a cubic regression
mod <- quap(
  flist = alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1 * weight_s + b2 * weight_s2 + b3 * weight_s3,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = data
)

# New weight entries
new_weights <- c(46.95, 43.72, 64.78, 32.59, 54.63)
new_weights_s <- (new_weights - mean_weight) / sd_weight

# Assemble in a data frame
new_data <- tibble(
  weight_s = new_weights_s,
  weight_s2 = weight_s^2,
  weight_s3 = weight_s^3
)

# Sample expected heights for the new weights from the posterior
pred <- link(mod, new_data)

# Get the mean and percentile interval for the expected height
new_data %>%
  mutate(mean_height = apply(pred, 2, mean)) %>%
  cbind(t(apply(pred, 2, PI, prob = 0.89)))

#### Exercise 4H2 ####

# Filter adults out from the Kalahari data
data <- Howell1
data <- data %>% filter(age < 18)

# Useful values
mean_weight <- mean(data$weight)
sd_weight <- sd(data$weight)

# Standardize weights
data <- data %>% mutate(weight_s = (weight - mean_weight) / sd_weight)

# Fit a linear regression
mod <- quap(
  flist = alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * weight_s,
    a ~ dnorm(100, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data = data
)

# Extract maximum a posteriori
map_a <- coef(mod)["a"]
map_b <- coef(mod)["b"]

# For an increase of 10kg the model predicts an increase in height of (in cm)...
10 * map_b / sd_weight

# Compute interval for the posterior of the mean
post_mean_heights <- as_tibble(t(apply(link(mod), 2, PI, prob = 0.89))) %>%
  rename(pi05_mean_height = "5%", pi94_mean_height = "94%")

# Compute interval for predicted heights
pred_heights <- as_tibble(t(apply(sim(mod), 2, PI, prob = 0.89))) %>%
  rename(pi05_pred_height = "5%", pi94_pred_height = "94%")

# Assemble all in one data frame
data <- data %>%
  cbind(post_mean_heights, pred_heights) %>%
  as_tibble()

# Plot all
data %>%
  ggplot(aes(x = weight, y = height)) +

  # The raw data
  geom_point() +

  # MAP regression line
  geom_abline(
    intercept = map_a - mean_weight * map_b / sd_weight,
    slope = map_b / sd_weight
  ) +

  # 89% PI of the posterior of mean heights as a dark shade
  geom_ribbon(
    aes(
      xmin = weight,
      xmax = weight,
      ymin = pi05_mean_height,
      ymax = pi94_mean_height
    ),
    alpha = 0.7
  ) +

  # 89% PI of simulated (predicted) heights as a lighter shade
  geom_ribbon(
    aes(
      xmin = weight,
      xmax = weight,
      ymin = pi05_pred_height,
      ymax = pi94_pred_height
    ),
    alpha = 0.2
  ) +
  xlab("Weight (kg)") +
  ylab("Height (cm)")

ggsave("results/exercise_4H2.png", width = 4, height = 3, dpi = 300)

# The model seems wrong in that the relationship between weight and height is
# not a straight line. Probably a quadratic or cubic regression would do a
# better job

#### Exercise 4H3 ####

# Load all the data
data <- Howell1

# Log-transform the weights
data <- data %>% mutate(log_weight = log(weight))

# Useful values
mean_log_weight <- mean(data$log_weight)
sd_log_weight <- sd(data$log_weight)

# Standardize the predictor
data <- data %>%
  mutate(log_weight_s = (log_weight - mean_log_weight) / sd_log_weight)

# Fit a regression
mod <- quap(
  flist = alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * log_weight_s,
    a ~ dnorm(150, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ),
  data = data
)

# Interpretation:
# a is the mean height of the population
# b is the increase in height predicted by the model when the standardized
# log-weight increases by one unit
# sigma is the standard deviation of the modelled normal distribution of heights
# around their predicted mean

# Compute interval for the posterior of the mean
post_mean_heights <- as_tibble(t(apply(link(mod), 2, PI, prob = 0.97))) %>%
  rename(pi05_mean_height = "2%", pi94_mean_height = "98%")

# Compute interval for predicted heights
pred_heights <- as_tibble(t(apply(sim(mod), 2, PI, prob = 0.97))) %>%
  rename(pi05_pred_height = "2%", pi94_pred_height = "98%")

# Assemble all in one data frame
data <- data %>%
  cbind(post_mean_heights, pred_heights) %>%
  as_tibble()

# Extract MAP coefficients
map_a <- coef(mod)[["a"]]
map_b <- coef(mod)[["b"]]

# Plot all (log-scale)
P1 <- data %>%
  ggplot(aes(x = log_weight, y = height)) +

  # The raw data
  geom_point() +

  # MAP regression line
  geom_abline(
    intercept = map_a - mean_log_weight * map_b / sd_log_weight,
    slope = map_b / sd_log_weight
  ) +

  # 89% PI of the posterior of mean heights as a dark shade
  geom_ribbon(
    aes(
      xmin = log_weight,
      xmax = log_weight,
      ymin = pi05_mean_height,
      ymax = pi94_mean_height
    ),
    alpha = 0.7
  ) +

  # 89% PI of simulated (predicted) heights as a lighter shade
  geom_ribbon(
    aes(
      xmin = log_weight,
      xmax = log_weight,
      ymin = pi05_pred_height,
      ymax = pi94_pred_height
    ),
    alpha = 0.2
  ) +
  xlab("Log-weight (log-kg)") +
  ylab("Height (cm)")

# Add a column for MAP mean height
data <- data %>%
  mutate(
    map_mean_height = map_a + map_b * log_weight_s
  )

# Plot all (normal scale)
P2 <- data %>%
  ggplot(aes(x = weight, y = height)) +

  # The raw data
  geom_point() +

  # MAP nonlinear regression line
  geom_line(aes(y = map_mean_height)) +

  # 89% PI of the posterior of mean heights as a dark shade
  geom_ribbon(
    aes(
      xmin = weight,
      xmax = weight,
      ymin = pi05_mean_height,
      ymax = pi94_mean_height
    ),
    alpha = 0.7
  ) +

  # 89% PI of simulated (predicted) heights as a lighter shade
  geom_ribbon(
    aes(
      xmin = weight,
      xmax = weight,
      ymin = pi05_pred_height,
      ymax = pi94_pred_height
    ),
    alpha = 0.2
  ) +
  xlab("Weight (kg)") +
  ylab("Height (cm)")

P1 / P2

ggsave("results/exercise_4H3.png", width = 4, height = 6, dpi = 300)

#### Exercise 4H4 ####

# Load data
data <- Howell1

# Standardize weights
data <- data %>%
  mutate(weight_s = (weight - mean(weight)) / sd(weight))

# Function to simulate height data from the prior
sim_from_prior <- function(
  x, n = 1e3, n_prior = 1e3, mean_a = 178, sd_a = 20, mean_b1 = 0, sd_b1 = 1,
  mean_b2 = 0, sd_b2 = 1, min_sigma = 0, max_sigma = 50
) {

  # x is a standardized weight value
  # n is the number of height values to simulate
  # n_prior is the number of estimates to sample from the prior distributions

  # Sample prior distributions
  a <- rnorm(n_prior, mean_a, sd_a)
  b1 <- rlnorm(n_prior, mean_b1, sd_b1)
  b2 <- rnorm(n_prior, mean_b2, sd_b2)
  sigma <- runif(n_prior, min_sigma, max_sigma)

  # Simulate heights
  rnorm(n = n, mean = a + b1 * x + b2 * x^2, sd = sigma)

}

# Predict 100 prior beliefs on heights for each of the weights in the dataset
pred_data <- data %>%
  as_tibble() %>%
  mutate(
    pred_height = purrr::map(weight_s, sim_from_prior, n = 100, n_prior = 1000)
  ) %>%
  unnest(cols = c(pred_height))

# Plot the data as predicted from the prior
P1 <- pred_data %>%
  ggplot(aes(x = weight, y = pred_height)) +
  geom_point(alpha = 0.2) +
  geom_point(mapping = aes(y = height), color = "red") +
  xlab("Weight (kg)") +
  ylab("Height (cm)") +
  ggtitle("Prior expectations versus real data")

# Another go with more reasonable priors
P2 <- data %>%
  as_tibble() %>%
  mutate(
    pred_height = purrr::map(
      weight_s, sim_from_prior, n = 100, n_prior = 1000, max_sigma = 5,
      sd_a = 10, mean_a = 150, mean_b1 = 0.5, sd_b1 = 0.1
    )
  ) %>%
  unnest(cols = c(pred_height)) %>%
  ggplot(aes(x = weight, y = pred_height)) +
  geom_point(alpha = 0.2) +
  geom_point(mapping = aes(y = height), color = "red") +
  xlab("Weight (kg)") +
  ylab("Height (cm)") +
  ggtitle("More reasonable priors")

P1 / P2

ggsave("results/exercise_4H4.png", width = 4, height = 6, dpi = 300)


#### Exercise 4H5 ####

# Load data
data <- as_tibble(cherry_blossoms)

# Remove missing data
data <- data %>% drop_na()

# Eyeball
data %>%
  ggplot(aes(x = temp, y = doy)) +
  geom_point()

# Useful values
mean_temp <- mean(data$temp)
sd_temp <- sd(data$temp)

# Standardize temperature
data <- data %>%
  mutate(temp_s = (temp - mean_temp) / sd_temp)

# Fit a linear regression
mod_linear <- quap(
  flist = alist(
    doy ~ dnorm(mu, sigma),
    mu <- a + b * temp_s,
    a ~ dnorm(100, 10),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 30)
  ),
  data = data
)

# Fit a quadratic regression
mod_polyn <- quap(
  flist = alist(
    doy ~ dnorm(mu, sigma),
    mu <- a + b1 * temp_s + b2 * temp_s^2,
    a ~ dnorm(100, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 1),
    sigma ~ dunif(0, 30)
  ),
  data = data
)

# Set-up knots for spline regression
num_knots <- 4
knots <- quantile(data$temp_s, probs = seq(0, 1, length.out = num_knots))

# Bases of the spline
bases <- bs(
  data$temp_s,
  knots = knots[-c(1, num_knots)],
  degree = 3,
  intercept = TRUE
)

# Fit a spline regression
mod_splines <- quap(
  alist(
    doy ~ dnorm(mu, sigma),
    mu <- a + base %*% w,
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ),
  data = list(doy = data$doy, base = bases),
  start = list(w = rep(0, ncol(bases)))
)

# Function to extract summary on the posterior of the expected response for each
# value of the predictor
extract_mu <- function(mod) {

  # Extract summaries of the posterior on the mean
  post_mu <- link(mod)
  means <- apply(post_mu, 2, mean)
  pis <- apply(post_mu, 2, PI, prob = 0.89)

  # Assemble in a tibble and return
  cbind(means, t(pis)) %>%
    as_tibble() %>%
    rename(mean_mu = "means", pi05_mu = "5%", pi94_mu = "94%")

}

# Function to extract intervals of the predicted response
predict_response <- function(mod) {

  # Extract predictions from the model
  preds <- sim(mod)
  means <- apply(preds, 2, mean)
  pis <- apply(preds, 2, PI, prob = 0.89)

  # Assemble in a tibble and return
  cbind(means, t(pis)) %>%
    as_tibble() %>%
    rename(mean_pred = "means", pi05_pred = "5%", pi94_pred = "94%")

}

# Plot each model with the raw data
plots <- purrr::map2(

  # Models
  list(mod_linear, mod_polyn, mod_splines),

  # Respective plot titles
  c(
    "Linear regression",
    "Quadratic regression",
    "Spline regression (4 knots of degree 3)"
  ),

  # Function to make each plot
  function(mod, title) {

    data %>%
      cbind(

        # Attach the posterior means and the predicted response
        extract_mu(mod),
        predict_response(mod)

      ) %>%
      ggplot(aes(x = temp, y = doy)) +
      geom_point() +

      # Regression line
      geom_line(mapping = aes(y = mean_mu)) +

      # Shade representing the posterior distribution of the mean
      geom_ribbon(
        mapping = aes(
          xmin = temp, xmax = temp,
          ymin = pi05_mu, ymax = pi94_mu
        ),
        alpha = 0.7
      ) +

      # Shade representing the distribution of the prediction
      geom_ribbon(
        mapping = aes(
          xmin = temp, xmax = temp,
          ymin = pi05_pred, ymax = pi94_pred
        ),
        alpha = 0.3
      ) +
      xlab("Temperature (C)") +
      ylab("Day of year") +
      ggtitle(title)

  }
)

# Assemble the plots
plots[[1]] / plots[[2]] / plots[[3]]

# Save them
ggsave("results/exercise_4H5.png", width = 4, height = 7, dpi = 300)

#### Exercise 4H6 ####

# Load the data
data <- as_tibble(cherry_blossoms)
data <- data %>% drop_na()

# Setup the knots
num_knots <- 15
knots <- quantile(data$year, probs = seq(0, 1, length.out = num_knots))

# Bases of the spline
bases <- bs(
  data$year,
  knots = knots[-c(1, num_knots)],
  degree = 3,
  intercept = TRUE
)

# Number of predictions per value of the predictor
n <- 100

# Function to predict and plot data from the prior given parameters for weights
predict_from_prior <- function(mean_w, sd_w) {

  # mean_w and sd_w are the mean and sd of the normal prior distribution of
  # knot weights in the spline regression model

  # Predict data from the prior
  pred_data <- data %>%
    mutate(
      pred_doy = purrr::map(seq(n()), function(i) {

        # Simulate n predictions from the prior
        map_dbl(seq(n), function(j) {

          # Sample model parameters from the prior
          a <- rnorm(1, 100, 10)
          w <- rnorm(ncol(bases), mean_w, sd_w)
          sigma <- rexp(1, 1)
          mu <- a + bases %*% w

          # Predict a value
          rnorm(1, mu[i], sigma)

        })

      })
    ) %>%
    unnest(cols = c(pred_doy))

  # Plot the prior predictions and the real data
  pred_data %>%
    ggplot(aes(x = year, y = pred_doy)) +
    geom_point(alpha = 0.1) +
    geom_point(mapping = aes(y = doy), color = "red") +
    xlab("Year") +
    ylab("Day of year")

}

# Try out different values
P1 <- predict_from_prior(0, 1) + ggtitle("Prior weights: mean = 0, s. d. = 1")
P2 <- predict_from_prior(0, 10) + ggtitle("Prior weights: mean = 0, s. d. = 10")
P3 <- predict_from_prior(5, 1) + ggtitle("Prior weights: mean = 5, s. d. = 1")
P4 <- predict_from_prior(5, 10) + ggtitle("Prior weights: mean = 5, s. d. = 10")

# Assemble the plots
(P1 | P2) / (P3 | P4)

ggsave("results/exercise_4H6.png", width = 8, height = 6, dpi = 300)

#### Exercise 4H7 ####

# Load the data
data <- as_tibble(cherry_blossoms)
data <- data %>% drop_na()

# Set-up knots for spline regression
num_knots <- 15
knots <- quantile(data$year, probs = seq(0, 1, length.out = num_knots))

# Bases of the spline (without intercept)
bases <- bs(
  data$year,
  knots = knots[-c(1, num_knots)],
  degree = 3,
  intercept = FALSE
)

# Fit a spline regression (omit the intercept)
mod <- quap(
  alist(
    doy ~ dnorm(mu, sigma),
    mu <- base %*% w,
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ),
  data = list(doy = data$doy, base = bases),
  start = list(w = rep(0, ncol(bases)))
)

# Sample posterior mean outcomes
post <- link(mod)

# Summarize posterior distributions
means <- apply(post, 2, mean)
pis <- apply(post, 2, PI, prob = 0.89)

# Assemble in a tibble
posteriors <- cbind(means, t(pis)) %>%
  as_tibble() %>%
  rename(mean_mu = "means", pi05_mu = "5%", pi94_mu = "94%")

# Predict outcomes from the model
preds <- sim(mod)

# Summarize predicted distributions
means <- apply(preds, 2, mean)
pis <- apply(preds, 2, PI, prob = 0.89)

# Assemble in a tibble
predictions <- cbind(means, t(pis)) %>%
  as_tibble() %>%
  rename(mean_pred = "means", pi05_pred = "5%", pi94_pred = "94%")

# Plot the data and the model
data %>%
  cbind(posteriors, predictions) %>%
  ggplot(aes(x = year, y = doy)) +

  # Observed data
  geom_point() +

  # Spline regression line
  geom_line(mapping = aes(y = mean_mu)) +

  # Add a ribbon for the posterior uncertainty
  geom_ribbon(
    mapping = aes(
      xmin = year, xmax = year,
      ymin = pi05_mu, ymax = pi94_mu
    ),
    alpha = 0.7
  ) +

  # Add a lighter ribbon for the predictions
  geom_ribbon(
    mapping = aes(
      xmin = year, xmax = year,
      ymin = pi05_pred, ymax = pi94_pred
    ),
    alpha = 0.3
  ) +
  xlab("Year") +
  ylab("Day of year") +
  ggtitle("Spline regression fitted without intercept")

ggsave("results/exercise_4H8.png", width = 4, height = 3, dpi = 300)

