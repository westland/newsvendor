# Install and load necessary packages
if (!require("mnonr")) install.packages("mnonr", dependencies = TRUE)
if (!require("moments")) install.packages("moments", dependencies = TRUE)
if (!require("survival")) install.packages("survival", dependencies = TRUE)
library(mnonr)
library(moments)
library(survival)
library(ggplot2)

# Data Generation Function
generate_data <- function(n, mu, variance, skewness, kurtosis) {
  sd <- sqrt(variance)
  Sigma <- matrix(c(variance, 0, 0, 1), nrow = 2, ncol = 2)
  generated_data <- unonr(
    n = n,
    mu = mu,
    Sigma = Sigma,
    skewness = skewness,
    kurtosis = kurtosis
  )
  return(generated_data[, 1])
}

# Kaplan-Meier Survival Function
km_survival <- function(x, x_vals, d_vals, n_vals, u) {
  numerator <- prod(1 - d_vals[x_vals <= x] / n_vals[x_vals <= x])
  denominator <- prod(1 - d_vals[x_vals <= u] / n_vals[x_vals <= u])
  return(numerator / denominator)
}

# GPD Survival Function
gpd_survival <- function(x, u, xi, beta) {
  if (x < u) {
    return(NA)
  }
  term <- 1 + xi * (x - u) / beta
  if (term <= 0) {
    stop("Invalid parameters: GPD survival function undefined.")
  }
  return(term^(-1 / xi))
}
