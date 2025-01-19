# Use the built-in 'optim' function to minimize the negative 
# log-likelihood (NLL). Note that real applications often rely on 
# specialized packages (ismev, extRemes, evd, etc.) for threshold 
# selection, stability checks, and robust fitting.

#############################
# 1. Generate Synthetic Data
#############################
set.seed(123)

# We'll simulate exceedances from a true GPD with known shape (xi_true) 
# and scale (beta_true) for demonstration. Typically, you'd extract 
# exceedances above a chosen threshold u from your dataset.

# For simulation, we can use the 'evd' package's rgpd() if available. 
# Alternatively, we simulate a uniform random variable and invert 
# the GPD CDF. Here, we demonstrate a direct method if 'evd' is installed.

# install.packages("evd")  # if not already installed
library(evd)

xi_true   <- 0.2   # True shape
beta_true <- 1.0   # True scale
n         <- 500   # Number of exceedances to simulate

# Simulate y ~ GPD(xi=xi_true, mu=0, beta=beta_true)
y_data <- rgpd(n, loc = 0, scale = beta_true, shape = xi_true)

# In real life, y_data would be the exceedances above threshold u:
#    y_data = D_i - u, for all D_i > u.
# For demonstration, we directly sampled from the GPD distribution.

###############################################
# 2. Define the Negative Log-Likelihood (NLL)
###############################################
# The GPD PDF (for y > 0, xi != 0) is:
#   f(y | xi, beta) = 1 / beta * (1 + xi*y / beta)^(-1/xi - 1),
# subject to:
#   beta > 0
#   1 + xi*y / beta > 0 for all data points.

gpd_negLL <- function(params, y) {
  xi   <- params[1]
  beta <- params[2]
  
  # Impose constraints: beta must be positive, and 
  # 1 + (xi*y)/beta must stay > 0 for all y.
  # If constraints are violated, return a large penalty.
  if (beta <= 0) return(1e10)
  
  # Check positivity inside (1 + xi*y/beta)
  inside <- 1 + (xi * y) / beta
  if (any(inside <= 0)) return(1e10)
  
  # Negative log-likelihood for GPD:
  #   - sum(log( f(y_i) ))
  #   f(y_i) = (1/beta) * [1 + (xi*y_i)/beta]^(-(1/xi) - 1)
  n <- length(y)
  
  # Log PDF of GPD:
  #   log f(y_i) = - log(beta) - (1/xi + 1)*log(1 + xi*y_i/beta)
  log_f <- -log(beta) - (1/xi + 1) * log(inside)
  
  # Return NEGATIVE log-likelihood:
  return(-sum(log_f))
}

###################################################
# 3. Use an Optimizer (optim) to Fit (xi, beta)
###################################################
# We'll provide an initial guess for (xi, beta). If you have a prior, 
# you might pick values more carefully. Using "method='BFGS'" or 
# "method='L-BFGS-B'" is typical. "L-BFGS-B" allows box constraints, 
# but here we just implement a penalty for invalid parameters.

start_params <- c(xi = 0.1, beta = 1.0)

# Minimize the negative log-likelihood:
fit <- optim(
  par    = start_params,
  fn     = gpd_negLL,
  y      = y_data,      # pass data via '...' 
  method = "BFGS"       # or "L-BFGS-B", "Nelder-Mead", etc.
)

# Check results:
fit$convergence        # 0 => converged
est_xi   <- fit$par[1]
est_beta <- fit$par[2]

cat("MLE via optim:\n")
cat("  xi   =", est_xi, "\n")
cat("  beta =", est_beta, "\n")

# Compare to the true values:
cat("\nTrue xi   =", xi_true, "\n")
cat("True beta =", beta_true, "\n")

# NOTE:
# 1. In practice, ensure thorough checks on convergence status, 
#    parameter boundaries, standard errors, etc.
# 2. Specialized extreme value packages handle threshold selection 
#    and advanced diagnostics automatically (ismev, evd, extRemes).


