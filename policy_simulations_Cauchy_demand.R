
###############################################################################
# (A) INSTALL / LOAD REQUIRED PACKAGES
###############################################################################
if (!require("survival")) install.packages("survival", dependencies = TRUE)
if (!require("evd"))      install.packages("evd",      dependencies = TRUE)
if (!require("dplyr"))    install.packages("dplyr",    dependencies = TRUE)

library(survival)
library(evd)
library(dplyr)

###############################################################################
# (B) DEMAND GENERATION VIA CAUCHY
#     location=100, scale in c(10,100,1000)
###############################################################################
# We'll simply generate via base R's rcauchy():
#   rcauchy(n, location=..., scale=...)
#
# Note: Cauchy can produce negative or extremely large values, but we'll
#       leave them as-is. If you prefer, you could clamp negative demands to 0.
###############################################################################
generate_cauchy_demand <- function(n = 50, 
                                   location = 100, 
                                   scale    = 10, 
                                   seed     = 123) {
  set.seed(seed)
  rcauchy(n, location = location, scale = scale)
}

###############################################################################
# (C) KAPLAN???MEIER FUNCTIONS (Similar to the original code)
###############################################################################
km_fit <- function(time, event) {
  survfit(Surv(time, event) ~ 1, type = "kaplan-meier")
}

km_to_stepfun <- function(sf) {
  timeVec <- sf$time
  survVec <- sf$surv
  cdfVec  <- 1 - survVec
  
  Sfun <- function(t) {
    if (length(timeVec) == 0) return(1)
    if (t < timeVec[1])      return(1)
    if (t >= timeVec[length(timeVec)]) {
      return(survVec[length(survVec)])
    }
    idx <- which(timeVec <= t)
    survVec[max(idx)]
  }
  
  Ffun <- function(t) 1 - Sfun(t)
  
  list(timeVec = timeVec, survVec = survVec, cdfVec = cdfVec,
       Sfun    = Sfun,    Ffun    = Ffun)
}

km_expected_cost <- function(Y, timeVec, cdfVec, b, g) {
  # b = shortage cost, g = holding cost
  if (length(timeVec) < 1) return(NA)
  
  pmf <- numeric(length(timeVec))
  pmf[1] <- cdfVec[1]
  if (length(timeVec) > 1) {
    for (i in 2:length(timeVec)) {
      pmf[i] <- cdfVec[i] - cdfVec[i-1]
    }
  }
  overage  <- pmax(Y - timeVec, 0)
  underage <- pmax(timeVec - Y, 0)
  E_over   <- sum(overage  * pmf)
  E_under  <- sum(underage * pmf)
  
  g * E_over + b * E_under
}

# We'll do a grid search from 0..(1.2 * max(timeVec))
km_optimal_cost <- function(time, event, b, g) {
  if (!any(event == 1)) {
    # no uncensored data => fallback
    return(list(Y_star = 0, cost_min = NA))
  }
  
  sf <- km_fit(time, event)
  if (length(sf$time) < 1) {
    return(list(Y_star = 0, cost_min = NA))
  }
  
  km_steps <- km_to_stepfun(sf)
  if (length(km_steps$timeVec) < 1) {
    return(list(Y_star = 0, cost_min = NA))
  }
  
  max_t  <- max(km_steps$timeVec)
  Y_grid <- seq(0, max_t * 1.2, length.out = 100)
  cost_grid <- sapply(Y_grid, function(yy) {
    km_expected_cost(yy, km_steps$timeVec, km_steps$cdfVec, b, g)
  })
  idx_min <- which.min(cost_grid)
  
  list(Y_star = Y_grid[idx_min],
       cost_min = cost_grid[idx_min])
}

###############################################################################
# (D) GPD FITTING & COST
###############################################################################
gpd_negLL <- function(params, y) {
  xi   <- params[1]
  beta <- params[2]
  if (beta <= 0) return(1e10)
  
  inside <- 1 + (xi * y) / beta
  if (any(inside <= 0)) return(1e10)
  
  log_f <- -log(beta) - (1/xi + 1)* log(inside)
  -sum(log_f)
}

gpd_fit <- function(y_data, start_params = c(xi = 0.1, beta = 1)) {
  optim(par = start_params, fn = gpd_negLL, y = y_data, method = "BFGS")
}

gpd_cdf <- function(y, xi, beta) {
  out <- numeric(length(y))
  for (i in seq_along(y)) {
    if (y[i] < 0) {
      out[i] <- 0
    } else {
      out[i] <- 1 - (1 + xi * y[i] / beta)^(-1/xi)
    }
  }
  out
}

gpd_expected_cost <- function(Y, data, threshold, b, g, xi, beta) {
  n <- length(data)
  d_leq <- data[data <= threshold]
  p_leq <- length(d_leq) / n
  
  over_leq   <- pmax(Y - d_leq, 0)
  under_leq  <- pmax(d_leq - Y, 0)
  E_over_leq  <- mean(over_leq)
  E_under_leq <- mean(under_leq)
  
  xMax <- max(data) - threshold
  if (xMax < 1e-6) {
    return(g * p_leq * E_over_leq + b * p_leq * E_under_leq)
  }
  
  nGrid   <- 200
  xGrid   <- seq(0, xMax, length.out = nGrid)
  cdfVals <- gpd_cdf(xGrid, xi, beta)
  pmfVals <- c(cdfVals[1], diff(cdfVals))
  if (sum(pmfVals) <= 1e-15) {
    return(g * p_leq * E_over_leq + b * p_leq * E_under_leq)
  }
  pmfVals <- pmfVals / sum(pmfVals)
  
  over_exceed   <- pmax(Y - (threshold + xGrid), 0)
  under_exceed  <- pmax((threshold + xGrid) - Y, 0)
  E_over_exceed  <- sum(over_exceed  * pmfVals)
  E_under_exceed <- sum(under_exceed * pmfVals)
  
  E_over  <- p_leq * E_over_leq  + (1 - p_leq)*E_over_exceed
  E_under <- p_leq * E_under_leq + (1 - p_leq)*E_under_exceed
  
  g * E_over + b * E_under
}

gpd_optimal_cost <- function(demands, b, g, thresh_prob = 0.8) {
  if (length(demands) < 2) {
    return(list(Y_star = mean(demands), cost_min = NA))
  }
  threshold <- quantile(demands, probs = thresh_prob)
  exceedances <- demands[demands > threshold] - threshold
  if (length(exceedances) < 2) {
    return(list(Y_star = mean(demands), cost_min = NA))
  }
  
  fit     <- gpd_fit(exceedances)
  xi_hat   <- fit$par[1]
  beta_hat <- fit$par[2]
  
  Y_grid <- seq(0, max(demands)*1.2, length.out = 100)
  cost_grid <- sapply(Y_grid, function(yy) {
    gpd_expected_cost(yy, demands, threshold, b, g, xi_hat, beta_hat)
  })
  idx_min <- which.min(cost_grid)
  
  list(Y_star = Y_grid[idx_min],
       cost_min = cost_grid[idx_min])
}

###############################################################################
# (E) ADDITIONAL POLICIES: (s,S), Base-Stock, EOQ
###############################################################################
# We'll fix some constants for demonstration:
s_const          <- 90
S_const          <- 120
base_stock_const <- 110
ordering_cost    <- 50

compute_sS_policy <- function(demands, s, S, b, g) {
  current_inventory <- tail(demands, 1)
  reorder_amount    <- ifelse(current_inventory < s, S - current_inventory, 0)
  
  overage  <- pmax(reorder_amount - demands, 0)
  underage <- pmax(demands - reorder_amount, 0)
  cost <- g * mean(overage) + b * mean(underage)
  
  list(Y_star = reorder_amount, cost_min = cost)
}

compute_base_stock <- function(demands, base_stock, b, g) {
  current_inventory <- tail(demands, 1)
  reorder_amount    <- max(0, base_stock - current_inventory)
  
  overage  <- pmax(reorder_amount - demands, 0)
  underage <- pmax(demands - reorder_amount, 0)
  cost <- g * mean(overage) + b * mean(underage)
  
  list(Y_star = reorder_amount, cost_min = cost)
}

compute_eoq <- function(demands, b, g, ordering_cost) {
  demand_rate    <- mean(demands)
  reorder_amount <- sqrt((2 * demand_rate * ordering_cost) / g) 
  # if g=0 then this formula is undefined, but let's leave it as-is
  # or do a safety check:
  if (g <= 0) {
    reorder_amount <- max(demands) # fallback
  }
  
  overage  <- pmax(reorder_amount - demands, 0)
  underage <- pmax(demands - reorder_amount, 0)
  cost <- g * mean(overage) + b * mean(underage)
  
  list(Y_star = reorder_amount, cost_min = cost)
}

###############################################################################
# (F) MAIN SIMULATION LOOP
# b in {0,1,2}, g in {0,1,2},
# location=100 (fixed),
# scale in {10,100,1000},
# n_sims=30, n_periods=50
###############################################################################
b_vals     <- c(1,2)
g_vals     <- c(1,2)
loc_par    <- 100
scale_vals <- c(100,1000)

n_sims    <- 30
n_periods <- 50

all_results <- data.frame()

# We form a parameter grid over b, g, scale
param_grid <- expand.grid(
  b     = b_vals,
  g     = g_vals,
  scale = scale_vals
)

for (i in seq_len(nrow(param_grid))) {
  b_par    <- param_grid$b[i]
  g_par    <- param_grid$g[i]
  scale_par<- param_grid$scale[i]
  
  # For storing the 30 x 50 results
  Ystar_KM_mat    <- matrix(NA, n_sims, n_periods)
  Cost_KM_mat     <- matrix(NA, n_sims, n_periods)
  
  Ystar_GPD_mat   <- matrix(NA, n_sims, n_periods)
  Cost_GPD_mat    <- matrix(NA, n_sims, n_periods)
  
  Ystar_sS_mat    <- matrix(NA, n_sims, n_periods)
  Cost_sS_mat     <- matrix(NA, n_sims, n_periods)
  
  Ystar_Base_mat  <- matrix(NA, n_sims, n_periods)
  Cost_Base_mat   <- matrix(NA, n_sims, n_periods)
  
  Ystar_EOQ_mat   <- matrix(NA, n_sims, n_periods)
  Cost_EOQ_mat    <- matrix(NA, n_sims, n_periods)
  
  for (sim in 1:n_sims) {
    # Generate Cauchy demands of length 50
    seed_offset <- 100000 * i + sim
    demands <- generate_cauchy_demand(
      n        = n_periods,
      location = loc_par,
      scale    = scale_par,
      seed     = seed_offset
    )
    
    # For each t=1..50, compute each policy
    for (t in 1:n_periods) {
      d_sub <- demands[1:t]
      
      # ========== 1) KM
      # We'll do the same censor threshold = 150
      time_t  <- ifelse(d_sub >= 150, 150, d_sub)
      event_t <- ifelse(d_sub >= 150, 0, 1)
      km_res  <- km_optimal_cost(time_t, event_t, b_par, g_par)
      Ystar_KM_mat[sim, t] <- km_res$Y_star
      Cost_KM_mat[sim, t]  <- km_res$cost_min
      
      # ========== 2) GPD
      gpd_res <- gpd_optimal_cost(d_sub, b_par, g_par, thresh_prob = 0.8)
      Ystar_GPD_mat[sim, t] <- gpd_res$Y_star
      Cost_GPD_mat[sim, t]  <- gpd_res$cost_min
      
      # ========== 3) (s,S)
      sS_res <- compute_sS_policy(d_sub, s_const, S_const, b_par, g_par)
      Ystar_sS_mat[sim, t] <- sS_res$Y_star
      Cost_sS_mat[sim, t]  <- sS_res$cost_min
      
      # ========== 4) Base-Stock
      base_res <- compute_base_stock(d_sub, base_stock_const, b_par, g_par)
      Ystar_Base_mat[sim, t] <- base_res$Y_star
      Cost_Base_mat[sim, t]  <- base_res$cost_min
      
      # ========== 5) EOQ
      eoq_res <- compute_eoq(d_sub, b_par, g_par, ordering_cost)
      Ystar_EOQ_mat[sim, t] <- eoq_res$Y_star
      Cost_EOQ_mat[sim, t]  <- eoq_res$cost_min
    }
  }
  
  # Average across the 30 simulations for each t=1..50
  for (t in 1:n_periods) {
    row_out <- data.frame(
      b       = b_par,
      g       = g_par,
      location= loc_par,
      scale   = scale_par,
      t       = t,
      
      Mean_Ystar_KM   = mean(Ystar_KM_mat[, t],   na.rm = TRUE),
      Mean_Cost_KM    = mean(Cost_KM_mat[, t],    na.rm = TRUE),
      
      Mean_Ystar_GPD  = mean(Ystar_GPD_mat[, t],  na.rm = TRUE),
      Mean_Cost_GPD   = mean(Cost_GPD_mat[, t],   na.rm = TRUE),
      
      Mean_Ystar_sS   = mean(Ystar_sS_mat[, t],   na.rm = TRUE),
      Mean_Cost_sS    = mean(Cost_sS_mat[, t],    na.rm = TRUE),
      
      Mean_Ystar_Base = mean(Ystar_Base_mat[, t], na.rm = TRUE),
      Mean_Cost_Base  = mean(Cost_Base_mat[, t],  na.rm = TRUE),
      
      Mean_Ystar_EOQ  = mean(Ystar_EOQ_mat[, t],  na.rm = TRUE),
      Mean_Cost_EOQ   = mean(Cost_EOQ_mat[, t],   na.rm = TRUE)
    )
    
    all_results <- rbind(all_results, row_out)
    print(all_results)
  }
}

# Finally, write to CSV
write.csv(all_results, file = "newsvendor_cauchy_results.csv", row.names = FALSE)

###############################################################################
# END
###############################################################################

