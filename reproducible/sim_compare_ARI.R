library(transport)
library(ggplot2)
library(gridExtra)

# Utils
na_lag_pairs = function (x, missing_idx, h = 1) {
  res = NULL
  if (is.matrix(x)) {
    n = nrow(x)
  } else {
    n = length(x)
  }
  counter = 1
  for (i in missing_idx) {
    if (i - h > 0) {
      res = rbind(res, c(x[i - h], x[i]))
    }
    if (i + h < n + 1) {
      res = rbind(res, c(x[i], x[i + h]))
    }
  }

  return(res)
}

model_estimate = function (x, model) {
  if (model == "ARI") {
    temp = arima(x, order = c(1, 0, 1), include.mean = F)
    return (unname(temp$coef))
  } else if (model == "Cyc") {
    temp = arima(x, order = c(1, 0, 0), include.mean = F)
    return (unname(temp$coef))
  }
  
  return (NULL)
}

# ARI: I

model_name = "ARI"
base_dir = paste0("./sim_data/", model_name, "/")

truth = as.matrix(read.csv(paste0(base_dir, model_name, "_original.csv"), header = F))
x_obs = as.matrix(read.csv(paste0(base_dir, model_name, "_miss1.csv"), header = F))
benchmarks = as.matrix(read.csv(paste0(base_dir, "benchmarks_miss1.csv")))
WI_lin = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_lin_miss1.csv"), header = F))
WI_Kalman = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_Kalman_miss1.csv"), header = F))

wass_d = matrix(0, nrow = 1000, ncol = 9)
model_coef = vector(mode = "list", length = 10)
acfs = matrix(0, nrow = 1000, ncol = 10 * 6)

for (sim in 1:1000) {
  dta = truth[,sim]
  
  lin = benchmarks[,sim]
  spl = benchmarks[,1000 + sim]
  Kalman = benchmarks[,2000 + sim]
  PT = benchmarks[,3000 + sim]
  SSA = benchmarks[,4000 + sim]
  wass = WI_lin[,sim]
  kwass = WI_lin[,1000 + sim]
  wass_Kalman = WI_Kalman[,sim]
  kwass_Kalman = WI_Kalman[,1000 + sim]
  
  dta = diff(dta, 1)
  lin = diff(lin, 1)
  spl = diff(spl, 1)
  Kalman = diff(Kalman, 1)
  PT = diff(PT, 1)
  SSA = diff(SSA, 1)
  wass = diff(wass, 1)
  kwass = diff(kwass, 1)
  wass_Kalman = diff(wass_Kalman, 1)
  kwass_Kalman = diff(kwass_Kalman, 1)
  
  wass_d[sim, 1] = wasserstein(pp(embed(dta, 3)), pp(embed(lin, 3)), p = 2)
  wass_d[sim, 2] = wasserstein(pp(embed(dta, 3)), pp(embed(spl, 3)), p = 2)
  wass_d[sim, 3] = wasserstein(pp(embed(dta, 3)), pp(embed(Kalman, 3)), p = 2)
  wass_d[sim, 4] = wasserstein(pp(embed(dta, 3)), pp(embed(PT, 3)), p = 2)
  wass_d[sim, 5] = wasserstein(pp(embed(dta, 3)), pp(embed(SSA, 3)), p = 2)
  wass_d[sim, 6] = wasserstein(pp(embed(dta, 3)), pp(embed(wass, 3)), p = 2)
  wass_d[sim, 7] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass, 3)), p = 2)
  wass_d[sim, 8] = wasserstein(pp(embed(dta, 3)), pp(embed(wass_Kalman, 3)), p = 2)
  wass_d[sim, 9] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass_Kalman, 3)), p = 2)
  
  if (sim == 1) {
    model_coef[[1]] = model_estimate(dta, model_name)
    model_coef[[2]] = model_estimate(lin, model_name)
    model_coef[[3]] = model_estimate(spl, model_name)
    model_coef[[4]] = model_estimate(Kalman, model_name)
    model_coef[[5]] = model_estimate(PT, model_name)
    model_coef[[6]] = model_estimate(SSA, model_name)
    model_coef[[7]] = model_estimate(wass, model_name)
    model_coef[[8]] = model_estimate(kwass, model_name)
    model_coef[[9]] = model_estimate(wass_Kalman, model_name)
    model_coef[[10]] = model_estimate(kwass_Kalman, model_name)
  } else {
    model_coef[[1]] = rbind(model_coef[[1]], model_estimate(dta, model_name))
    model_coef[[2]] = rbind(model_coef[[2]], model_estimate(lin, model_name))
    model_coef[[3]] = rbind(model_coef[[3]], model_estimate(spl, model_name))
    model_coef[[4]] = rbind(model_coef[[4]], model_estimate(Kalman, model_name))
    model_coef[[5]] = rbind(model_coef[[5]], model_estimate(PT, model_name))
    model_coef[[6]] = rbind(model_coef[[6]], model_estimate(SSA, model_name))
    model_coef[[7]] = rbind(model_coef[[7]], model_estimate(wass, model_name))
    model_coef[[8]] = rbind(model_coef[[8]], model_estimate(kwass, model_name))
    model_coef[[9]] = rbind(model_coef[[9]], model_estimate(wass_Kalman, model_name))
    model_coef[[10]] = rbind(model_coef[[10]], model_estimate(kwass_Kalman, model_name))
  }
  
  acfs[sim, 1:6] = c(acf(dta, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 7:12] = c(acf(lin, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 13:18] = c(acf(spl, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 19:24] = c(acf(Kalman, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 25:30] = c(acf(PT, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 31:36] = c(acf(SSA, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 37:42] = c(acf(wass, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 43:48] = c(acf(kwass, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 49:54] = c(acf(wass_Kalman, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 55:60] = c(acf(kwass_Kalman, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  
  if (sim == 1) {
    plot_range = range(dta)
    plot_range = c(floor(plot_range[1]), ceiling(plot_range[2]))
    
    df1 = data.frame(x = embed(dta, 2)[,2], y = embed(dta, 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = "x_{t-1}", y = "x_t")
    
    temp = na_lag_pairs(lin, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation",
           x = "x_{t-1}",
           y = "x_t")
    
    # temp = na_lag_pairs(spl, which(is.na(x_obs[,1])))
    # temp = data.frame(x = temp[,1],
    #                   y = temp[,2])
    # p3 = ggplot(temp, aes(x = x, y = y)) +
    #   geom_point(color = "steelblue1", size = 1) +
    #   coord_cartesian(xlim = plot_range, ylim = plot_range) +
    #   theme_minimal() +
    #   labs(title = "Spline smoothing",
    #        x = "x_{t-1}",
    #        y = "x_t")
    temp = na_lag_pairs(SSA, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) +
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "iSSA",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(Kalman, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(PT, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(wass, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(kwass, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "k-WI (linear)",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(wass_Kalman, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(kwass_Kalman, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "k-WI (Kalman)",
           x = "x_{t-1}",
           y = "x_t")
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
  }
  cat("iteration", sim, "\n")
}

par(mfrow = c(3, 3))
if (is.vector(model_coef[[1]])) {
  for (i in 1:10) {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
      next
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "iSSA"
    } else if (i == 7) {
      method = "WI (linear)"
    } else if (i == 8) {
      method = "kWI (linear)"
    } else if (i == 9){
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    
    if (i == 1) {
      hist(model_coef[[1]], prob = TRUE, col = "lightblue", 
           main = "Ground truth", xlab = "")
      lines(density(model_coef[[1]][,j]), col = "blue", lwd = 2)
    } else {
      hist(model_coef[[i]], prob = TRUE, col = "lightblue", 
           main = method,
           xlab = "")
      lines(density(model_coef[[1]]), col = "blue", lwd = 2)
    }
  }
} else {
  for (j in 1:dim(model_coef[[1]])[2]) {
    for (i in 1:10) {
      if (i == 2) {
        method = "Linear"
      } else if (i == 3) {
        method = "Spline"
        next
      } else if (i == 4) {
        method = "Kalman"
      } else if (i == 5) {
        method = "Scalar filter"
      } else if (i == 6) {
        method = "iSSA"
      } else if (i == 7) {
        method = "WI (linear)"
      } else if (i == 8) {
        method = "kWI (linear)"
      } else if (i == 9){
        method = "WI (Kalman)"
      } else {
        method = "kWI (Kalman)"
      }
      
      if (i == 1) {
        hist(model_coef[[1]][,j], prob = TRUE, col = "lightblue", 
             main = "Ground truth", xlab = "")
        lines(density(model_coef[[1]][,j]), col = "blue", lwd = 2)
      } else {
        hist(model_coef[[i]][,j], prob = TRUE, col = "lightblue", 
             main = method,
             xlab = "")
        lines(density(model_coef[[1]][,j]), col = "blue", lwd = 2)
      }
    }
  }
}
for (j in 1:3) {
  for (i in 1:10) {
    if (i == 1) {
      hist(acfs[,j + 1], prob = TRUE, col = "lightblue", 
           main = paste0("Ground truth: ACF Lag", j),
           xlab = "")
      lines(density(acfs[,j + 1]), col = "blue", lwd = 2)
    } else {
      if (i == 2) {
        method = "Linear"
      } else if (i == 3) {
        method = "Spline"
        next
      } else if (i == 4) {
        method = "Kalman"
      } else if (i == 5) {
        method = "Scalar filter"
      } else if (i == 6) {
        method = "iSSA"
      } else if (i == 7) {
        method = "WI (linear)"
      } else if (i == 8) {
        method = "kWI (linear)"
      } else if (i == 9) {
        method = "WI (Kalman)"
      } else {
        method = "kWI (Kalman)"
      }
      hist(acfs[,((i - 1) * 6) + j + 1], prob = TRUE, col = "lightblue", 
           main = paste0(method, ": ACF Lag ", j),
           xlab = "")
      lines(density(acfs[,j + 1]), col = "blue", lwd = 2)
    }
  }
}

for (j in 1:9) {
  if (j == 1) {
    cat("Wasserstein loss for", "Linear:", round(mean(wass_d[,1]), 2), "\n")
  } else if (j == 2) {
    cat("Wasserstein loss for", "Spline:", round(mean(wass_d[,2]), 2), "\n")
  } else if (j == 3) {
    cat("Wasserstein loss for", "Kalman:", round(mean(wass_d[,3]), 2), "\n")
  } else if (j == 4) {
    cat("Wasserstein loss for", "PT:", round(mean(wass_d[,4]), 2), "\n")
  } else if (j == 5) {
    cat("Wasserstein loss for", "iSSA:", round(mean(wass_d[,5]), 2), "\n")
  } else if (j == 6) {
    cat("Wasserstein loss for", "WI (linear):", round(mean(wass_d[,6]), 2), "\n")
  } else if (j == 7) {
    cat("Wasserstein loss for", "kWI (linear):", round(mean(wass_d[,7]), 2), "\n")
  } else if (j == 8) {
    cat("Wasserstein loss for", "WI (Kalman):", round(mean(wass_d[,8]), 2), "\n")
  } else {
    cat("Wasserstein loss for", "kWI (Kalman):", round(mean(wass_d[,9]), 2), "\n")
  }
}
print(matrix(round(colMeans(acfs), 2), ncol = 10))
GT_acf = colMeans(acfs[,1:6])
for (j in 1:9) {
  sub_acfs = acfs[,(j * 6 + 1):((j + 1) * 6)]
  temp = t(t(sub_acfs) - GT_acf)
  if (j == 1) {
    cat("ACF loss for", "Linear:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 2) {
    cat("ACF loss for", "Spline:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 3) {
    cat("ACF loss for", "Kalman:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 4) {
    cat("ACF loss for", "PT:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 5) {
    cat("ACF loss for", "iSSA:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 6) {
    cat("ACF loss for", "WI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 7) {
    cat("ACF loss for", "kWI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 8) {
    cat("ACF loss for", "WI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else {
    cat("ACF loss for", "kWI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  }
}

GT_arcoef = colMeans(model_coef[[1]])
for (j in 1:9) {
  sub_arcoef = model_coef[[j+1]]
  temp = t(t(sub_arcoef) - GT_arcoef)
  if (j == 1) {
    cat("estimation error for", "Linear:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 2) {
    cat("estimation error for", "Spline:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 3) {
    cat("estimation error for", "Kalman:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 4) {
    cat("estimation error for", "PT:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 5) {
    cat("estimation error for", "iSSA:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 6) {
    cat("estimation error for", "WI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 7) {
    cat("estimation error for", "kWI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 8) {
    cat("estimation error for", "WI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else {
    cat("estimation error for", "kWI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  }
}

################################################################################
# ARI: II

model_name = "ARI"
base_dir = paste0("./sim_data/", model_name, "/")

truth = as.matrix(read.csv(paste0(base_dir, model_name, "_original.csv"), header = F))
x_obs = as.matrix(read.csv(paste0(base_dir, model_name, "_miss2.csv"), header = F))
benchmarks = as.matrix(read.csv(paste0(base_dir, "benchmarks_miss2.csv")))
WI_lin = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_lin_miss2.csv"), header = F))
WI_Kalman = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_Kalman_miss2.csv"), header = F))

wass_d = matrix(0, nrow = 1000, ncol = 9)
model_coef = vector(mode = "list", length = 10)
acfs = matrix(0, nrow = 1000, ncol = 10 * 6)

for (sim in 1:1000) {
  dta = truth[,sim]
  
  lin = benchmarks[,sim]
  spl = benchmarks[,1000 + sim]
  Kalman = benchmarks[,2000 + sim]
  PT = benchmarks[,3000 + sim]
  SSA = benchmarks[,4000 + sim]
  wass = WI_lin[,sim]
  kwass = WI_lin[,1000 + sim]
  wass_Kalman = WI_Kalman[,sim]
  kwass_Kalman = WI_Kalman[,1000 + sim]
  
  dta = diff(dta, 1)
  lin = diff(lin, 1)
  spl = diff(spl, 1)
  Kalman = diff(Kalman, 1)
  PT = diff(PT, 1)
  SSA = diff(SSA, 1)
  wass = diff(wass, 1)
  kwass = diff(kwass, 1)
  wass_Kalman = diff(wass_Kalman, 1)
  kwass_Kalman = diff(kwass_Kalman, 1)
  
  wass_d[sim, 1] = wasserstein(pp(embed(dta, 3)), pp(embed(lin, 3)), p = 2)
  wass_d[sim, 2] = wasserstein(pp(embed(dta, 3)), pp(embed(spl, 3)), p = 2)
  wass_d[sim, 3] = wasserstein(pp(embed(dta, 3)), pp(embed(Kalman, 3)), p = 2)
  wass_d[sim, 4] = wasserstein(pp(embed(dta, 3)), pp(embed(PT, 3)), p = 2)
  wass_d[sim, 5] = wasserstein(pp(embed(dta, 3)), pp(embed(SSA, 3)), p = 2)
  wass_d[sim, 6] = wasserstein(pp(embed(dta, 3)), pp(embed(wass, 3)), p = 2)
  wass_d[sim, 7] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass, 3)), p = 2)
  wass_d[sim, 8] = wasserstein(pp(embed(dta, 3)), pp(embed(wass_Kalman, 3)), p = 2)
  wass_d[sim, 9] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass_Kalman, 3)), p = 2)
  
  if (sim == 1) {
    model_coef[[1]] = model_estimate(dta, model_name)
    model_coef[[2]] = model_estimate(lin, model_name)
    model_coef[[3]] = model_estimate(spl, model_name)
    model_coef[[4]] = model_estimate(Kalman, model_name)
    model_coef[[5]] = model_estimate(PT, model_name)
    model_coef[[6]] = model_estimate(SSA, model_name)
    model_coef[[7]] = model_estimate(wass, model_name)
    model_coef[[8]] = model_estimate(kwass, model_name)
    model_coef[[9]] = model_estimate(wass_Kalman, model_name)
    model_coef[[10]] = model_estimate(kwass_Kalman, model_name)
  } else {
    model_coef[[1]] = rbind(model_coef[[1]], model_estimate(dta, model_name))
    model_coef[[2]] = rbind(model_coef[[2]], model_estimate(lin, model_name))
    model_coef[[3]] = rbind(model_coef[[3]], model_estimate(spl, model_name))
    model_coef[[4]] = rbind(model_coef[[4]], model_estimate(Kalman, model_name))
    model_coef[[5]] = rbind(model_coef[[5]], model_estimate(PT, model_name))
    model_coef[[6]] = rbind(model_coef[[6]], model_estimate(SSA, model_name))
    model_coef[[7]] = rbind(model_coef[[7]], model_estimate(wass, model_name))
    model_coef[[8]] = rbind(model_coef[[8]], model_estimate(kwass, model_name))
    model_coef[[9]] = rbind(model_coef[[9]], model_estimate(wass_Kalman, model_name))
    model_coef[[10]] = rbind(model_coef[[10]], model_estimate(kwass_Kalman, model_name))
  }
  
  acfs[sim, 1:6] = c(acf(dta, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 7:12] = c(acf(lin, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 13:18] = c(acf(spl, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 19:24] = c(acf(Kalman, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 25:30] = c(acf(PT, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 31:36] = c(acf(SSA, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 37:42] = c(acf(wass, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 43:48] = c(acf(kwass, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 49:54] = c(acf(wass_Kalman, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  acfs[sim, 55:60] = c(acf(kwass_Kalman, lag.max = 5, type = "covariance", plot = FALSE)$acf)
  
  if (sim == 1) {
    plot_range = range(dta)
    plot_range = c(floor(plot_range[1]), ceiling(plot_range[2]))
    
    df1 = data.frame(x = embed(dta, 2)[,2], y = embed(dta, 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = "x_{t-1}", y = "x_t")
    
    temp = na_lag_pairs(lin, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation",
           x = "x_{t-1}",
           y = "x_t")
    
    # temp = na_lag_pairs(spl, which(is.na(x_obs[,1])))
    # temp = data.frame(x = temp[,1],
    #                   y = temp[,2])
    # p3 = ggplot(temp, aes(x = x, y = y)) +
    #   geom_point(color = "steelblue1", size = 1) +
    #   coord_cartesian(xlim = plot_range, ylim = plot_range) +
    #   theme_minimal() +
    #   labs(title = "Spline smoothing",
    #        x = "x_{t-1}",
    #        y = "x_t")
    temp = na_lag_pairs(SSA, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) +
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "iSSA",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(Kalman, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(PT, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(wass, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(kwass, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "k-WI (linear)",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(wass_Kalman, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)",
           x = "x_{t-1}",
           y = "x_t")
    
    temp = na_lag_pairs(kwass_Kalman, which(is.na(x_obs[,1])))
    temp = data.frame(x = temp[,1],
                      y = temp[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "k-WI (Kalman)",
           x = "x_{t-1}",
           y = "x_t")
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
  }
  cat("iteration", sim, "\n")
}

par(mfrow = c(3, 3))
if (is.vector(model_coef[[1]])) {
  for (i in 1:10) {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
      next
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "iSSA"
    } else if (i == 7) {
      method = "WI (linear)"
    } else if (i == 8) {
      method = "kWI (linear)"
    } else if (i == 9){
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    
    if (i == 1) {
      hist(model_coef[[1]], prob = TRUE, col = "lightblue", 
           main = "Ground truth", xlab = "")
      lines(density(model_coef[[1]][,j]), col = "blue", lwd = 2)
    } else {
      hist(model_coef[[i]], prob = TRUE, col = "lightblue", 
           main = method,
           xlab = "")
      lines(density(model_coef[[1]]), col = "blue", lwd = 2)
    }
  }
} else {
  for (j in 1:dim(model_coef[[1]])[2]) {
    for (i in 1:10) {
      if (i == 2) {
        method = "Linear"
      } else if (i == 3) {
        method = "Spline"
        next
      } else if (i == 4) {
        method = "Kalman"
      } else if (i == 5) {
        method = "Scalar filter"
      } else if (i == 6) {
        method = "iSSA"
      } else if (i == 7) {
        method = "WI (linear)"
      } else if (i == 8) {
        method = "kWI (linear)"
      } else if (i == 9){
        method = "WI (Kalman)"
      } else {
        method = "kWI (Kalman)"
      }
      
      if (i == 1) {
        hist(model_coef[[1]][,j], prob = TRUE, col = "lightblue", 
             main = "Ground truth", xlab = "")
        lines(density(model_coef[[1]][,j]), col = "blue", lwd = 2)
      } else {
        hist(model_coef[[i]][,j], prob = TRUE, col = "lightblue", 
             main = method,
             xlab = "")
        lines(density(model_coef[[1]][,j]), col = "blue", lwd = 2)
      }
    }
  }
}
for (j in 1:3) {
  for (i in 1:10) {
    if (i == 1) {
      hist(acfs[,j + 1], prob = TRUE, col = "lightblue", 
           main = paste0("Ground truth: ACF Lag", j),
           xlab = "")
      lines(density(acfs[,j + 1]), col = "blue", lwd = 2)
    } else {
      if (i == 2) {
        method = "Linear"
      } else if (i == 3) {
        method = "Spline"
        next
      } else if (i == 4) {
        method = "Kalman"
      } else if (i == 5) {
        method = "Scalar filter"
      } else if (i == 6) {
        method = "iSSA"
      } else if (i == 7) {
        method = "WI (linear)"
      } else if (i == 8) {
        method = "kWI (linear)"
      } else if (i == 9) {
        method = "WI (Kalman)"
      } else {
        method = "kWI (Kalman)"
      }
      hist(acfs[,((i - 1) * 6) + j + 1], prob = TRUE, col = "lightblue", 
           main = paste0(method, ": ACF Lag ", j),
           xlab = "")
      lines(density(acfs[,j + 1]), col = "blue", lwd = 2)
    }
  }
}

for (j in 1:9) {
  if (j == 1) {
    cat("Wasserstein loss for", "Linear:", round(mean(wass_d[,1]), 2), "\n")
  } else if (j == 2) {
    cat("Wasserstein loss for", "Spline:", round(mean(wass_d[,2]), 2), "\n")
  } else if (j == 3) {
    cat("Wasserstein loss for", "Kalman:", round(mean(wass_d[,3]), 2), "\n")
  } else if (j == 4) {
    cat("Wasserstein loss for", "PT:", round(mean(wass_d[,4]), 2), "\n")
  } else if (j == 5) {
    cat("Wasserstein loss for", "iSSA:", round(mean(wass_d[,5]), 2), "\n")
  } else if (j == 6) {
    cat("Wasserstein loss for", "WI (linear):", round(mean(wass_d[,6]), 2), "\n")
  } else if (j == 7) {
    cat("Wasserstein loss for", "kWI (linear):", round(mean(wass_d[,7]), 2), "\n")
  } else if (j == 8) {
    cat("Wasserstein loss for", "WI (Kalman):", round(mean(wass_d[,8]), 2), "\n")
  } else {
    cat("Wasserstein loss for", "kWI (Kalman):", round(mean(wass_d[,9]), 2), "\n")
  }
}
print(matrix(round(colMeans(acfs), 2), ncol = 10))
GT_acf = colMeans(acfs[,1:6])
for (j in 1:9) {
  sub_acfs = acfs[,(j * 6 + 1):((j + 1) * 6)]
  temp = t(t(sub_acfs) - GT_acf)
  if (j == 1) {
    cat("ACF loss for", "Linear:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 2) {
    cat("ACF loss for", "Spline:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 3) {
    cat("ACF loss for", "Kalman:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 4) {
    cat("ACF loss for", "PT:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 5) {
    cat("ACF loss for", "iSSA:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 6) {
    cat("ACF loss for", "WI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 7) {
    cat("ACF loss for", "kWI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 8) {
    cat("ACF loss for", "WI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else {
    cat("ACF loss for", "kWI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  }
}

GT_arcoef = colMeans(model_coef[[1]])
for (j in 1:9) {
  sub_arcoef = model_coef[[j+1]]
  temp = t(t(sub_arcoef) - GT_arcoef)
  if (j == 1) {
    cat("estimation error for", "Linear:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 2) {
    cat("estimation error for", "Spline:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 3) {
    cat("estimation error for", "Kalman:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 4) {
    cat("estimation error for", "PT:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 5) {
    cat("estimation error for", "iSSA:", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 6) {
    cat("estimation error for", "WI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 7) {
    cat("estimation error for", "kWI (linear):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else if (j == 8) {
    cat("estimation error for", "WI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  } else {
    cat("estimation error for", "kWI (Kalman):", round(sqrt(colMeans(temp^2)), 2), "\n")
  }
}
