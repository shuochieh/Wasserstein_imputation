# Analysis using groundwater data

AR_loglik = function (x, max_order = 10) {
  logliks = rep(NA, max_order)
  bics = rep(NA, max_order)
  n = length(x)
  
  for (p in 1:max_order) {
    fit = ar(x, aic = FALSE, order.max = p, na.action = na.pass)
    bics[p] = n * log(mean(fit$resid^2, na.rm = T)) + log(n) * p
  }

  
  return (list("logliks" = NULL, "BIC" = bics))
}
AR_rolling = function (x, x_truth, order, training_n) {
  n = length(x)
  if (training_n >= (n - order)) {
    stop("training sample size too large")
  }
  
  model = ar(x[1:training_n], aic = FALSE, order.max = order, demean = FALSE)
  
  counter = 0
  res = 0
  for (i in ((training_n + 1):n)) {
    x_aux = x_truth[(i - 1):(i - order)]
    if (sum(is.na(x_aux)) > 0) {
      next
    }
    x_hat = c(x_aux %*% model$ar)
    
    if (is.na(x_truth[i])) {
      next
    } else {
      counter = counter + 1
      res = res + unname((x_truth[i] - x_hat)^2)
    }
  }
  
  # cat("N. obs =", counter, "\n")
  
  return (sqrt(res / counter))
}

ground_truth = as.matrix(read.csv("./real_data/small_series2.csv", header = FALSE))
ground_truth = t(ground_truth)
d = ncol(ground_truth)

for (pct in c(0.4)) {
  
  x_obs = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_masked.csv"), header = FALSE))
  x_obs = t(x_obs)
  
  LP_trend = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_LPfit.csv"), 
                                header = FALSE))
  LP_trend = t(LP_trend)
  
  par(mfrow = c(5, 6), xpd = FALSE)
  for (j in 1:15) {
    plot(ground_truth[,j], type = "l", col = "darkblue",
         xlab = "", ylab = "", bty = "L", xaxt = "n",
         main = paste("Series", j))
    plot(x_obs[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L", xaxt = "n",
         main = paste("Series", j, "observed"))
    lines(LP_trend[,j], col = "lightblue")
  }

  par(mfrow = c(3, 5), xpd = FALSE)
  for (j in 1:15) {
    plot(x_obs[,j] - LP_trend[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L", xaxt = "n",
         main = paste("Series", j, "residuals"))
  }
  
}

dev.off()
all_pred_err = array(NA, dim = c(6, 7, d))
for (pct in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  
  max_order = 15
  ar_loglik = array(NA, dim = c(8, d, max_order))
  ar_ic = array(NA, dim = c(8, d, max_order))
  selected_order = array(NA, dim = c(8, d))
  acf_hat = array(NA, dim = c(8, d, 3))
  pred_err = array(NA, dim = c(7, d))
  
  x_obs = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_masked.csv"), header = FALSE))
  x_obs = t(x_obs)
  
  linimp = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_lin.csv"), 
                              header = FALSE))
  linimp = t(linimp)
  
  SSA = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SSA.csv"), 
                           header = FALSE))
  SSA = t(SSA)
  
  KS = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_KS.csv"), 
                          header = FALSE))
  KS = t(KS)
  
  SF = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SF.csv"), 
                          header = FALSE))
  SF = t(SF)
  
  mean_interpolation = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_mean.csv"), 
                                          header = FALSE))
  mean_interpolation = t(mean_interpolation)
  
  TWI = as.matrix(read.csv(paste0("./real_data/WI_small", 10 * pct, ".csv"), 
                           header = FALSE))
  kTWI = as.matrix(read.csv(paste0("./real_data/kWI_small", 10 * pct, ".csv"), 
                            header = FALSE))
  LP_trend = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_LPfit.csv"), 
                                header = FALSE))
  LP_trend = t(LP_trend)
  
  TWI = TWI + LP_trend
  kTWI = kTWI + LP_trend

  for (j in 1:d) {
    # cat("Series", j, "\n")
    
    model = AR_loglik(ground_truth[,j], max_order)
    ar_ic[1,j,] = model$BIC
    acf_hat[1,j,] = acf(ground_truth[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]

    model = AR_loglik(linimp[,j], max_order)
    ar_ic[2,j,] = model$BIC
    acf_hat[2,j,] = acf(linimp[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[1,j] = AR_rolling(linimp[,j], ground_truth[,j], max_order, 200)
    
    model = AR_loglik(SSA[,j], max_order)
    ar_ic[3,j,] = model$BIC
    acf_hat[3,j,] = acf(SSA[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[2,j] = AR_rolling(SSA[,j], ground_truth[,j], max_order, 200)
    
    model = AR_loglik(KS[,j], max_order)
    ar_ic[4,j,] = model$BIC
    acf_hat[4,j,] = acf(KS[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[3,j] = AR_rolling(KS[,j], ground_truth[,j], max_order, 200)
    
    model = AR_loglik(SF[,j], max_order)
    ar_ic[5,j,] = model$BIC
    acf_hat[5,j,] = acf(SF[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[4,j] = AR_rolling(SF[,j], ground_truth[,j], max_order, 200)
    
    model = AR_loglik(mean_interpolation[,j], max_order)
    ar_ic[6,j,] = model$BIC
    acf_hat[6,j,] = acf(mean_interpolation[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[5,j] = AR_rolling(mean_interpolation[,j], ground_truth[,j], max_order, 200)
    
    model = AR_loglik(TWI[,j], max_order)
    ar_ic[7,j,] = model$BIC
    acf_hat[7,j,] = acf(TWI[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[6,j] = AR_rolling(TWI[,j], ground_truth[,j], max_order, 200)
    
    model = AR_loglik(kTWI[,j], max_order)
    ar_ic[8,j,] = model$BIC
    acf_hat[8,j,] = acf(kTWI[,j], lag.max = 3, na.action = na.pass, plot = F,
                        type = "covariance")$acf[2:4]
    pred_err[7,j] = AR_rolling(kTWI[,j], ground_truth[,j], max_order, 200)

    for (k in 1:8) {
      selected_order[k,j] = which.min(ar_ic[k,j,])
    }
  }

  cat("Linear:", length(which(abs(selected_order[2,] - selected_order[1,]) < 2)), "\n")
  cat("iSSA:", length(which(abs(selected_order[3,] - selected_order[1,]) < 2)), "\n")
  cat("Kalman smoothing:", length(which(abs(selected_order[4,] - selected_order[1,]) < 2)), "\n")
  cat("Scalar filter:", length(which(abs(selected_order[5,] - selected_order[1,]) < 2)), "\n")
  cat("Mean interpolation:", length(which(abs(selected_order[6,] - selected_order[1,]) < 2)), "\n")
  cat("TWI:", length(which(abs(selected_order[7,] - selected_order[1,]) < 2)), "\n")
  cat("k-TWI:", length(which(abs(selected_order[8,] - selected_order[1,]) < 2)), "\n\n")

  all_pred_err[pct * 10,,] = pred_err
}
par(mfrow = c(3, 5))
for (j in 1:d) {
  for (i in 1:7) {
    if (i == 1) {
      plot(x = 0.1 * c(1:6),
           y = all_pred_err[,i,j], type = "b", pch = 19, col = 2 + i, lty = 2,
           ylab = "", xlab = "missing ratio",
           ylim = range(all_pred_err[,,j]),
           main = paste("series", j))
      if (j == 1) {
        legend("topleft",
               legend = c("Linear", "iSSA", "Kalman", "SF", "Mean", "TWI", "k-TWI"),
               lty = c(2, 2, 2, 2, 2, 1, 1),
               col = c(3, 4, 5, 6, 7, 1, 2), 
               lwd = c(1, 1, 1, 1, 1, 2, 2),
               ncol = 2, cex = 1, bty = "n")
      }
    } else if (i == 6) {
      lines(x = 0.1 * c(1:6),
            y = all_pred_err[,i,j], type = "b", pch = 19, col = 1, lwd = 2)
    } else if (i == 7) {
      lines(x = 0.1 * c(1:6),
            y = all_pred_err[,i,j], type = "b", pch = 19, col = 2, lwd = 2)
    } else {
      lines(x = 0.1 * c(1:6),
            y = all_pred_err[,i,j], type = "b", pch = 19, col = 2 + i, lty = 2)
    }
  }
}

r_sq = function (x, xhat) {
  return (1 - sum((x - xhat)^2) / sum((x - mean(x))^2))
}
mae = function (x, xhat, type = "mean") {
  if (type == "mean") {
    return (mean(abs(x - xhat)))
  } else if (type == "median") {
    return (median(abs(x - xhat)))
  }
}
acf_type = "covariance"
for (pct in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  
  cat("\n\npercentage:", pct, "\n")
  max_order = 3    
  if (acf_type == "partial") {
    acf_hat = array(NA, dim = c(8, d, max_order))
  } else {
    acf_hat = array(NA, dim = c(8, d, max_order))
  }
  
  x_obs = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_masked.csv"), header = FALSE))
  x_obs = t(x_obs)
  
  linimp = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_lin.csv"), 
                              header = FALSE))
  linimp = t(linimp)
  
  SSA = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SSA.csv"), 
                           header = FALSE))
  SSA = t(SSA)
  
  KS = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_KS.csv"), 
                          header = FALSE))
  KS = t(KS)
  
  SF = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SF.csv"), 
                          header = FALSE))
  SF = t(SF)
  
  mean_interpolation = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_mean.csv"), 
                                          header = FALSE))
  mean_interpolation = t(mean_interpolation)
  
  TWI = as.matrix(read.csv(paste0("./real_data/WI_small", 10 * pct, ".csv"), 
                           header = FALSE))
  kTWI = as.matrix(read.csv(paste0("./real_data/kWI_small", 10 * pct, ".csv"), 
                            header = FALSE))
  LP_trend = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_LPfit.csv"), 
                                header = FALSE))
  LP_trend = t(LP_trend)
  
  TWI = TWI + LP_trend
  kTWI = kTWI + LP_trend
  

  for (j in 1:d) {
    # cat("Series", j, "\n")
    
    model = AR_loglik(ground_truth[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[1,j,] = acf(ground_truth[,j], lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[1,j,] = acf(ground_truth[,j], lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(linimp[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[2,j,] = acf(linimp[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[2,j,] = acf(linimp[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(SSA[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[3,j,] = acf(SSA[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[3,j,] = acf(SSA[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(KS[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[4,j,] = acf(KS[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[4,j,] = acf(KS[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(SF[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[5,j,] = acf(SF[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[5,j,] = acf(SF[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(mean_interpolation[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[6,j,] = acf(mean_interpolation[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[6,j,] = acf(mean_interpolation[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(TWI[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[7,j,] = acf(TWI[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[7,j,] = acf(TWI[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

    model = AR_loglik(kTWI[,j], max_order)
    if (acf_type != "partial") {
      acf_hat[8,j,] = acf(kTWI[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[8,j,] = acf(kTWI[,j], lag.max = max_order, na.action = na.fail, plot = F,
                          type = acf_type)$acf
    }

  }
  
  # ACF analysis
  methods = c("Ground truth", "Linear", "iSSA", "Kalman smoothing", "Scalar filter",
              "Mean interpolation", "TWI", "k-TWI")
  for (j in 2:8) {
    if (j < 8) {
      if (j == 2) {
        cat("Mean: ")
      }
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,])) * 10, 2), "& ")
    } else {
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,])) * 10, 2), "\n")
    }
  }
  
  for (j in 2:8) {
    if (j < 8) {
      if (j == 2) {
        cat("Median: ")
      }
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,]), "median") * 10, 2), "& ")
    } else {
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,]), "median") * 10, 2), "\n")
    }
  }
  
  # for (j in 2:8) {
  #   for (k in 1:d) {
  #     if (k == 1) {
  #       # plot(x = acf_hat[1,k,], y = acf_hat[j,k,], xlab = "True ACF", ylab = "Imputed ACF",
  #       #      main = paste0(methods[j], " (MAE: ", round(mae(c(acf_hat[1,,]), c(acf_hat[j,,])) * 10, 3), ")")
  #       #      , xlim = range(acf_hat[1,,]), ylim = range(acf_hat[1,,]),
  #       #      pch = 19, col = "steelblue", cex = 1.5)
  #       # grid(col = "gray50")
  #     } else {
  #       # points(x = acf_hat[1,k,], y = acf_hat[j,k,], pch = 19, col = "steelblue", cex = 1.5)
  #     }
  #   }
  #   # abline(a = 0, b = 1, col = "lightblue", lwd = 2, lty = 2)
  # }
  
  for (h in 1:max_order) {
    cat("lag =", h, "\n")
    for (k in 2:8) {
      cat(round(mae(acf_hat[1,,h], acf_hat[k,,h]) * 10, 2), "& ")
      # cat("...", methods[k], round(mae(acf_hat[1,,h], acf_hat[k,,h]) * 10, 2), "& ")
      if (k == 8) {
        cat("\n")
      }
    }
  }

}


dev.off()
for (pct in c(0.2, 0.4, 0.6)) {
  
  x_obs = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_masked.csv"), header = FALSE))
  x_obs = t(x_obs)
  
  TWI = as.matrix(read.csv(paste0("./real_data/WI_small", 10 * pct, ".csv"), 
                           header = FALSE))
  kTWI = as.matrix(read.csv(paste0("./real_data/kWI_small", 10 * pct, ".csv"), 
                            header = FALSE))
  LP_trend = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_LPfit.csv"), 
                                header = FALSE))
  LP_trend = t(LP_trend)
  
  TWI = TWI + LP_trend
  kTWI = kTWI + LP_trend
  
  par(mfrow = c(5, 9), mar = c(0.5, 0.5, 1.5, 0.5), xpd = FALSE)
  for (j in 1:d) {
    ylim = range(ground_truth[,j], na.rm = T)
    plot(ground_truth[,j], type = "l", col = "darkblue",
         xlab = "", ylab = "", bty = "L", xaxt = "n",
         main = paste("Series", j), lwd = 1.5,
         ylim = ylim, yaxt = "n")
    plot(TWI[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L", xaxt = "n",
         main = paste("Series", j, "(TWI)"), lwd = 1.5,
         ylim = ylim,
         yaxt = "n")
    plot(kTWI[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L", xaxt = "n",
         main = paste("Series", j, "(k-TWI)"), lwd = 1.5,
         ylim = ylim,
         yaxt = "n")
    # lines(TWI[,j], col = "lightblue")
  }
}
for (pct in c(0.4)) {
  
  max_order = 3
  acf_hat = array(NA, dim = c(8, 15, max_order + 1))
  
  x_obs = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_masked.csv"), header = FALSE))
  x_obs = t(x_obs)
  
  linimp = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_lin.csv"), 
                              header = FALSE))
  linimp = t(linimp)
  
  SSA = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SSA.csv"), 
                           header = FALSE))
  SSA = t(SSA)
  
  KS = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_KS.csv"), 
                          header = FALSE))
  KS = t(KS)
  
  SF = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SF.csv"), 
                          header = FALSE))
  SF = t(SF)
  
  mean_interpolation = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_mean.csv"), 
                                          header = FALSE))
  mean_interpolation = t(mean_interpolation)
  
  TWI = as.matrix(read.csv(paste0("./real_data/WI_small", 10 * pct, ".csv"), 
                           header = FALSE))
  kTWI = as.matrix(read.csv(paste0("./real_data/kWI_small", 10 * pct, ".csv"), 
                            header = FALSE))
  LP_trend = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_LPfit.csv"), 
                                header = FALSE))
  LP_trend = t(LP_trend)
  
  TWI = TWI + LP_trend
  kTWI = kTWI + LP_trend
  
  
  # kTWI = as.matrix(read.csv(paste0("./real_data/kWI_small", 10 * pct, ".csv"), 
  #                           header = FALSE))
  
  par(mfrow = c(3, 7), xpd = FALSE)
  for (j in 1:15) {
    plot(x_obs[,j], type = "l", col = "darkblue",
         xlab = "", ylab = "", bty = "L",
         main = "Observed data")
    grid(col = "gray50", lty = "dotted")
    plot(ground_truth[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L",
         main = "Truth")
    grid(col = "gray50", lty = "dotted")
    plot(linimp[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L",
         main = "Linear")
    grid(col = "gray50", lty = "dotted")
    plot(KS[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L",
         main = "Kalman smoothing")
    grid(col = "gray50", lty = "dotted")
    plot(SSA[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L",
         main = "iSSA")
    grid(col = "gray50", lty = "dotted")
    # plot(mean_interpolation[,j], type = "l", col = "steelblue",
    #      xlab = "", ylab = "", bty = "L",
    #      main = "mean")
    # grid(col = "gray50", lty = "dotted")
    plot(TWI[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L",
         main = "TWI")
    grid(col = "gray50", lty = "dotted")
    plot(kTWI[,j], type = "l", col = "steelblue",
         xlab = "", ylab = "", bty = "L",
         main = "k-TWI")
    grid(col = "gray50", lty = "dotted")
  }
  
}


# Analysis using residuals
ar_res = function (x, order = 12) {
  model = ar(x, aic = F, order.max = order, na.action = na.pass)
  
  return (model$resid)
}

acf_type = "covariance"
for (pct in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  
  cat("\n\npercentage:", pct, "\n")
  max_order = 3    
  if (acf_type == "partial") {
    acf_hat = array(NA, dim = c(8, d, max_order))
  } else {
    acf_hat = array(NA, dim = c(8, d, max_order))
  }
  
  ground_truth = as.matrix(read.csv("./real_data/small_series2.csv", header = FALSE))
  ground_truth = t(ground_truth)
  
  linimp = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_lin.csv"), 
                              header = FALSE))
  linimp = t(linimp)
  
  SSA = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SSA.csv"), 
                           header = FALSE))
  SSA = t(SSA)
  
  KS = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_KS.csv"), 
                          header = FALSE))
  KS = t(KS)
  
  SF = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_SF.csv"), 
                          header = FALSE))
  SF = t(SF)
  
  mean_interpolation = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_mean.csv"), 
                                          header = FALSE))
  mean_interpolation = t(mean_interpolation)
  
  TWI = as.matrix(read.csv(paste0("./real_data/WI_small", 10 * pct, ".csv"), 
                           header = FALSE))
  kTWI = as.matrix(read.csv(paste0("./real_data/kWI_small", 10 * pct, ".csv"), 
                            header = FALSE))
  LP_trend = as.matrix(read.csv(paste0("./real_data/small_series", 10 * pct, "_LPfit.csv"), 
                                header = FALSE))
  LP_trend = t(LP_trend)
  
  TWI = TWI + LP_trend
  kTWI = kTWI + LP_trend
  
  
  for (j in 1:d) {
    if (acf_type != "partial") {
      acf_hat[1,j,] = acf(ar_res(ground_truth[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[1,j,] = acf(ar_res(ground_truth[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[2,j,] = acf(ar_res(linimp[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[2,j,] = acf(ar_res(linimp[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[3,j,] = acf(ar_res(SSA[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[3,j,] = acf(ar_res(SSA[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[4,j,] = acf(ar_res(KS[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[4,j,] = acf(ar_res(KS[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[5,j,] = acf(ar_res(SF[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[5,j,] = acf(ar_res(SF[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[6,j,] = acf(ar_res(mean_interpolation[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[6,j,] = acf(ar_res(mean_interpolation[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[7,j,] = acf(ar_res(TWI[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[7,j,] = acf(ar_res(TWI[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
    if (acf_type != "partial") {
      acf_hat[8,j,] = acf(ar_res(kTWI[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf[2:(1 + max_order)]
    } else {
      acf_hat[8,j,] = acf(ar_res(kTWI[,j], 1), 
                          lag.max = max_order, na.action = na.pass, plot = F,
                          type = acf_type)$acf
    }
    
  }
  
  # ACF analysis
  methods = c("Ground truth", "Linear", "iSSA", "Kalman smoothing", "Scalar filter",
              "Mean interpolation", "TWI", "k-TWI")
  for (j in 2:8) {
    if (j < 8) {
      if (j == 2) {
        cat("Mean: ")
      }
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,])) * 10, 2), "& ")
    } else {
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,])) * 10, 2), "\n")
    }
  }
  
  for (j in 2:8) {
    if (j < 8) {
      if (j == 2) {
        cat("Median: ")
      }
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,]), "median") * 10, 2), "& ")
    } else {
      cat(round(mae(c(acf_hat[1,,]), c(acf_hat[j,,]), "median") * 10, 2), "\n")
    }
  }
  

  for (h in 1:max_order) {
    cat("lag =", h, "\n")
    for (k in 2:8) {
      cat(round(mae(acf_hat[1,,h], acf_hat[k,,h]) * 10, 2), "& ")
      if (k == 8) {
        cat("\n")
      }
    }
  }
  
}
