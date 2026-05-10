# Analysis using groundwater data

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

ground_truth = as.matrix(read.csv("./real_data/all_series/all_series.csv", header = FALSE))
ground_truth = t(ground_truth)
d = ncol(ground_truth)

all_pred_err = array(NA, dim = c(6, 7, d))
x_obs = ground_truth
linimp = as.matrix(read.csv(paste0("./real_data/all_series/all_series_lin.csv"), 
                            header = FALSE))
linimp = t(linimp)
SSA = as.matrix(read.csv(paste0("./real_data/all_series/all_series_SSA.csv"), 
                         header = FALSE))
SSA = t(SSA)
KS = as.matrix(read.csv(paste0("./real_data/all_series/all_series_KS.csv"), 
                        header = FALSE))
KS = t(KS)
SF = as.matrix(read.csv(paste0("./real_data/all_series/all_series_SF.csv"), 
                        header = FALSE))
SF = t(SF)
mean_interpolation = as.matrix(read.csv(paste0("./real_data/all_series/all_series_mean.csv"), 
                                        header = FALSE))
mean_interpolation = t(mean_interpolation)
TWI = as.matrix(read.csv(paste0("./real_data/all_series/WI_all_series.csv"), 
                         header = FALSE))
kTWI = as.matrix(read.csv(paste0("./real_data/all_series/kWI_all_series.csv"), 
                          header = FALSE))
mTWI = as.matrix(read.csv(paste0("./real_data/all_series/WI_all_series_multi.csv"), 
                          header = FALSE))
mkTWI = as.matrix(read.csv(paste0("./real_data/all_series/kWI_all_series_multi.csv"), 
                            header = FALSE))

LP_trend = as.matrix(read.csv(paste0("./real_data/all_series/all_series_LPfit.csv"), 
                              header = FALSE))
LP_trend = t(LP_trend)

pred_err = array(NA, dim = c(9, d))
training_n = 200
order = 12
for (j in 1:d) {
  pred_err[1,j] = AR_rolling(linimp[,j] - LP_trend[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[2,j] = AR_rolling(SSA[,j]- LP_trend[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[3,j] = AR_rolling(KS[,j] - LP_trend[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[4,j] = AR_rolling(SF[,j] - LP_trend[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[5,j] = AR_rolling(mean_interpolation[,j], ground_truth[,j], order, training_n)
  pred_err[6,j] = AR_rolling(TWI[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[7,j] = AR_rolling(kTWI[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[8,j] = AR_rolling(mTWI[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
  pred_err[9,j] = AR_rolling(mkTWI[,j], ground_truth[,j] - LP_trend[,j], order, training_n)
}

missing_percentage = rep(NA, ncol(ground_truth))
for (j in 1:d) {
  missing_percentage[j] = mean(is.na(ground_truth[,j]))
}

par(mfrow = c(3, 3))
plot(x = missing_percentage * 100, y = pred_err[1,], xlab = "Percent", ylab = "", 
     main = "Linear", ylim = c(0, max(pred_err) * 1.1))
plot(x = missing_percentage * 100, y = pred_err[2,], xlab = "Percent", ylab = "", 
     main = "SSA", ylim = c(0, max(pred_err) * 1.1))
plot(x = missing_percentage * 100, y = pred_err[3,], xlab = "Percent", ylab = "", 
     main = "Kalman smoothing", ylim = c(0, max(pred_err) * 1.1))
# plot(x = missing_percentage * 100, y = pred_err[4,], xlab = "Percent", ylab = "",
#      main = "Scalar filtering", ylim = c(0, max(pred_err) * 1.1))
# plot(x = missing_percentage * 100, y = pred_err[5,], xlab = "Percent", ylab = "",
#      main = "mean_interpolation", ylim = c(0, max(pred_err) * 1.1))
plot(x = missing_percentage * 100, y = pred_err[6,], xlab = "Percent", ylab = "", 
     main = "TWI (linear)", ylim = c(0, max(pred_err) * 1.1))
plot(x = missing_percentage * 100, y = pred_err[7,], xlab = "Percent", ylab = "", 
     main = "k-TWI (linear)", ylim = c(0, max(pred_err) * 1.1))
plot(x = missing_percentage * 100, y = pred_err[8,], xlab = "Percent", ylab = "", 
     main = "multivariate TWI (linear)", ylim = c(0, max(pred_err) * 1.1))
plot(x = missing_percentage * 100, y = pred_err[9,], xlab = "Percent", ylab = "", 
     main = "multivariate k-TWI (linear)", ylim = c(0, max(pred_err) * 1.1))

par(mfrow = c(1, 1))
boxplot(t(pred_err)[,c(-5, -4)], names = c("Linear", "SSA", "Kalman smoothing", "TWI (linear)", "k-TWI (linear)",
                                           "multi-TWI (linear)", "multi-kTWI (linear)"))

par(mfcol = c(3, 4))
plot(x = missing_percentage * 100, y = pred_err[1,] - pred_err[6,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Linear interpolation - TWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[2,] - pred_err[6,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "SSA - TWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[3,] - pred_err[6,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Kalman smoothing - TWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")

plot(x = missing_percentage * 100, y = pred_err[1,] - pred_err[7,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Linear interpolation - kTWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[2,] - pred_err[7,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "SSA - kTWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[3,] - pred_err[7,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Kalman smoothing - kTWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")

plot(x = missing_percentage * 100, y = pred_err[1,] - pred_err[8,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Linear interpolation - multi-TWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[2,] - pred_err[8,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "SSA - multi-TWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[3,] - pred_err[8,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Kalman smoothing - multi-TWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")

plot(x = missing_percentage * 100, y = pred_err[1,] - pred_err[9,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Linear interpolation - multi-kTWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[2,] - pred_err[9,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "SSA - multi-kTWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
plot(x = missing_percentage * 100, y = pred_err[3,] - pred_err[9,], xlab = "Percentage of Missing Data",
     ylab = "", col = "steelblue", main = "Kalman smoothing - multi-kTWI", ylim = c(-0.06, 0.18), pch = 19)
abline(h = 0, col = "firebrick4")
