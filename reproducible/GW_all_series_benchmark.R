### Groundwater analysis: All series
library(glmnet)
library(imputeTS)
library(Rssa)
library(locfit)

s_filter = function (x, ar_order = 2, lambda = 0.001) {
  x = as.matrix(x)
  k = ncol(x)
  
  temp = rowSums(is.na(x))
  miss_idx = which(temp > 0)
  
  iv = matrix(0, nrow(x), length(miss_idx))
  for (i in 1:length(miss_idx)) {
    iv[miss_idx[i], i] = 1
  }
  
  xmis = x
  xmis[miss_idx,] = 0 # rnorm(length(miss_idx))
  
  y_aug = embed(xmis, ar_order + 1)
  y = y_aug[,1:k]
  ylag = y_aug[,-c(1:k)]
  
  model = glmnet(x = cbind(ylag, iv[(ar_order + 1):nrow(iv),]),
                 y = y, 
                 alpha = 0,
                 family = "mgaussian",
                 lambda = lambda)
  
  if (k > 1) {
    est_imp = matrix(0, nrow = length(miss_idx), ncol = k)
    for (i in 1:length(miss_idx)) {
      for (j in 1:k) {
        if (j == 1) {
          temp = c(as.matrix(model$beta[[j]]))[(ar_order * k) + i]
        } else {
          temp = c(temp, c(as.matrix(model$beta[[j]]))[(ar_order * k) + i])
        }
      }
      est_imp[i,] = -temp
    }
    xmis[miss_idx,] = est_imp
  } else {
    xmis[miss_idx,] = -c(model$beta[(ar_order + 1):length(model$beta)])
  }
  
  
  return (xmis)
}

issa = function (x, d, L = 6, r = 3) {
  if (d > 1) {
    x = as.matrix(x)
    model = ssa(x, L = L, kind = "mssa")
  } else {
    x = as.vector(x)
    model = ssa(x, L = L)
  }
  
  res = igapfill(model, groups = list(c(1:r)), maxiter = 100,
                 fill = na_mean(x, "mean"))
  
  return (res)
}

ssa_cv = function (x, max_rank = 6, L = 6, mask_prop = 0.1, nrep = 10,
                   preserve = NULL) {
  # x must be vector
  x_obs = x
  
  obs_idx = which(!is.na(x_obs))
  n_mask = ceiling(length(obs_idx) * mask_prop)
  
  errors = rep(NA, max_rank)
  
  for (r in 1:max_rank) {
    cv_rmse = rep(NA, nrep)
    for (j in 1:nrep) {
      if (is.null(preserve)) {
        mask_idx = sample(obs_idx, n_mask)
      } else {
        mask_idx = sample(obs_idx[which(obs_idx < preserve)], n_mask)
      }
      x_masked = x_obs
      x_masked[mask_idx] = NA
      
      x_filled = issa(x_masked, 1, L = L, r = r)
      
      cv_rmse[j] = sqrt(mean((x_obs[mask_idx] - x_filled[mask_idx])^2, na.rm = TRUE))
    }
    
    errors[r] = mean(cv_rmse)
  }
  
  optimal_r = which.min(errors)
  
  x_filled = issa(x, 1, L = L, r = optimal_r)
  
  return (list("errors" = errors, "x_filled" = x_filled, "optimal_r" = optimal_r))
}

dta = unname(as.matrix(read.csv("./real_data/GW_select.csv", header = F)))
dta = (dta - rowMeans(dta, na.rm = TRUE)) / apply(dta, 1, sd, na.rm = TRUE)

# Clip data larger than 5 standard deviations
for (j in 1:nrow(dta)) {
  for (i in 1:ncol(dta)) {
    if (is.na(dta[j,i])) {
      next
    } else if (abs(dta[j,i]) > 5) {
      dta[j,i] = NA
    }
  }
}

# de-trend each row
detrend = function (x) {
  n = length(x)
  model = lm(x~c(1:n))
  b0 = model$coefficients[1]
  b1 = model$coefficients[2]
  
  return (x - b0 - b1 * c(1:n))
}
for (j in 1:nrow(dta)) {
  dta[j,] = detrend(dta[j,])
}

# Focus on data after 1995:01
# New data dimension: (176, 308)
dta = dta[, 28:ncol(dta)]

# h <- hist(rowMeans(is.na(dta)), breaks = 12, plot = F)
# h$counts <- 100 * h$counts / sum(h$counts)
# plot(h, ylab = "Percent of all sites", xlab = "Percentage of missing data in the sample",
#      main = "")
# 
# plot(x = seq(from = as.Date("1995-01-01"), to = as.Date("2020-08-01"), by = "month"),
#      y = colMeans(is.na(dta)) * 100,
#      type = "l", col = "steelblue", lwd = 1.5, xaxt = "n", 
#      xlab = "Time", ylab = "Percent",
#      main = "Percent of all cites that record a missing data")
# axis.Date(1, at = seq(from = as.Date("1995-01-01"), to = as.Date("2020-08-01"), by = "5 years"), format = "%Y")

### Benchmark methods
lin_imp = array(NA, dim = dim(dta))
KS_imp = array(NA, dim = dim(dta))
SF_imp = array(NA, dim = dim(dta))
SSA_imp = array(NA, dim = dim(dta))
mean_imp = array(NA, dim = dim(dta))
lp_trend = array(NA, dim = dim(dta))
lp_resd = array(NA, dim = dim(dta))
lp_linimp = array(NA, dim = dim(dta))

for (i in 1:nrow(dta)) {
  cat("Imputing series", i, "\n")
  
  dta_slice = dta[i,]
  
  lin_imp[i,] = na_interpolation(dta_slice, "linear")
  KS_imp[i,] = na_kalman(dta_slice, "auto.arima")
  SF_imp[i,] = s_filter(dta_slice, 6)
  SSA_imp[i,] = ssa_cv(dta_slice, 6, 6, 0.1, 10, ncol(dta) - 35)$x_filled
  mean_imp[i,] = na_mean(dta_slice, "mean")
  
  # Remove local polynomial
  model = locfit(dta_slice[which(!is.na(dta_slice))] ~ lp(which(!is.na(dta_slice)), deg = 0, h = 12))
  res = c(dta_slice - predict(model, newdata = 1:length(dta_slice)))
  
  lp_trend[i,] = predict(model, newdata = 1:length(dta_slice))
  lp_resd[i,] = res
  lp_linimp[i,] = na_interpolation(res, "linear")
  
  # save results
  write.table(dta, paste0("./real_data/all_series/all_series", ".csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(lin_imp, paste0("./real_data/all_series/all_series", "_lin.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(KS_imp, paste0("./real_data/all_series/all_series", "_KS.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(SF_imp, paste0("./real_data/all_series/all_series", "_SF.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(SSA_imp, paste0("./real_data/all_series/all_series", "_SSA.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(mean_imp, paste0("./real_data/all_series/all_series", "_mean.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(lp_trend, paste0("./real_data/all_series/all_series", "_LPfit.csv"),
              row.names = F, col.names = F, sep = ",")
  write.table(lp_resd, paste0("./real_data/all_series/all_series", "_LPresd.csv"),
              row.names = F, col.names = F, sep = ",")
  write.table(lp_linimp, paste0("./real_data/all_series/all_series", "_LPlinimp.csv"),
              row.names = F, col.names = F, sep = ",")
}

