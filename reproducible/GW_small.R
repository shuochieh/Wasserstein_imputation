### Groundwater analysis: Small dataset benchmark
library(glmnet)
library(imputeTS)
library(Rssa)

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

# Focus on data after 1995:01, and on the subset of stations with less than < 10 missing data;
# New data dimension: (15, 308)
subset_idx = which(rowSums(is.na(dta[,28:ncol(dta)])) < 10)
dta = dta[subset_idx, 28:ncol(dta)]

#subset_idx = which(rowSums(is.na(dta)) < 30) # stations with less than < 30 missing data;
# New data dimension: (19, 308)
#dta = dta[subset_idx, 28:ncol(dta)]

par(mfrow = c(3,5))
dates = seq(from = as.Date("1995-01-01"), to = as.Date("2020-08-01"), by = "month")
for (i in 1:15) {
  plot(x = dates, y = dta[i,], type = "l", col = "steelblue", lwd = 1.5,
       xaxt = "n", yaxt = "n",
       main = paste("Series", i), xlab = "", ylab = "")
  axis.Date(1, at = seq(as.Date("1995-01-01"), as.Date("2020-01-01"), by = "5 years"), format = "%Y")
  axis(2, las = 1, cex.axis = 0.7)
  
  # Add grid
  grid(col = "gray90", lty = "dotted")
}

### Artificially create missing data
set.seed(1)
missing_pct = c(0.2, 0.3, 0.4, 0.5)

for (pct in missing_pct) {
  dta_artificial = dta
  lin_imp = array(NA, dim = dim(dta))
  KS_imp = array(NA, dim = dim(dta))
  SF_imp = array(NA, dim = dim(dta))
  SSA_imp = array(NA, dim = dim(dta))
  mean_imp = array(NA, dim = dim(dta))
  
  for (i in 1:15) {
    cat("Imputing...missing percentage", 100 * pct, "%. Series", i, "\n")
    dta_artificial[i,sample(ncol(dta) - 36, floor(pct * ncol(dta)))] = NA
    lin_imp[i,] = na_interpolation(dta_artificial[i,], "linear")
    KS_imp[i,] = na_kalman(dta_artificial[i,], "auto.arima")
    SF_imp[i,] = s_filter(dta_artificial[i,], 6)
    SSA_imp[i,] = ssa_cv(dta_artificial[i,], 6, 6, 0.1, 10, ncol(dta) - 35)$x_filled
    mean_imp[i,] = na_mean(dta_artificial[i,], "mean")
  }
  
  write.table(dta, paste0("./real_data/small_series", pct * 10, ".csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(dta_artificial, paste0("./real_data/small_series", pct * 10, "_masked.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(lin_imp, paste0("./real_data/small_series", pct * 10, "_lin.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(KS_imp, paste0("./real_data/small_series", pct * 10, "_KS.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(SF_imp, paste0("./real_data/small_series", pct * 10, "_SF.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(SSA_imp, paste0("./real_data/small_series", pct * 10, "_SSA.csv"), 
              row.names = F, col.names = F, sep = ",")
  write.table(mean_imp, paste0("./real_data/small_series", pct * 10, "_mean.csv"), 
              row.names = F, col.names = F, sep = ",")
}