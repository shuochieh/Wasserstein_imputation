library(glmnet)
library(imputeTS)

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

# linear interpolation
lin_imp = as.matrix(na_interpolation(t(dta), "linear"))

# spline interpolation
spl_imp = as.matrix(na_interpolation(t(dta), "spline"))

# Kalman smoothing
KS_imp = as.matrix(na_kalman(t(dta), "auto.arima"))

# scalar filtering
SF_imp = matrix(NA, ncol = nrow(dta), nrow = ncol(dta))
for (j in 1:nrow(dta)) {
  SF_imp[,j] = s_filter(dta[j,], 6)
}

# Save
write.table(lin_imp, "./real_data/lin_imp.csv", row.names = F, col.names = F, sep = ",")
write.table(spl_imp, "./real_data/spl_imp.csv", row.names = F, col.names = F, sep = ",")
write.table(KS_imp, "./real_data/KS_imp.csv", row.names = F, col.names = F, sep = ",")
write.table(SF_imp, "./real_data/SF_imp.csv", row.names = F, col.names = F, sep = ",")
write.table(dta, "./real_data/GW_clip.csv", row.names = F, col.names = F, sep = ",")