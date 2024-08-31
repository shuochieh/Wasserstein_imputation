library(imputeTS)
library(MTS)
library(ggplot2)
library(gridExtra)
library(glmnet)

# Some utilities
na_lag_pairs = function (x, missing_idx) {
  res = matrix(NA, nrow = length(missing_idx) * 2, ncol = 2)
  counter = 1
  for (i in missing_idx) {
    res[counter,] = c(x[i - 1], x[i])
    res[counter + 1,] = c(x[i], x[i + 1])
    counter = counter + 2
  }
  if (counter != (2 * length(missing_idx) + 1)) {
    cat("Counter:", counter, "; num missing", length(missing_idx), "\n")
  }
  
  return(res)
}

### scalar filter algorithm of Pena and Tsay (2021)
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
                 lambda = lambda)

  xmis[miss_idx,] = -model$beta[-c(1:ar_order)]
  

  # return (list("xmis" = xmis, "iv" = iv)) 
  return (xmis)
}

### Missing pattern 1
dta_com = read.csv("./sim_data/AR/AR_original.csv", header = F)
dta = read.csv("./sim_data/AR/AR_miss1.csv", header = F)

# Exploratory
df = data.frame(x = c(dta_com[1:(nrow(dta_com) - 1),1]), 
                y = c(dta_com[2:nrow(dta_com),1]))
p1 = ggplot(df, aes(x = x, y = y)) + 
  geom_point(color = "steelblue", size = 1) +
  theme_minimal() +
  labs(title = "Scatterplot of original data",
       x = "x_{t-1}",
       y = "x_t")

# linear interpolation
lin_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  lin_imp[,sim] = na_interpolation(dta[,sim], "linear")
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(lin_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p2 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of linear imputation",
       x = "x_{t-1}",
       y = "x_t")

# spline interpolation
spl_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  spl_imp[,sim] = na_interpolation(dta[,sim], "spline")
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(spl_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p3 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of spline imputation",
       x = "x_{t-1}",
       y = "x_t")

# Kalman smoothing (auto.arima)
KS_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  KS_imp[,sim] = na_kalman(dta[,sim], "auto.arima")
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(KS_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p4 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of Kalman smoothing",
       x = "x_{t-1}",
       y = "x_t")

# scalar filtering
SF_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  SF_imp[,sim] = s_filter(dta[,sim], 6, lambda = 0.001) # AR order may not be too big
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(SF_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p5 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of PT filter",
       x = "x_{t-1}",
       y = "x_t")

grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 3)

output = cbind(lin_imp, spl_imp, KS_imp, SF_imp)
write.csv(output, "./sim_data/AR/benchmarks_miss1.csv", row.names = F)











### Missing pattern 2
dta_com = read.csv("./sim_data/AR/AR_original.csv", header = F)
dta = read.csv("./sim_data/AR/AR_miss2.csv", header = F)

# Exploratory
df = data.frame(x = c(dta_com[1:(nrow(dta_com) - 1),1]), 
                y = c(dta_com[2:nrow(dta_com),1]))
p1 = ggplot(df, aes(x = x, y = y)) + 
  geom_point(color = "steelblue", size = 1) +
  theme_minimal() +
  labs(title = "Scatterplot of original data",
       x = "x_{t-1}",
       y = "x_t")

# linear interpolation
lin_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  lin_imp[,sim] = na_interpolation(dta[,sim], "linear")
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(lin_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p2 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of linear imputation",
       x = "x_{t-1}",
       y = "x_t")

# spline interpolation
spl_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  spl_imp[,sim] = na_interpolation(dta[,sim], "spline")
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(spl_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p3 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of spline imputation",
       x = "x_{t-1}",
       y = "x_t")

# Kalman smoothing (auto.arima)
KS_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  KS_imp[,sim] = na_kalman(dta[,sim], "auto.arima")
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(KS_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p4 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of Kalman smoothing",
       x = "x_{t-1}",
       y = "x_t")

# scalar filtering
SF_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:ncol(dta)) {
  SF_imp[,sim] = s_filter(dta[,sim], 6, lambda = 0.001) # AR order may not be too big
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(SF_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p5 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_minimal() +
  labs(title = "Scatterplot of PT filter",
       x = "x_{t-1}",
       y = "x_t")

grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 3)

output = cbind(lin_imp, spl_imp, KS_imp, SF_imp)
write.csv(output, "./sim_data/AR/benchmarks_miss2.csv", row.names = F)






