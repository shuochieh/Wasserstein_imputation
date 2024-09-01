library(imputeTS)
library(MTS)
library(ggplot2)
library(gridExtra)
library(glmnet)

# Remember to specify model_name
# model_name = "ARMA"
# d = 1 (dimension of the target series)
# n_exper = 1000 (# of monte carlo simulations)


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

### Missing pattern 1
dta_com = read.csv(paste0("./sim_data/", model_name, "/", model_name, "_original.csv"), header = F)
dta = read.csv(paste0("./sim_data/", model_name, "/", model_name, "_miss1.csv"), header = F)

# Exploratory
df = data.frame(x = c(dta_com[1:(nrow(dta_com) - 1),1]), 
                y = c(dta_com[2:nrow(dta_com),1]))
p1 = ggplot(df, aes(x = x, y = y)) + 
  geom_point(color = "steelblue", size = 1) +
  theme_minimal() +
  labs(title = "Scatterplot of original data",
       x = "x_{t-1}",
       y = "x_t")

plot_range = range(dta_com[,1])
plot_range = c(floor(plot_range[1]), ceiling(plot_range[2]))

# linear interpolation
cat("\n======= Linear interpolation ==========================\n")
lin_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  lin_imp[,((sim - 1) * d + 1):(sim * d)] = as.matrix(na_interpolation(dta[,((sim - 1) * d + 1):(sim * d)], "linear"))
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

# spline interpolation
cat("\n======= Spline interpolation ==========================\n")
spl_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  spl_imp[,((sim - 1) * d + 1):(sim * d)] = as.matrix(na_interpolation(dta[,((sim - 1) * d + 1):(sim * d)], "spline"))
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

# Kalman smoothing (auto.arima)
cat("\n======= Kalman smoothing ==========================\n")
KS_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  KS_imp[,((sim - 1) * d + 1):(sim * d)] = as.matrix(na_kalman(dta[,((sim - 1) * d + 1):(sim * d)], "auto.arima"))
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

# scalar filtering
cat("\n======= Pena-Tsay filter ==========================\n")
SF_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  SF_imp[,((sim - 1) * d + 1):(sim * d)] = s_filter(dta[,((sim - 1) * d + 1):(sim * d)], 6, lambda = 0.0001) 
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(lin_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p2 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of linear imputation",
       x = "x_{t-1}",
       y = "x_t")

temp = na_lag_pairs(spl_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p3 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of spline imputation",
       x = "x_{t-1}",
       y = "x_t")

temp = na_lag_pairs(KS_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p4 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of Kalman smoothing",
       x = "x_{t-1}",
       y = "x_t")

temp = na_lag_pairs(SF_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p5 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of PT filter",
       x = "x_{t-1}",
       y = "x_t")

grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 3)

output = cbind(lin_imp, spl_imp, KS_imp, SF_imp)
write.csv(output, paste0("./sim_data/", model_name, "/benchmarks_miss1.csv"), row.names = F)

par(mfrow = c(1,3))
plot(spl_imp[1:300,1], type = "l")
plot(KS_imp[1:300,1], type = "l")
plot(SF_imp[1:300,1], type = "l")







### Missing pattern 2
dta_com = read.csv(paste0("./sim_data/", model_name, "/", model_name, "_original.csv"), header = F)
dta = read.csv(paste0("./sim_data/", model_name, "/", model_name, "_miss2.csv"), header = F)

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
cat("\n======= Linear interpolation ==========================\n")
lin_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  lin_imp[,((sim - 1) * d + 1):(sim * d)] = as.matrix(na_interpolation(dta[,((sim - 1) * d + 1):(sim * d)], "linear"))
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

# spline interpolation
cat("\n======= Spline interpolation ==========================\n")
spl_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  spl_imp[,((sim - 1) * d + 1):(sim * d)] = as.matrix(na_interpolation(dta[,((sim - 1) * d + 1):(sim * d)], "spline"))
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

# Kalman smoothing (auto.arima)
cat("\n======= Kalman smoothing ==========================\n")
KS_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  KS_imp[,((sim - 1) * d + 1):(sim * d)] = as.matrix(na_kalman(dta[,((sim - 1) * d + 1):(sim * d)], "auto.arima"))
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

# scalar filtering
cat("\n======= Pena-Tsay filter ==========================\n")
SF_imp = matrix(NA, nrow = nrow(dta), ncol = ncol(dta))
for (sim in 1:n_exper) {
  SF_imp[,((sim - 1) * d + 1):(sim * d)] = s_filter(dta[,((sim - 1) * d + 1):(sim * d)], 6, lambda = 0.001) # AR order may not be too big
  if (sim %% 100 == 0 && sim > 99) {
    cat("Iteration", sim, "\n")
  }
}

temp = na_lag_pairs(lin_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p2 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of linear imputation",
       x = "x_{t-1}",
       y = "x_t")

temp = na_lag_pairs(spl_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p3 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of spline imputation",
       x = "x_{t-1}",
       y = "x_t")

temp = na_lag_pairs(KS_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p4 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of Kalman smoothing",
       x = "x_{t-1}",
       y = "x_t")

temp = na_lag_pairs(SF_imp[,1], which(is.na(dta[,1])))
temp = data.frame(x = temp[,1],
                  y = temp[,2])
p5 = ggplot(temp, aes(x = x, y = y)) + 
  geom_point(color = "steelblue1", size = 1) +
  coord_cartesian(xlim = plot_range, ylim = plot_range) +
  theme_minimal() +
  labs(title = "Scatterplot of PT filter",
       x = "x_{t-1}",
       y = "x_t")

grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 3)

par(mfrow = c(1,3))
plot(spl_imp[,1], type = "l")
plot(KS_imp[,1], type = "l")
plot(SF_imp[,1], type = "l")

output = cbind(lin_imp, spl_imp, KS_imp, SF_imp)
write.csv(output, paste0("./sim_data/", model_name, "/benchmarks_miss2.csv"), row.names = F)






