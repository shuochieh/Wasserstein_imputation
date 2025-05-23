library(TSA)
library(glmnet)
library(reshape2)
library(ggplot2)
library(fastICA)

LYB_FM = function (x, r, h, demean = TRUE) {
  # x: (n by d) matrix
  # r: number of factors
  # h: maximum lags to use
  
  if (demean) {
    x = t(t(x) - colMeans(x))
  }
  
  pd = 0
  H = h + 1
  n = nrow(x)
  for (i in 1:h) {
    temp = (t(x[H:n,]) %*% x[(H - i):(n - i),]) / n
    pd = pd + temp %*% t(temp)
  }
  
  # Eigenanalysis
  Evec = eigen(pd)$vectors
  V = Evec[,1:r]
  
  # Extract factors and residuals
  f_hat = x %*% V
  e_hat = x - f_hat %*% t(V)
  
  return (list("V" = V, "f_hat" = f_hat, "e_hat" = e_hat))
}

AR_test = function (x, n1, missing_idx, p = 6) {
  # x: a vector of time series of length n
  # missing_idx: a vector of length (n - n1); equals one if missing in 
  #              the original data
  train_x = x[1:n1]
  test_x = x[-c(1:n1)]
  
  aux = embed(train_x, p + 1)
  model = lm(aux[,1]~aux[,-1])
  pred = rep(NA, length(test_x))
  for (i in 1:length(test_x)) {
    if (missing_idx[i]) {
      next
    } else {
      pred[i] = c(1, x[(n1 + i - 1):(n1 + i - p)]) %*% c(model$coefficients)
      if (is.na(pred[i])) {
        pred[i] = mean(train_x)
      }
    }
  }
  
  return (pred)
}

raw_clip = t(unname(as.matrix(read.csv("./real_data/GW_clip.csv", header = F))))

lin_imp = unname(as.matrix(read.csv("./real_data/lin_imp.csv", header = F)))
spl_imp = unname(as.matrix(read.csv("./real_data/spl_imp.csv", header = F)))
KS_imp = unname(as.matrix(read.csv("./real_data/KS_imp.csv", header = F)))
SF_imp = unname(as.matrix(read.csv("./real_data/SF_imp.csv", header = F)))

WI_lin = unname(as.matrix(read.csv("./real_data/WI_lin.csv", header = F)))
kWI_lin = unname(as.matrix(read.csv("./real_data/kWI_lin.csv", header = F)))
WI_Kalman = unname(as.matrix(read.csv("./real_data/WI_Kalman.csv", header = F)))
kWI_Kalman = unname(as.matrix(read.csv("./real_data/kWI_Kalman.csv", header = F)))

# Factor analysis
par(mfrow = c(4, 2), mar = c(2.5, 2, 1, 0.5))
r = 2

# linear 
lin_fac = LYB_FM(lin_imp, r, 6)
temp = lin_fac$f_hat
temp[,1] = -temp[,1]
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(linear)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# spline
spl_fac = LYB_FM(spl_imp, r, 6)
temp = spl_fac$f_hat
temp[,1] = -temp[,1]
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(spline)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# Kalman
KS_fac = LYB_FM(KS_imp, r, 6)
temp = KS_fac$f_hat
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(Kalman)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# Scalar filter
SF_fac = LYB_FM(SF_imp, r, 6)
temp = SF_fac$f_hat
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(ScalarF)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# WI (linear)
WI_lin_fac = LYB_FM(WI_lin, r, 6)
temp = WI_lin_fac$f_hat
temp[,2] = -temp[,2]
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(TWI - linear)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# kWI (linear)
kWI_lin_fac = LYB_FM(kWI_lin, r, 6)
temp = kWI_lin_fac$f_hat
temp = -temp
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(k-TWI - linear)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# WI (Kalman)
WI_Kalman_fac = LYB_FM(WI_Kalman, r, 6)
temp = WI_Kalman_fac$f_hat
temp[,1] = -temp[,1]
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(TWI - Kalman)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# kWI (Kalman)
kWI_Kalman_fac = LYB_FM(kWI_Kalman, r, 6)
temp = kWI_Kalman_fac$f_hat
temp = -temp
for (ii in 1:r) {
  plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
       y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
       main = paste("Estimated factor", ii, "(k-TWI - Kalman)"),
       ylim = range(temp, na.rm = TRUE))
  abline(v = as.Date("2002-01-01"), col = "red", lty = 2)
}

# ICA analysis
dev.off()
par(mfrow = c(8, 2), mar = c(2.5, 2, 1, 0.5))

add_shading <- function(start_month, end_month, x_vals, color = rgb(0.8, 0.8, 0.8, 0.5)) {
  years <- unique(format(x_vals, "%Y"))
  for (y in years) {
    rect(xleft = as.Date(paste0(y, "-", start_month, "-01")),
         xright = as.Date(paste0(y, "-", end_month, "-30")),
         ybottom = par("usr")[3], 
         ytop = par("usr")[4], 
         col = color, border = NA)
  }
}

set.seed(5566)
# Kalman
ica = fastICA(KS_imp, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp = -temp
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(Kalman)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(Kalman)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))

# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(Kalman)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }

# WI (Kalman)
ica = fastICA(WI_Kalman, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(TWI - Kalman)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(TWI - Kalman)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))

# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(TWI - Kalman)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }

# kWI (Kalman)
ica = fastICA(kWI_Kalman, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp[,1] = -temp[,1]
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(kTWI - Kalman)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(kTWI - Kalman)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))

# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(kTWI - Kalman)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }

# linear 
ica = fastICA(lin_imp, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp[,1] = -temp[,1]
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(linear)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(linear)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
# add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#             color = rgb(0.5, 0.5, 1, alpha = 0.3))

# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(linear)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }

# TWI (linear)
ica = fastICA(WI_lin, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp[,1] = -temp[,1]
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(TWI - linear)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(TWI - linear)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))

# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(TWI - linear)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }

# kTWI (linear)
ica = fastICA(kWI_lin, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp = -temp
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(kTWI - linear)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(kTWI - linear)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))
# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(kTWI - linear)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }


# spline 
ica = fastICA(spl_imp, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
#temp[,2] = -temp[,2]
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(spline)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(spline)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(spline)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }

# scalar filter
ica = fastICA(SF_imp, n.comp = 2)
temp = ica$S[,order(rowSums(ica$A^2), decreasing = T)]
temp[,1] = -temp[,1]
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,1], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 1, "(scalar filter)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(1, 0.5, 0.5, alpha = 0.3))
plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
     y = temp[,2], type = "l", lwd = 1, xlab = "", ylab = "",
     main = paste("Estimated IC", 2, "(scalar filter)"),
     ylim = range(temp, na.rm = TRUE))
abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
            color = rgb(0.5, 0.5, 1, alpha = 0.3))
# for (ii in 1:2) {
#   plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#        y = temp[,ii], type = "l", lwd = 1, xlab = "", ylab = "",
#        main = paste("Estimated IC", ii, "(scalar filter)"),
#        ylim = range(temp, na.rm = TRUE))
#   abline(v = as.Date("2002-01-01"), col = "red", lty = 2, lwd = 2)
#   add_shading("08", "09", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(1, 0.5, 0.5, alpha = 0.3))
#   add_shading("05", "06", seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
#               color = rgb(0.5, 0.5, 1, alpha = 0.3))
# }



# Prediction analysis
n_training = 120

# linear
pred_linear = matrix(0, nrow = nrow(lin_imp) - n_training, ncol = ncol(lin_imp))
for (j in 1:ncol(lin_imp)) {
  temp = AR_test(lin_imp[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_linear[,j] = temp
}

# spline
pred_spline = matrix(0, nrow = nrow(spl_imp) - n_training, ncol = ncol(spl_imp))
for (j in 1:ncol(spl_imp)) {
  temp = AR_test(spl_imp[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_spline[,j] = temp
}

# Kalman
pred_Kalman = matrix(0, nrow = nrow(KS_imp) - n_training, ncol = ncol(KS_imp))
for (j in 1:ncol(KS_imp)) {
  temp = AR_test(KS_imp[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_Kalman[,j] = temp
}

# Scalar filter
pred_SF = matrix(0, nrow = nrow(SF_imp) - n_training, ncol = ncol(SF_imp))
for (j in 1:ncol(SF_imp)) {
  temp = AR_test(SF_imp[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_SF[,j] = temp
}

# WI - linear
pred_WI_lin = matrix(0, nrow = nrow(WI_lin) - n_training, ncol = ncol(WI_lin))
for (j in 1:ncol(WI_lin)) {
  temp = AR_test(WI_lin[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_WI_lin[,j] = temp
}

# kWI - linear
pred_kWI_lin = matrix(0, nrow = nrow(kWI_lin) - n_training, ncol = ncol(kWI_lin))
for (j in 1:ncol(kWI_lin)) {
  temp = AR_test(kWI_lin[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_kWI_lin[,j] = temp
}

# WI - Kalman
pred_WI_Kalman = matrix(0, nrow = nrow(WI_Kalman) - n_training, ncol = ncol(WI_Kalman))
for (j in 1:ncol(WI_Kalman)) {
  temp = AR_test(WI_Kalman[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_WI_Kalman[,j] = temp
}

# kWI - Kalman
pred_kWI_Kalman = matrix(0, nrow = nrow(kWI_Kalman) - n_training, ncol = ncol(kWI_Kalman))
for (j in 1:ncol(kWI_Kalman)) {
  temp = AR_test(kWI_Kalman[,j], n_training, is.na(raw_clip[-c(1:n_training),j]))
  pred_kWI_Kalman[,j] = temp
}


errors = matrix(0, nrow = 8, ncol = 176)

# linear
errors[1,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_linear)^2, na.rm = T))

# spline
errors[2,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_spline)^2, na.rm = T))

# Kalman
errors[3,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_Kalman)^2, na.rm = T))

# SF
errors[4,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_SF)^2, na.rm = T))

# WI - lin
errors[5,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_WI_lin)^2, na.rm = T))

# kWI - lin
errors[6,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_kWI_lin)^2, na.rm = T))

# WI - Kalman
errors[7,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_WI_Kalman)^2, na.rm = T))

# kWI - Kalman
errors[8,] = sqrt(colMeans((raw_clip[-c(1:n_training),] - pred_kWI_Kalman)^2, na.rm = T))



method_names = c("Linear", "Spline", "Kalman", "Scalar filter", "TWI (lin)", "kTWI (lin)",
                 "TWI (Kal)", "kTWI (Kal)")
method_names = method_names[order(apply(errors, 1, median), decreasing = T)]
errors = errors[order(apply(errors, 1, median), decreasing = T),]
data_long <- melt(log(errors))
colnames(data_long) <- c("Method", "Observation", "RMSE")

# Convert Method to a factor with appropriate labels
data_long$Method <- factor(data_long$Method, 
                           labels = method_names)

# Create the violin plot
p <- ggplot(data_long, aes(x = Method, y = RMSE, fill = Method)) +
  geom_violin(trim = FALSE, color = "black") +  # Violin plot
  geom_boxplot(width = 0.1, position = position_dodge(0.9), fill = "white") +  
  scale_fill_brewer(palette = "Set3") +  
  labs(title = "", x = "", y = "log RMSE") +  
  theme_minimal(base_size = 18) +  
  theme(legend.position = "none")  

# Display the plot
print(p)

round(apply(t(errors), 2, median), 3)
round(apply(t(errors), 2, mean), 3)

summary(t(errors))
