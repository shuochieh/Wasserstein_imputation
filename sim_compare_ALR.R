library(transport)
library(ggplot2)
library(gridExtra)
library(imputeTS)
library(tidyr)

# Utils
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

ccf_summary = function (x) {
  temp11 = ccf(x[,1], x[,1], lag.max = 2, plot = F)
  temp22 = ccf(x[,2], x[,2], lag.max = 2, plot = F)
  temp33 = ccf(x[,3], x[,3], lag.max = 2, plot = F)
  temp12 = ccf(x[,1], x[,2], lag.max = 2, plot = F)
  temp13 = ccf(x[,1], x[,3], lag.max = 2, plot = F)
  temp23 = ccf(x[,2], x[,3], lag.max = 2, plot = F)
  
  res = array(0, dim = c(3, 3, 5))
  
  res[1,1,1] = c(temp11[0]$acf)
  res[2,2,1] = c(temp22[0]$acf)
  res[3,3,1] = c(temp33[0]$acf)
  res[1,2,1] = c(temp12[0]$acf)
  res[1,3,1] = c(temp13[0]$acf)
  res[2,3,1] = c(temp23[0]$acf)
  res[2,1,1] = c(temp12[0]$acf)
  res[3,1,1] = c(temp13[0]$acf)
  res[3,2,1] = c(temp23[0]$acf)
  
  res[1,1,2] = c(temp11[1]$acf)
  res[2,2,2] = c(temp22[1]$acf)
  res[3,3,2] = c(temp33[1]$acf)
  res[1,2,2] = c(temp12[1]$acf)
  res[1,3,2] = c(temp13[1]$acf)
  res[2,3,2] = c(temp23[1]$acf)
  res[2,1,2] = c(temp12[1]$acf)
  res[3,1,2] = c(temp13[1]$acf)
  res[3,2,2] = c(temp23[1]$acf)
  
  res[1,1,3] = c(temp11[2]$acf)
  res[2,2,3] = c(temp22[2]$acf)
  res[3,3,3] = c(temp33[2]$acf)
  res[1,2,3] = c(temp12[2]$acf)
  res[1,3,3] = c(temp13[2]$acf)
  res[2,3,3] = c(temp23[2]$acf)
  res[2,1,3] = c(temp12[2]$acf)
  res[3,1,3] = c(temp13[2]$acf)
  res[3,2,3] = c(temp23[2]$acf)
  
  res[1,1,4] = c(temp11[-1]$acf)
  res[2,2,4] = c(temp22[-1]$acf)
  res[3,3,4] = c(temp33[-1]$acf)
  res[1,2,4] = c(temp12[-1]$acf)
  res[1,3,4] = c(temp13[-1]$acf)
  res[2,3,4] = c(temp23[-1]$acf)
  res[2,1,4] = c(temp12[-1]$acf)
  res[3,1,4] = c(temp13[-1]$acf)
  res[3,2,4] = c(temp23[-1]$acf)
  
  res[1,1,5] = c(temp11[-2]$acf)
  res[2,2,5] = c(temp22[-2]$acf)
  res[3,3,5] = c(temp33[-2]$acf)
  res[1,2,5] = c(temp12[-2]$acf)
  res[1,3,5] = c(temp13[-2]$acf)
  res[2,3,5] = c(temp23[-2]$acf)
  res[2,1,5] = c(temp12[-2]$acf)
  res[3,1,5] = c(temp13[-2]$acf)
  res[3,2,5] = c(temp23[-2]$acf)
  
  return (res)
  # res0 = matrix(0, nrow = 3, ncol = 3)
  # res1 = matrix(0, nrow = 3, ncol = 3)
  # res2 = matrix(0, nrow = 3, ncol = 3)
  # res_1 = matrix(0, nrow = 3, ncol = 3)
  # res_2 = matrix(0, nrow = 3, ncol = 3)
  # 
  # res0[1,1] = c(temp11[0]$acf)
  # res0[2,2] = c(temp22[0]$acf)
  # res0[3,3] = c(temp33[0]$acf)
  # res0[1,2] = c(temp12[0]$acf)
  # res0[1,3] = c(temp13[0]$acf)
  # res0[2,3] = c(temp23[0]$acf)
  # res0[2,1] = c(temp12[0]$acf)
  # res0[3,1] = c(temp13[0]$acf)
  # res0[3,2] = c(temp23[0]$acf)
  # 
  # res1[1,1] = c(temp11[1]$acf)
  # res1[2,2] = c(temp22[1]$acf)
  # res1[3,3] = c(temp33[1]$acf)
  # res1[1,2] = c(temp12[1]$acf)
  # res1[1,3] = c(temp13[1]$acf)
  # res1[2,3] = c(temp23[1]$acf)
  # res1[2,1] = c(temp12[1]$acf)
  # res1[3,1] = c(temp13[1]$acf)
  # res1[3,2] = c(temp23[1]$acf)
  # 
  # res2[1,1] = c(temp11[2]$acf)
  # res2[2,2] = c(temp22[2]$acf)
  # res2[3,3] = c(temp33[2]$acf)
  # res2[1,2] = c(temp12[2]$acf)
  # res2[1,3] = c(temp13[2]$acf)
  # res2[2,3] = c(temp23[2]$acf)
  # res2[2,1] = c(temp12[2]$acf)
  # res2[3,1] = c(temp13[2]$acf)
  # res2[3,2] = c(temp23[2]$acf)
  # 
  # res_1[1,1] = c(temp11[-1]$acf)
  # res_1[2,2] = c(temp22[-1]$acf)
  # res_1[3,3] = c(temp33[-1]$acf)
  # res_1[1,2] = c(temp12[-1]$acf)
  # res_1[1,3] = c(temp13[-1]$acf)
  # res_1[2,3] = c(temp23[-1]$acf)
  # res_1[2,1] = c(temp12[-1]$acf)
  # res_1[3,1] = c(temp13[-1]$acf)
  # res_1[3,2] = c(temp23[-1]$acf)
  # 
  # res_2[1,1] = c(temp11[-2]$acf)
  # res_2[2,2] = c(temp22[-2]$acf)
  # res_2[3,3] = c(temp33[-2]$acf)
  # res_2[1,2] = c(temp12[-2]$acf)
  # res_2[1,3] = c(temp13[-2]$acf)
  # res_2[2,3] = c(temp23[-2]$acf)
  # res_2[2,1] = c(temp12[-2]$acf)
  # res_2[3,1] = c(temp13[-2]$acf)
  # res_2[3,2] = c(temp23[-2]$acf)
  # 
  # return (list(res_2, res_1, res0, res1, res2))
}

imp_plot = function (imp, x_obs, method) {
  imp = imp[1:300,]
  x_obs = x_obs[1:300,]
  
  data <- data.frame(
    index = 1:nrow(imp),
    y1 = imp[,1],
    y2 = imp[,2],
    y3 = imp[,3],
    x1_is_known = !is.na(x_obs[,1]),
    x2_is_known = !is.na(x_obs[,2]),
    x3_is_known = !is.na(x_obs[,3])
  )
  
  data_long <- pivot_longer(data, cols = c(y1, y2, y3), names_to = "series", values_to = "y_value")
  known_values <- pivot_longer(data, cols = c(x1_is_known, x2_is_known, x3_is_known), names_to = "series", values_to = "is_known")
  
  data_long$is_known <- known_values$is_known
  
  p = ggplot(data_long, aes(x = index, y = y_value, color = series, group = series)) +
    geom_line(size = 0.5) +  # Plot lines for each series
    # geom_point(data = subset(data_long, is_known == TRUE), shape = 16, size = 1.5, color = "steelblue") +  # Known values as blue points
    geom_point(data = subset(data_long, is_known == FALSE), shape = 18, size = 1.5, color = "indianred") +  # Imputed values as red diamonds
    scale_color_manual(values = c("y1" = "lightslategray", "y2" = "lightslategray", "y3" = "lightslategray")) +  # Custom colors for the lines
    theme_minimal() +  # Minimal theme
    labs(title = method, x = "", y = "Values") +  
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0.1, 0.6))
  
  return (p)
  
}

# ALR: I

model_name = "ALR"
base_dir = paste0("./sim_data/", model_name, "/")

truth = as.matrix(read.csv(paste0(base_dir, model_name, "_original.csv"), header = F))
x_obs = as.matrix(read.csv(paste0(base_dir, model_name, "_miss1.csv"), header = F))
benchmarks = as.matrix(read.csv(paste0(base_dir, "benchmarks_miss1.csv")))
WI_lin = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_lin_miss1.csv"), header = F))
WI_Kalman = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_Kalman_miss1.csv"), header = F))

wass_d = matrix(0, nrow = 1000, ncol = 8)
ccf_gt = array(0, dim = c(3, 3, 5))
ccfs = array(0, dim = c(3, 3, 5, 8))
ccf_loss = matrix(0, nrow = 5, ncol = 8)

for (sim in 1:1000) {
  dta = truth[,((sim - 1) * 3 + 1):(sim * 3)]
  ccf_gt = ccf_gt + (ccf_summary(dta) / 1000)
}

for (sim in 1:1000) {
  dta = truth[,((sim - 1) * 3 + 1):(sim * 3)]
  
  lin = benchmarks[,((sim - 1) * 3 + 1):(sim * 3)]
  spl = benchmarks[,3000 + ((sim - 1) * 3 + 1):(sim * 3)]
  Kalman = benchmarks[,6000 + ((sim - 1) * 3 + 1):(sim * 3)]
  PT = benchmarks[,9000 + ((sim - 1) * 3 + 1):(sim * 3)]
  wass = WI_lin[,((sim - 1) * 3 + 1):(sim * 3)]
  kwass = WI_lin[,3000 + ((sim - 1) * 3 + 1):(sim * 3)]
  wass_Kalman = WI_Kalman[,((sim - 1) * 3 + 1):(sim * 3)]
  kwass_Kalman = WI_Kalman[,3000 + ((sim - 1) * 3 + 1):(sim * 3)]
  
  wass_d[sim, 1] = wasserstein(pp(embed(dta, 3)), pp(embed(lin, 3)), p = 2)
  wass_d[sim, 2] = wasserstein(pp(embed(dta, 3)), pp(embed(spl, 3)), p = 2)
  wass_d[sim, 3] = wasserstein(pp(embed(dta, 3)), pp(embed(Kalman, 3)), p = 2)
  wass_d[sim, 4] = wasserstein(pp(embed(dta, 3)), pp(embed(PT, 3)), p = 2)
  wass_d[sim, 5] = wasserstein(pp(embed(dta, 3)), pp(embed(wass, 3)), p = 2)
  wass_d[sim, 6] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass, 3)), p = 2)
  wass_d[sim, 7] = wasserstein(pp(embed(dta, 3)), pp(embed(wass_Kalman, 3)), p = 2)
  wass_d[sim, 8] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass_Kalman, 3)), p = 2)
  
  ccfs[,,,1] = ccfs[,,,1] + ccf_summary(lin) / 1000
  ccfs[,,,2] = ccfs[,,,2] + ccf_summary(spl) / 1000
  ccfs[,,,3] = ccfs[,,,3] + ccf_summary(Kalman) / 1000
  ccfs[,,,4] = ccfs[,,,4] + ccf_summary(PT) / 1000
  ccfs[,,,5] = ccfs[,,,5] + ccf_summary(wass) / 1000
  ccfs[,,,6] = ccfs[,,,6] + ccf_summary(kwass) / 1000
  ccfs[,,,7] = ccfs[,,,7] + ccf_summary(wass_Kalman) / 1000
  ccfs[,,,8] = ccfs[,,,8] + ccf_summary(kwass_Kalman) / 1000
  
  ccf_loss[,1] = ccf_loss[,1] + apply((ccf_summary(lin) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,2] = ccf_loss[,2] + apply((ccf_summary(spl) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,3] = ccf_loss[,3] + apply((ccf_summary(Kalman) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,4] = ccf_loss[,4] + apply((ccf_summary(PT) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,5] = ccf_loss[,5] + apply((ccf_summary(wass) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,6] = ccf_loss[,6] + apply((ccf_summary(kwass) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,7] = ccf_loss[,7] + apply((ccf_summary(wass_Kalman) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,8] = ccf_loss[,8] + apply((ccf_summary(kwass_Kalman) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  
  if (sim == 1) {
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = lin[which(is.na(x_obs[,1])),1],
                      y = lin[which(is.na(x_obs[,2])),2])
    p1 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], lin[,1])), 
                      ylim = range(c(dta[,2], lin[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs linear",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = lin[which(is.na(x_obs[,1])),1],
                      y = lin[which(is.na(x_obs[,3])),3])
    p2 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], lin[,1])), 
                      ylim = range(c(dta[,3], lin[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs linear",
           x = "Component 1",
           y = "Component 3")

    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = lin[which(is.na(x_obs[,2])),2],
                      y = lin[which(is.na(x_obs[,3])),3])
    p3 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], lin[,2])), 
                      ylim = range(c(dta[,3], lin[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs linear",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = spl[which(is.na(x_obs[,1])),1],
                      y = spl[which(is.na(x_obs[,2])),2])
    p4 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], lin[,1])), 
                      ylim = range(c(dta[,2], lin[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs spline",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = spl[which(is.na(x_obs[,1])),1],
                      y = spl[which(is.na(x_obs[,3])),3])
    p5 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], spl[,1])), 
                      ylim = range(c(dta[,3], spl[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs spline",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = spl[which(is.na(x_obs[,2])),2],
                      y = spl[which(is.na(x_obs[,3])),3])
    p6 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], spl[,2])), 
                      ylim = range(c(dta[,3], spl[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs spline",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = Kalman[which(is.na(x_obs[,1])),1],
                      y = Kalman[which(is.na(x_obs[,2])),2])
    p7 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], Kalman[,1])), 
                      ylim = range(c(dta[,2], Kalman[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs Kalman",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = Kalman[which(is.na(x_obs[,1])),1],
                      y = Kalman[which(is.na(x_obs[,3])),3])
    p8 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], Kalman[,1])), 
                      ylim = range(c(dta[,3], Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs Kalman",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = Kalman[which(is.na(x_obs[,2])),2],
                      y = Kalman[which(is.na(x_obs[,3])),3])
    p9 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], Kalman[,2])), 
                      ylim = range(c(dta[,3], Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs Kalman",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = PT[which(is.na(x_obs[,1])),1],
                      y = PT[which(is.na(x_obs[,2])),2])
    p10 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], PT[,1])), 
                      ylim = range(c(dta[,2], PT[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs ScalarF",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = PT[which(is.na(x_obs[,1])),1],
                      y = PT[which(is.na(x_obs[,3])),3])
    p11 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], PT[,1])), 
                      ylim = range(c(dta[,3], PT[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs ScalarF",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = PT[which(is.na(x_obs[,2])),2],
                      y = PT[which(is.na(x_obs[,3])),3])
    p12 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], PT[,2])), 
                      ylim = range(c(dta[,3], PT[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs ScalarF",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = wass[which(is.na(x_obs[,1])),1],
                      y = wass[which(is.na(x_obs[,2])),2])
    p13 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass[,1])), 
                      ylim = range(c(dta[,2], wass[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (linear)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = wass[which(is.na(x_obs[,1])),1],
                      y = wass[which(is.na(x_obs[,3])),3])
    p14 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass[,1])), 
                      ylim = range(c(dta[,3], wass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (linear)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = wass[which(is.na(x_obs[,2])),2],
                      y = wass[which(is.na(x_obs[,3])),3])
    p15 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], wass[,2])), 
                      ylim = range(c(dta[,3], wass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (linear)",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = kwass[which(is.na(x_obs[,1])),1],
                      y = kwass[which(is.na(x_obs[,2])),2])
    p16 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass[,1])), 
                      ylim = range(c(dta[,2], kwass[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (linear)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = kwass[which(is.na(x_obs[,1])),1],
                      y = kwass[which(is.na(x_obs[,3])),3])
    p17 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass[,1])), 
                      ylim = range(c(dta[,3], kwass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (linear)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = kwass[which(is.na(x_obs[,2])),2],
                      y = kwass[which(is.na(x_obs[,3])),3])
    p18 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], kwass[,2])), 
                      ylim = range(c(dta[,3], kwass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (linear)",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = wass_Kalman[which(is.na(x_obs[,1])),1],
                      y = wass_Kalman[which(is.na(x_obs[,2])),2])
    p19 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass_Kalman[,1])), 
                      ylim = range(c(dta[,2], wass_Kalman[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (Kalman)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = wass_Kalman[which(is.na(x_obs[,1])),1],
                      y = wass_Kalman[which(is.na(x_obs[,3])),3])
    p20 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass_Kalman[,1])), 
                      ylim = range(c(dta[,3], wass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (Kalman)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = wass_Kalman[which(is.na(x_obs[,2])),2],
                      y = wass_Kalman[which(is.na(x_obs[,3])),3])
    p21 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], wass_Kalman[,2])), 
                      ylim = range(c(dta[,3], wass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (Kalman)",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = kwass_Kalman[which(is.na(x_obs[,1])),1],
                      y = kwass_Kalman[which(is.na(x_obs[,2])),2])
    p22 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass_Kalman[,1])), 
                      ylim = range(c(dta[,2], kwass_Kalman[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (Kalman)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = kwass_Kalman[which(is.na(x_obs[,1])),1],
                      y = kwass_Kalman[which(is.na(x_obs[,3])),3])
    p23 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass_Kalman[,1])), 
                      ylim = range(c(dta[,3], kwass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (Kalman)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = kwass_Kalman[which(is.na(x_obs[,2])),2],
                      y = kwass_Kalman[which(is.na(x_obs[,3])),3])
    p24 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], kwass_Kalman[,2])), 
                      ylim = range(c(dta[,3], kwass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (Kalman)",
           x = "Component 2",
           y = "Component 3")
    
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                 nrow = 4, ncol = 3)
    grid.arrange(p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24,
                 nrow = 4, ncol = 3)
  }
  cat("iteration", sim, "\n")
}

temp = as.data.frame(as.table(ccf_gt[,,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p1 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "")

temp = as.data.frame(as.table(ccf_gt[,,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p2 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccf_gt[,,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p3 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p4 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Linear (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p5 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Linear (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p6 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Linear (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p7 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Spline (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p8 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Spline (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p9 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Spline (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,3]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p10 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Kalman (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,3]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p11 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Kalman (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,3]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p12 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Kalman (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p13 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("ScalarF (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p14 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("ScalarF (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p15 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("ScalarF (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15,
             nrow = 5, ncol = 3)

temp = as.data.frame(as.table(ccf_gt[,,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p1 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "")

temp = as.data.frame(as.table(ccf_gt[,,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p2 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccf_gt[,,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p3 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,5]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p4 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-lin (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,5]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p5 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-lin (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,5]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p6 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-lin (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,6]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p7 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-lin (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,6]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p8 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-lin (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,6]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p9 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-lin (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,7]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p10 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-Kalman (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,2,7]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p11 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-Kalman (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,7]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p12 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-Kalman (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,8]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p13 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-Kalman (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,2,8]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p14 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-Kalman (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,8]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p15 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-Kalman (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15,
             nrow = 5, ncol = 3)


print(round(colMeans(wass_d), 4))
print(round(sqrt(ccf_loss / 1000), 4))


################################################################################
################################################################################
# ALR: II

model_name = "ALR"
base_dir = paste0("./sim_data/", model_name, "/")

truth = as.matrix(read.csv(paste0(base_dir, model_name, "_original.csv"), header = F))
x_obs = as.matrix(read.csv(paste0(base_dir, model_name, "_miss2.csv"), header = F))
benchmarks = as.matrix(read.csv(paste0(base_dir, "benchmarks_miss2.csv")))
WI_lin = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_lin_miss2.csv"), header = F))
WI_Kalman = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_Kalman_miss2.csv"), header = F))

wass_d = matrix(0, nrow = 1000, ncol = 8)
ccf_gt = array(0, dim = c(3, 3, 5))
ccfs = array(0, dim = c(3, 3, 5, 8))
ccf_loss = matrix(0, nrow = 5, ncol = 8)

for (sim in 1:1000) {
  dta = truth[,((sim - 1) * 3 + 1):(sim * 3)]
  ccf_gt = ccf_gt + (ccf_summary(dta) / 1000)
}

for (sim in 1:1000) {
  dta = truth[,((sim - 1) * 3 + 1):(sim * 3)]
  
  lin = benchmarks[,((sim - 1) * 3 + 1):(sim * 3)]
  spl = benchmarks[,3000 + ((sim - 1) * 3 + 1):(sim * 3)]
  Kalman = benchmarks[,6000 + ((sim - 1) * 3 + 1):(sim * 3)]
  PT = benchmarks[,9000 + ((sim - 1) * 3 + 1):(sim * 3)]
  wass = WI_lin[,((sim - 1) * 3 + 1):(sim * 3)]
  kwass = WI_lin[,3000 + ((sim - 1) * 3 + 1):(sim * 3)]
  wass_Kalman = WI_Kalman[,((sim - 1) * 3 + 1):(sim * 3)]
  kwass_Kalman = WI_Kalman[,3000 + ((sim - 1) * 3 + 1):(sim * 3)]
  
  wass_d[sim, 1] = wasserstein(pp(embed(dta, 3)), pp(embed(lin, 3)), p = 2)
  wass_d[sim, 2] = wasserstein(pp(embed(dta, 3)), pp(embed(spl, 3)), p = 2)
  wass_d[sim, 3] = wasserstein(pp(embed(dta, 3)), pp(embed(Kalman, 3)), p = 2)
  wass_d[sim, 4] = wasserstein(pp(embed(dta, 3)), pp(embed(PT, 3)), p = 2)
  wass_d[sim, 5] = wasserstein(pp(embed(dta, 3)), pp(embed(wass, 3)), p = 2)
  wass_d[sim, 6] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass, 3)), p = 2)
  wass_d[sim, 7] = wasserstein(pp(embed(dta, 3)), pp(embed(wass_Kalman, 3)), p = 2)
  wass_d[sim, 8] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass_Kalman, 3)), p = 2)
  
  ccfs[,,,1] = ccfs[,,,1] + ccf_summary(lin) / 1000
  ccfs[,,,2] = ccfs[,,,2] + ccf_summary(spl) / 1000
  ccfs[,,,3] = ccfs[,,,3] + ccf_summary(Kalman) / 1000
  ccfs[,,,4] = ccfs[,,,4] + ccf_summary(PT) / 1000
  ccfs[,,,5] = ccfs[,,,5] + ccf_summary(wass) / 1000
  ccfs[,,,6] = ccfs[,,,6] + ccf_summary(kwass) / 1000
  ccfs[,,,7] = ccfs[,,,7] + ccf_summary(wass_Kalman) / 1000
  ccfs[,,,8] = ccfs[,,,8] + ccf_summary(kwass_Kalman) / 1000
  
  ccf_loss[,1] = ccf_loss[,1] + apply((ccf_summary(lin) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,2] = ccf_loss[,2] + apply((ccf_summary(spl) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,3] = ccf_loss[,3] + apply((ccf_summary(Kalman) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,4] = ccf_loss[,4] + apply((ccf_summary(PT) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,5] = ccf_loss[,5] + apply((ccf_summary(wass) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,6] = ccf_loss[,6] + apply((ccf_summary(kwass) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,7] = ccf_loss[,7] + apply((ccf_summary(wass_Kalman) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  ccf_loss[,8] = ccf_loss[,8] + apply((ccf_summary(kwass_Kalman) - ccf_gt)^2, MARGIN = 3, FUN = sum)
  
  if (sim == 1) {
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = lin[which(is.na(x_obs[,1])),1],
                      y = lin[which(is.na(x_obs[,2])),2])
    p1 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], lin[,1])), 
                      ylim = range(c(dta[,2], lin[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs linear",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = lin[which(is.na(x_obs[,1])),1],
                      y = lin[which(is.na(x_obs[,3])),3])
    p2 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], lin[,1])), 
                      ylim = range(c(dta[,3], lin[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs linear",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = lin[which(is.na(x_obs[,2])),2],
                      y = lin[which(is.na(x_obs[,3])),3])
    p3 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], lin[,2])), 
                      ylim = range(c(dta[,3], lin[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs linear",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = spl[which(is.na(x_obs[,1])),1],
                      y = spl[which(is.na(x_obs[,2])),2])
    p4 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], lin[,1])), 
                      ylim = range(c(dta[,2], lin[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs spline",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = spl[which(is.na(x_obs[,1])),1],
                      y = spl[which(is.na(x_obs[,3])),3])
    p5 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], spl[,1])), 
                      ylim = range(c(dta[,3], spl[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs spline",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = spl[which(is.na(x_obs[,2])),2],
                      y = spl[which(is.na(x_obs[,3])),3])
    p6 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], spl[,2])), 
                      ylim = range(c(dta[,3], spl[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs spline",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = Kalman[which(is.na(x_obs[,1])),1],
                      y = Kalman[which(is.na(x_obs[,2])),2])
    p7 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], Kalman[,1])), 
                      ylim = range(c(dta[,2], Kalman[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs Kalman",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = Kalman[which(is.na(x_obs[,1])),1],
                      y = Kalman[which(is.na(x_obs[,3])),3])
    p8 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], Kalman[,1])), 
                      ylim = range(c(dta[,3], Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs Kalman",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = Kalman[which(is.na(x_obs[,2])),2],
                      y = Kalman[which(is.na(x_obs[,3])),3])
    p9 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], Kalman[,2])), 
                      ylim = range(c(dta[,3], Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs Kalman",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = PT[which(is.na(x_obs[,1])),1],
                      y = PT[which(is.na(x_obs[,2])),2])
    p10 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], PT[,1])), 
                      ylim = range(c(dta[,2], PT[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs ScalarF",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = PT[which(is.na(x_obs[,1])),1],
                      y = PT[which(is.na(x_obs[,3])),3])
    p11 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], PT[,1])), 
                      ylim = range(c(dta[,3], PT[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs ScalarF",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = PT[which(is.na(x_obs[,2])),2],
                      y = PT[which(is.na(x_obs[,3])),3])
    p12 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], PT[,2])), 
                      ylim = range(c(dta[,3], PT[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs ScalarF",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = wass[which(is.na(x_obs[,1])),1],
                      y = wass[which(is.na(x_obs[,2])),2])
    p13 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass[,1])), 
                      ylim = range(c(dta[,2], wass[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (linear)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = wass[which(is.na(x_obs[,1])),1],
                      y = wass[which(is.na(x_obs[,3])),3])
    p14 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass[,1])), 
                      ylim = range(c(dta[,3], wass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (linear)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = wass[which(is.na(x_obs[,2])),2],
                      y = wass[which(is.na(x_obs[,3])),3])
    p15 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], wass[,2])), 
                      ylim = range(c(dta[,3], wass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (linear)",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = kwass[which(is.na(x_obs[,1])),1],
                      y = kwass[which(is.na(x_obs[,2])),2])
    p16 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass[,1])), 
                      ylim = range(c(dta[,2], kwass[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (linear)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = kwass[which(is.na(x_obs[,1])),1],
                      y = kwass[which(is.na(x_obs[,3])),3])
    p17 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass[,1])), 
                      ylim = range(c(dta[,3], kwass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (linear)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = kwass[which(is.na(x_obs[,2])),2],
                      y = kwass[which(is.na(x_obs[,3])),3])
    p18 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], kwass[,2])), 
                      ylim = range(c(dta[,3], kwass[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (linear)",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = wass_Kalman[which(is.na(x_obs[,1])),1],
                      y = wass_Kalman[which(is.na(x_obs[,2])),2])
    p19 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass_Kalman[,1])), 
                      ylim = range(c(dta[,2], wass_Kalman[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (Kalman)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = wass_Kalman[which(is.na(x_obs[,1])),1],
                      y = wass_Kalman[which(is.na(x_obs[,3])),3])
    p20 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], wass_Kalman[,1])), 
                      ylim = range(c(dta[,3], wass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (Kalman)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = wass_Kalman[which(is.na(x_obs[,2])),2],
                      y = wass_Kalman[which(is.na(x_obs[,3])),3])
    p21 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], wass_Kalman[,2])), 
                      ylim = range(c(dta[,3], wass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs WI (Kalman)",
           x = "Component 2",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,1], y = dta[,2])
    temp = data.frame(x = kwass_Kalman[which(is.na(x_obs[,1])),1],
                      y = kwass_Kalman[which(is.na(x_obs[,2])),2])
    p22 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass_Kalman[,1])), 
                      ylim = range(c(dta[,2], kwass_Kalman[,2]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (Kalman)",
           x = "Component 1",
           y = "Component 2")
    
    df1 = data.frame(x = dta[,1], y = dta[,3])
    temp = data.frame(x = kwass_Kalman[which(is.na(x_obs[,1])),1],
                      y = kwass_Kalman[which(is.na(x_obs[,3])),3])
    p23 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,1], kwass_Kalman[,1])), 
                      ylim = range(c(dta[,3], kwass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (Kalman)",
           x = "Component 1",
           y = "Component 3")
    
    df1 = data.frame(x = dta[,2], y = dta[,3])
    temp = data.frame(x = kwass_Kalman[which(is.na(x_obs[,2])),2],
                      y = kwass_Kalman[which(is.na(x_obs[,3])),3])
    p24 = ggplot() +
      geom_point(data = df1, aes(x = x, y = y), color = "steelblue", size = 1) + 
      geom_point(data = temp, aes(x = x, y = y), color = "red3", size = 1, shape = 17) +
      coord_cartesian(xlim = range(c(dta[,2], kwass_Kalman[,2])), 
                      ylim = range(c(dta[,3], kwass_Kalman[,3]))) +
      theme_minimal() +
      labs(title = "Ground truth vs kWI (Kalman)",
           x = "Component 2",
           y = "Component 3")
    
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                 nrow = 4, ncol = 3)
    grid.arrange(p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24,
                 nrow = 4, ncol = 3)
    
    p0 = imp_plot(dta, x_obs, "Ground truth")
    p1 = imp_plot(lin, x_obs, "Linear")
    p2 = imp_plot(spl, x_obs, "Spline")
    p3 = imp_plot(Kalman, x_obs, "Kalman")
    p4 = imp_plot(PT, x_obs, "ScalarF")
    p5 = imp_plot(wass, x_obs, "WI (linear)")
    p6 = imp_plot(kwass, x_obs, "kWI (linear)")
    p7 = imp_plot(wass_Kalman, x_obs, "WI (Kalman)")
    p8 = imp_plot(kwass_Kalman, x_obs, "kWI (Kalman)")
    
    grid.arrange(p0, p1, p2, p3, p4, p5, p6, p7, p8, nrow=3, ncol=3)
  }
  cat("iteration", sim, "\n")
}

temp = as.data.frame(as.table(ccf_gt[,,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p1 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "")

temp = as.data.frame(as.table(ccf_gt[,,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p2 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccf_gt[,,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p3 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p4 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Linear (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p5 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Linear (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p6 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Linear (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p7 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Spline (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p8 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Spline (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p9 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Spline (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,3]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p10 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Kalman (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,3]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p11 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Kalman (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,3]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p12 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Kalman (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p13 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("ScalarF (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p14 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("ScalarF (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p15 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("ScalarF (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15,
             nrow = 5, ncol = 3)

temp = as.data.frame(as.table(ccf_gt[,,1]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p1 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "")

temp = as.data.frame(as.table(ccf_gt[,,2]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p2 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccf_gt[,,4]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p3 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("Ground truth (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,5]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p4 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-lin (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,5]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p5 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-lin (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,5]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p6 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-lin (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,6]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p7 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-lin (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")


temp = as.data.frame(as.table(ccfs[,,2,6]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p8 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-lin (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,6]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p9 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-lin (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,7]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p10 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-Kalman (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,2,7]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p11 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-Kalman (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,7]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p12 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("WI-Kalman (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,1,8]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p13 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-Kalman (cross-correlation Lag ", 0, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,2,8]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p14 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-Kalman (cross-correlation Lag ", 1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

temp = as.data.frame(as.table(ccfs[,,4,8]))
temp$Var1 <- factor(temp$Var1, levels = rev(unique(temp$Var1)), labels = c("x3", "x2", "x1"))
temp$Var2 <- factor(temp$Var2, levels = unique(temp$Var2), labels = c("x1", "x2", "x3"))

p15 = ggplot(temp, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "steelblue", limits = c(-1, 1)) +
  theme_minimal()  +
  labs(title = paste0("kWI-Kalman (cross-correlation Lag ", -1, ")"),
       x = "",
       y = "") +
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15,
             nrow = 5, ncol = 3)


print(round(colMeans(wass_d), 4))
print(round(sqrt(ccf_loss / 1000), 4))








