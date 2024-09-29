library(transport)
library(ggplot2)
library(gridExtra)

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

model_estimate = function (x, model_name) {
  
  if (model_name == "NLVAR") {
    dta = embed(x, 2)
    y1 = dta[,1]
    y2 = dta[,2]
    x1 = dta[,3]
    x2 = dta[,4]
    x2_trans = (1 / (1 + exp(-3 * x2))) - 0.5
    
    model1 = lm(y1~0+x1+x2_trans)
    model2 = lm(y2~0+x1+x2)
  }

  return (unname(c(model1$coefficients, model2$coefficients)))
}

# NLVAR: I

model_name = "NLVAR"
base_dir = paste0("./sim_data/", model_name, "/")

truth = as.matrix(read.csv(paste0(base_dir, model_name, "_original.csv"), header = F))
x_obs = as.matrix(read.csv(paste0(base_dir, model_name, "_miss1.csv"), header = F))
benchmarks = as.matrix(read.csv(paste0(base_dir, "benchmarks_miss1.csv")))
WI_lin = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_lin_miss1.csv"), header = F))
WI_Kalman = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_Kalman_miss1.csv"), header = F))

wass_d = matrix(0, nrow = 1000, ncol = 8)
model_coef = vector(mode = "list", length = 9)

for (sim in 1:1000) {
  dta = truth[,((sim - 1) * 2 + 1):(sim * 2)]
  
  lin = benchmarks[,((sim - 1) * 2 + 1):(sim * 2)]
  spl = benchmarks[,2000 + ((sim - 1) * 2 + 1):(sim * 2)]
  Kalman = benchmarks[,4000 + ((sim - 1) * 2 + 1):(sim * 2)]
  PT = benchmarks[,6000 + ((sim - 1) * 2 + 1):(sim * 2)]
  wass = WI_lin[,((sim - 1) * 2 + 1):(sim * 2)]
  kwass = WI_lin[,2000 + ((sim - 1) * 2 + 1):(sim * 2)]
  wass_Kalman = WI_Kalman[,((sim - 1) * 2 + 1):(sim * 2)]
  kwass_Kalman = WI_Kalman[,2000 + ((sim - 1) * 2 + 1):(sim * 2)]
  
  wass_d[sim, 1] = wasserstein(pp(embed(dta, 3)), pp(embed(lin, 3)), p = 2)
  wass_d[sim, 2] = wasserstein(pp(embed(dta, 3)), pp(embed(spl, 3)), p = 2)
  wass_d[sim, 3] = wasserstein(pp(embed(dta, 3)), pp(embed(Kalman, 3)), p = 2)
  wass_d[sim, 4] = wasserstein(pp(embed(dta, 3)), pp(embed(PT, 3)), p = 2)
  wass_d[sim, 5] = wasserstein(pp(embed(dta, 3)), pp(embed(wass, 3)), p = 2)
  wass_d[sim, 6] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass, 3)), p = 2)
  wass_d[sim, 7] = wasserstein(pp(embed(dta, 3)), pp(embed(wass_Kalman, 3)), p = 2)
  wass_d[sim, 8] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass_Kalman, 3)), p = 2)
  
  if (sim == 1) {
    model_coef[[1]] = model_estimate(dta, model_name)
    model_coef[[2]] = model_estimate(lin, model_name)
    model_coef[[3]] = model_estimate(spl, model_name)
    model_coef[[4]] = model_estimate(Kalman, model_name)
    model_coef[[5]] = model_estimate(PT, model_name)
    model_coef[[6]] = model_estimate(wass, model_name)
    model_coef[[7]] = model_estimate(kwass, model_name)
    model_coef[[8]] = model_estimate(wass_Kalman, model_name)
    model_coef[[9]] = model_estimate(kwass_Kalman, model_name)
  } else {
    model_coef[[1]] = rbind(model_coef[[1]], model_estimate(dta, model_name))
    model_coef[[2]] = rbind(model_coef[[2]], model_estimate(lin, model_name))
    model_coef[[3]] = rbind(model_coef[[3]], model_estimate(spl, model_name))
    model_coef[[4]] = rbind(model_coef[[4]], model_estimate(Kalman, model_name))
    model_coef[[5]] = rbind(model_coef[[5]], model_estimate(PT, model_name))
    model_coef[[6]] = rbind(model_coef[[6]], model_estimate(wass, model_name))
    model_coef[[7]] = rbind(model_coef[[7]], model_estimate(kwass, model_name))
    model_coef[[8]] = rbind(model_coef[[8]], model_estimate(wass_Kalman, model_name))
    model_coef[[9]] = rbind(model_coef[[9]], model_estimate(kwass_Kalman, model_name))
  }
  
  if (sim == 1) {
    plot_range = range(dta)
    plot_range = c(floor(plot_range[1]), ceiling(plot_range[2]))
    
    x_coord = 1
    y_coord = 1
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
                                   y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
    x_coord = 1
    y_coord = 2
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
    x_coord = 2
    y_coord = 1
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
    x_coord = 2
    y_coord = 2
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
  }
  cat("iteration", sim, "\n")
}

par(mfrow = c(3, 3))
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,1], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(0.2, 0.65))
    lines(density(model_coef[[1]][,1]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,1], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi11"),
         xlab = "",
         xlim = c(0.2, 0.65))
    lines(density(model_coef[[1]][,1]), col = "blue", lwd = 2) 
  }
}
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,2], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(4.2, 8.1))
    lines(density(model_coef[[1]][,2]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,2], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi12"),
         xlab = "", 
         xlim = c(4.2, 8.1))
    lines(density(model_coef[[1]][,2]), col = "blue", lwd = 2)
  }
}
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,3], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(-0.1, 0.1))
    lines(density(model_coef[[1]][,3]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,3], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi21"),
         xlab = "",
         xlim = c(-0.1, 0.1))
    lines(density(model_coef[[1]][,3]), col = "blue", lwd = 2)
  }
}
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,4], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(0.2, 0.75))
    lines(density(model_coef[[1]][,4]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,4], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi22"),
         xlab = "",
         xlim = c(0.2, 0.75))
    lines(density(model_coef[[1]][,4]), col = "blue", lwd = 2)
  }
}

print(round(colMeans(wass_d), 4))
GT_arcoef = colMeans(model_coef[[1]])
for (j in 1:8) {
  sub_arcoef = model_coef[[j+1]]
  temp = t(t(sub_arcoef) - GT_arcoef)
  if (j == 1) {
    cat("estimation error for", "Linear:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 2) {
    cat("estimation error for", "Spline:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 3) {
    cat("estimation error for", "Kalman:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 4) {
    cat("estimation error for", "PT:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 5) {
    cat("estimation error for", "WI (linear):", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 6) {
    cat("estimation error for", "kWI (linear):", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 7) {
    cat("estimation error for", "WI (Kalman):", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else {
    cat("estimation error for", "kWI (Kalman):", round(sqrt(colMeans(temp^2)), 3), "\n")
  }
}

################################################################################
# NLVAR: II

model_name = "NLVAR"
base_dir = paste0("./sim_data/", model_name, "/")

truth = as.matrix(read.csv(paste0(base_dir, model_name, "_original.csv"), header = F))
x_obs = as.matrix(read.csv(paste0(base_dir, model_name, "_miss2.csv"), header = F))
benchmarks = as.matrix(read.csv(paste0(base_dir, "benchmarks_miss2.csv")))
WI_lin = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_lin_miss2.csv"), header = F))
WI_Kalman = as.matrix(read.csv(paste0(base_dir, model_name, "_WI_Kalman_miss2.csv"), header = F))

wass_d = matrix(0, nrow = 1000, ncol = 8)
model_coef = vector(mode = "list", length = 9)

for (sim in 1:1000) {
  dta = truth[,((sim - 1) * 2 + 1):(sim * 2)]
  
  lin = benchmarks[,((sim - 1) * 2 + 1):(sim * 2)]
  spl = benchmarks[,2000 + ((sim - 1) * 2 + 1):(sim * 2)]
  Kalman = benchmarks[,4000 + ((sim - 1) * 2 + 1):(sim * 2)]
  PT = benchmarks[,6000 + ((sim - 1) * 2 + 1):(sim * 2)]
  wass = WI_lin[,((sim - 1) * 2 + 1):(sim * 2)]
  kwass = WI_lin[,2000 + ((sim - 1) * 2 + 1):(sim * 2)]
  wass_Kalman = WI_Kalman[,((sim - 1) * 2 + 1):(sim * 2)]
  kwass_Kalman = WI_Kalman[,2000 + ((sim - 1) * 2 + 1):(sim * 2)]
  
  wass_d[sim, 1] = wasserstein(pp(embed(dta, 3)), pp(embed(lin, 3)), p = 2)
  wass_d[sim, 2] = wasserstein(pp(embed(dta, 3)), pp(embed(spl, 3)), p = 2)
  wass_d[sim, 3] = wasserstein(pp(embed(dta, 3)), pp(embed(Kalman, 3)), p = 2)
  wass_d[sim, 4] = wasserstein(pp(embed(dta, 3)), pp(embed(PT, 3)), p = 2)
  wass_d[sim, 5] = wasserstein(pp(embed(dta, 3)), pp(embed(wass, 3)), p = 2)
  wass_d[sim, 6] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass, 3)), p = 2)
  wass_d[sim, 7] = wasserstein(pp(embed(dta, 3)), pp(embed(wass_Kalman, 3)), p = 2)
  wass_d[sim, 8] = wasserstein(pp(embed(dta, 3)), pp(embed(kwass_Kalman, 3)), p = 2)
  
  if (sim == 1) {
    model_coef[[1]] = model_estimate(dta, model_name)
    model_coef[[2]] = model_estimate(lin, model_name)
    model_coef[[3]] = model_estimate(spl, model_name)
    model_coef[[4]] = model_estimate(Kalman, model_name)
    model_coef[[5]] = model_estimate(PT, model_name)
    model_coef[[6]] = model_estimate(wass, model_name)
    model_coef[[7]] = model_estimate(kwass, model_name)
    model_coef[[8]] = model_estimate(wass_Kalman, model_name)
    model_coef[[9]] = model_estimate(kwass_Kalman, model_name)
  } else {
    model_coef[[1]] = rbind(model_coef[[1]], model_estimate(dta, model_name))
    model_coef[[2]] = rbind(model_coef[[2]], model_estimate(lin, model_name))
    model_coef[[3]] = rbind(model_coef[[3]], model_estimate(spl, model_name))
    model_coef[[4]] = rbind(model_coef[[4]], model_estimate(Kalman, model_name))
    model_coef[[5]] = rbind(model_coef[[5]], model_estimate(PT, model_name))
    model_coef[[6]] = rbind(model_coef[[6]], model_estimate(wass, model_name))
    model_coef[[7]] = rbind(model_coef[[7]], model_estimate(kwass, model_name))
    model_coef[[8]] = rbind(model_coef[[8]], model_estimate(wass_Kalman, model_name))
    model_coef[[9]] = rbind(model_coef[[9]], model_estimate(kwass_Kalman, model_name))
  }
  
  if (sim == 1) {
    plot_range = range(dta)
    plot_range = c(floor(plot_range[1]), ceiling(plot_range[2]))
    
    x_coord = 1
    y_coord = 1
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
    x_coord = 1
    y_coord = 2
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
    x_coord = 2
    y_coord = 1
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
    x_coord = 2
    y_coord = 2
    
    df1 = data.frame(x = embed(dta[,x_coord], 2)[,2], y = embed(dta[,y_coord], 2)[,1])
    p1 = ggplot(df1, aes(x = x, y = y)) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      geom_point(color = "steelblue", size = 1) +
      theme_minimal() +
      labs(title = "Ground truth", x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(lin[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(lin[,y_coord], which(is.na(x_obs[,1])))[,2])
    p2 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Linear interpolation", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(spl[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(spl[,y_coord], which(is.na(x_obs[,1])))[,2])
    p3 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Spline smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p4 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Kalman smoothing", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(PT[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(PT[,y_coord], which(is.na(x_obs[,1])))[,2])
    p5 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "Scalar filter", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p6 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass[,y_coord], which(is.na(x_obs[,1])))[,2])
    p7 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (linear)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(wass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(wass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p8 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "WI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    temp = data.frame(x = na_lag_pairs(kwass_Kalman[,x_coord], which(is.na(x_obs[,1])))[,1],
                      y = na_lag_pairs(kwass_Kalman[,y_coord], which(is.na(x_obs[,1])))[,2])
    p9 = ggplot(temp, aes(x = x, y = y)) + 
      geom_point(color = "steelblue1", size = 1) +
      coord_cartesian(xlim = plot_range, ylim = plot_range) +
      theme_minimal() +
      labs(title = "kWI (Kalman)", 
           x = paste0("x_{t-1,", x_coord, "}"), 
           y = paste0("x_{t,", y_coord,"}"))
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
    
  }
  cat("iteration", sim, "\n")
}

par(mfrow = c(3, 3))
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,1], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(0.2, 0.65))
    lines(density(model_coef[[1]][,1]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,1], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi11"),
         xlab = "",
         xlim = c(0.2, 0.65))
    lines(density(model_coef[[1]][,1]), col = "blue", lwd = 2) 
  }
}
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,2], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(4.2, 8.1))
    lines(density(model_coef[[1]][,2]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,2], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi12"),
         xlab = "", 
         xlim = c(4.2, 8.1))
    lines(density(model_coef[[1]][,2]), col = "blue", lwd = 2)
  }
}
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,3], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(-0.1, 0.1))
    lines(density(model_coef[[1]][,3]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,3], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi21"),
         xlab = "",
         xlim = c(-0.1, 0.1))
    lines(density(model_coef[[1]][,3]), col = "blue", lwd = 2)
  }
}
for (i in 1:9) {
  if (i == 1) {
    hist(model_coef[[1]][,4], prob = TRUE, col = "lightblue", 
         main = "Ground truth", xlab = "", xlim = c(0.2, 0.75))
    lines(density(model_coef[[1]][,4]), col = "blue", lwd = 2)
  } else {
    if (i == 2) {
      method = "Linear"
    } else if (i == 3) {
      method = "Spline"
    } else if (i == 4) {
      method = "Kalman"
    } else if (i == 5) {
      method = "Scalar filter"
    } else if (i == 6) {
      method = "WI (linear)"
    } else if (i == 7) {
      method = "kWI (linear)"
    } else if (i == 8) {
      method = "WI (Kalman)"
    } else {
      method = "kWI (Kalman)"
    }
    hist(model_coef[[i]][,4], prob = TRUE, col = "lightblue", 
         main = paste0(method, ": phi22"),
         xlab = "",
         xlim = c(0.2, 0.75))
    lines(density(model_coef[[1]][,4]), col = "blue", lwd = 2)
  }
}

print(round(colMeans(wass_d), 4))
GT_arcoef = colMeans(model_coef[[1]])
for (j in 1:8) {
  sub_arcoef = model_coef[[j+1]]
  temp = t(t(sub_arcoef) - GT_arcoef)
  if (j == 1) {
    cat("estimation error for", "Linear:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 2) {
    cat("estimation error for", "Spline:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 3) {
    cat("estimation error for", "Kalman:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 4) {
    cat("estimation error for", "PT:", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 5) {
    cat("estimation error for", "WI (linear):", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 6) {
    cat("estimation error for", "kWI (linear):", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else if (j == 7) {
    cat("estimation error for", "WI (Kalman):", round(sqrt(colMeans(temp^2)), 3), "\n")
  } else {
    cat("estimation error for", "kWI (Kalman):", round(sqrt(colMeans(temp^2)), 3), "\n")
  }
}

