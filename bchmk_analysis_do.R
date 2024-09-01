model_names = c("ALR", "NLVAR", "ARI", "Cyc", "AR", "ARMA", "TAR", "GARCH")
ds = c(3, 2, 1, 1, 1, 1, 1, 1)
n_exper = 1000
show_plot = FALSE

for (i in 1:length(model_names)) {
  model_name = model_names[i]
  d = ds[i]
  cat("Starting: ", model_name, "\n")
  source("./bchmk_analysis.R")
  rm(list = setdiff(ls(), c("model_names", "ds", "n_exper", "show_plot")))
}