model_names = c("AR", "ARMA", "TAR", "GARCH", "ALR", "ARI", "Cyc", "NLVAR")
ds = c(1, 1, 1, 1, 3, 1, 1, 2)

for (i in 1:length(model_names)) {
  model_name = model_names[i]
  d = ds[i]
  n_exper = 1000
  source("./bchmk_analysis.R")
  rm(list = setdiff(ls(), c("model_names", "ds")))
}