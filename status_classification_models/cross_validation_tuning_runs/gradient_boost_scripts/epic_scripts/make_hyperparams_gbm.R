#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

# GBM Hyperparams

#------------------------------------------------------------------------------#

library(tidymodels)

# number of levels
n_level1 <- 20 # ntree values
n_level2 <- 10 # tree depth values
n_level3 <- 5  # nCpG valuess

vals1 <- (seq(50, 1000, length=n_level1))

vals2 <- as.integer(seq(1, 20, length=n_level2))

vals3 <- as.integer(seq(50, 1000, length=n_level3))

# create the grid of hyperparameters
param_vals <- expand.grid(vals1, vals2, vals3)
# create a seed column
param_vals$seed <- c(1:nrow(param_vals))

# export
write.table(param_vals,
    "gbm_hyperparam_args_file",
    sep = ":",
    col.names = F,
    row.names = F,
    quote = F
)
