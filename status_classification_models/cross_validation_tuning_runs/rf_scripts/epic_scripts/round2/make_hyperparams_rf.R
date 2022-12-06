#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

# Random Forest hyperparams

# create a colon-separated file of hyperparam values for launching SLURM jobs

# mtry (20 levels)
# min_n (10 levels)
#------------------------------------------------------------------------------#

library(tidymodels)

# number of levels
n_level1 <- 20
n_level2 <- 10

# create values of param1
vals1 <- as.integer(seq(1, 100, length=n_level1)) # mtry with 100 nCpG

# create values of param2
vals2 <- as.integer(seq(1, 15, length=n_level2)) # min_n

# create the grid of hyperparameters
param_vals <- expand.grid(vals1, vals2)

# create a seed column
param_vals$seed <- c(1:nrow(param_vals))

# export
write.table(param_vals,
    "rf_hyperparam_args_file",
    sep = ":",
    col.names = F,
    row.names = F,
    quote = F
)
