#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

# Random Forest hyperparams

# create a colon-separated file of hyperparam values for launching SLURM jobs

# ntree (10 levels)
# nCpG (10 levels)
#------------------------------------------------------------------------------#

library(tidymodels)

# number of levels
n_level1 <- 10
n_level2 <- 10

# create values of param1
vals1 <- as.integer(seq(10, 1000, length=n_level1)) # nCpG

# create values of param2
vals2 <- as.integer(seq(50, 1000, length=n_level2)) # ntree

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
