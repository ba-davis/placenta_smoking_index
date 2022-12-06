#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

# LASSO hyperparams

# create a colon-separated file of hyperparam values for launching SLURM jobs
# penalty values (use grid_regular with 20 values)
# nCpG values (use seq from 10-1000 with 10 values)
# random seed values 

#------------------------------------------------------------------------------#

library(tidymodels)

# number of levels per hyperparameter
n_level1 <- 20 # penalty values
n_level2 <- 10 # nCpG values

# create values for first hyperparameter
vals1 <- as.data.frame(grid_regular(penalty(), levels = n_level1))$penalty

# create values for second hyperparameter
# numeric vector of number of top CpGs to use during variable pre-selection
vals2 <- as.integer(seq(10, 1000, length = n_level2))

# create the grid of hyperparameters
param_vals <- expand.grid(vals1, vals2)
# create a seed column
param_vals$seed <- c(1:nrow(param_vals))

print(paste0("Number of runs: ", nrow(param_vals)))

# export
write.table(param_vals,
    "lasso_hyperparam_args_file",
    sep = ":",
    col.names = F,
    row.names = F,
    quote = F
)