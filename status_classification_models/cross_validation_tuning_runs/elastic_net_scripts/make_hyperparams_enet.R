#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

# Elastic Net hyperparams

# create a colon-separated file of hyperparam values for launching SLURM jobs
# penalty values (use grid_regular with 10 values)
# mixture values (use grid_regular with 10 values)
# nCpG values (use seq from 50-1000 with 5 values)
# random seed values 

#------------------------------------------------------------------------------#

library(tidymodels)

# number of levels
n_level1 <- 10 # penalty values
n_level2 <- 10 # mixture values
n_level3 <- 5  # nCpG values

# create reasonable penalty values, store in a numeric vector
vals1 <- as.data.frame(grid_regular(penalty(), levels=n_level1))$penalty

# create reasonable mixture values, store in a numeric vector
vals2 <- as.data.frame(grid_regular(mixture(), levels=n_level2))$mixture

# create numeric vector of number of top CpGs to use during variable pre-selection
vals3 <- as.integer(seq(50, 1000, length=n_level3))

# create the grid of hyperparameters
param_vals <- expand.grid(vals1, vals2, vals3)
# create a seed column
param_vals$seed <- c(1:nrow(param_vals))

# export
write.table(param_vals,
    "enet_hyperparam_args_file",
    sep = ":",
    col.names = F,
    row.names = F,
    quote = F
)
