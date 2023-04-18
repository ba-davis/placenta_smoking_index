
# For each cpg in one of the three final models,
# add a mean beta value to be the weighted average

library(tidyverse)

#--------------------------------------------------------#
# FUNCTIONS:

# function to calculate average
calc_avg <- function(df) {
    mean_vec <- vector()
    for (i in seq_len(ncol(df))) {
        mean_val <- mean(df[, i])
        mean_vec <- c(mean_vec, mean_val)
    }

    newdf <- data.frame(cpg = colnames(df),
        mean = mean_vec)

    return(newdf)
}

# calculate weighted average per cpg
# assumes the two input dfs have cpgs in same order
calc_wt_avg <- function(nonsmoker_df, smoker_df, smoke_rate = 0.11) {
    wt_mean <- vector()
    for (i in seq_len(nrow(nonsmoker_df))) {
        smoker_wt_avg <- smoker_df$mean[i] * smoke_rate
        nonsmoker_wt_avg <- nonsmoker_df$mean[i] * (1 - smoke_rate)
        weighted_avg <- smoker_wt_avg + nonsmoker_wt_avg

        wt_mean <- c(wt_mean, weighted_avg)
    }

    newdf <- data.frame(cpg = nonsmoker_df$cpg,
        weighted_mean = wt_mean)

    return(newdf)
}

#--------------------------------------------------------#

# read in metadata table
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
    "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))
table(pd$status)
# our VCSIP data has 24 nonsmokers and 72 smokers
# 25% of VCSIP samples are nonsmokers, 75% are smokers
# national average is 89% nonsmokers, 11% smokers according to cdc
# we will use a rate of 12.3% sourced by Lyndsey

#--------------------------------------------------------#
# Model 1 create reference table

inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model1/")

# read in the raw beta matrix for model 1 cpgs of 96 VCSIP samples
beta <- read.delim(paste0(inpath, "model1.vcsip.raw_beta_matrix.txt"),
    header = TRUE)

# subset beta matrix to nonsmokers only
beta_nonsmokers <- beta[rownames(beta) %in%
    rownames(pd[pd$status == "nonsmoker", ]), ]
# subset beta matrix to smokers only
beta_smokers <- beta[rownames(beta) %in%
    rownames(pd[pd$status == "smoker", ]), ]

# calculate average beta value for smokers and nonsmokers separately
nonsmoker_mean_df <- calc_avg(beta_nonsmokers)
smoker_mean_df <- calc_avg(beta_smokers)

# calculate weighted average
model1_wt_avg <- calc_wt_avg(nonsmoker_mean_df,
    smoker_mean_df,
    smoke_rate = 0.123)

# export
write.table(model1_wt_avg,
    "model1_vcsip_weighted_avg_table.txt",
    sep="\t",
    col.names=T,
    row.names=T,
    quote=F)

#--------------------------------------------------------#
# Model 2 create reference table

inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model2/ncpg120/")

# read in the raw beta matrix for model 1 cpgs of 96 VCSIP samples
beta <- read.delim(paste0(inpath, "model2.vcsip.raw_beta_matrix.txt"),
    header = TRUE)

# subset beta matrix to nonsmokers only
beta_nonsmokers <- beta[rownames(beta) %in%
    rownames(pd[pd$status == "nonsmoker", ]), ]
# subset beta matrix to smokers only
beta_smokers <- beta[rownames(beta) %in%
    rownames(pd[pd$status == "smoker", ]), ]

# calculate average beta value for smokers and nonsmokers separately
nonsmoker_mean_df <- calc_avg(beta_nonsmokers)
smoker_mean_df <- calc_avg(beta_smokers)

# calculate weighted average
model2_wt_avg <- calc_wt_avg(nonsmoker_mean_df,
    smoker_mean_df,
    smoke_rate = 0.123)

# export
write.table(model2_wt_avg,
    "model2_vcsip_weighted_avg_table.txt",
    sep="\t",
    col.names=T,
    row.names=T,
    quote=F)
    
#--------------------------------------------------------#
# Model 3 create reference table

inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model3/")

# read in the raw beta matrix for model 1 cpgs of 96 VCSIP samples
beta <- read.delim(paste0(inpath, "model3.vcsip.raw_beta_matrix.txt"),
    header = TRUE)

# subset beta matrix to nonsmokers only
beta_nonsmokers <- beta[rownames(beta) %in%
    rownames(pd[pd$status == "nonsmoker", ]), ]
# subset beta matrix to smokers only
beta_smokers <- beta[rownames(beta) %in%
    rownames(pd[pd$status == "smoker", ]), ]

# calculate average beta value for smokers and nonsmokers separately
nonsmoker_mean_df <- calc_avg(beta_nonsmokers)
smoker_mean_df <- calc_avg(beta_smokers)

# calculate weighted average
model3_wt_avg <- calc_wt_avg(nonsmoker_mean_df,
    smoker_mean_df,
    smoke_rate = 0.123)

# export
write.table(model3_wt_avg,
    "model3_vcsip_weighted_avg_table.txt",
    sep="\t",
    col.names=T,
    row.names=T,
    quote=F)
