

# conda activate tidymodels
library(tidyverse)
library(ggrepel)
library(readxl)

source("my_fxns.R")

# read in the metadata table to get status
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))
pd.keep <- pd[ ,c("Sample_Name", "status")]

# read in the beta matrix of all cpgs from all 3 models
beta.sub <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "summarize_models/psi_scores/three_models_nonzero_cpgs_vcsip_beta_matrix_with_rownames.txt"),
  header=T)

# read in the nonzero coef cpgs summary file for all 5 models
res <- read_xlsx(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "summarize_models/cpg_tables/all_model_nonzero_cpgs_annot.xlsx"))
colnames(res)[1] <- "term"

# make new res where cpgs with NA in col of interest are removed
df <- res[!(is.na(res$model1_EPIC_coef)), c(1,2)]
res.model1 <- df[df$model1_EPIC_coef != 0, ]
df <- res[!(is.na(res$model2_450k_coef)), c(1,3)]
res.model2 <- df[df$model2_450k_coef != 0, ]
df <- res[!(is.na(res$model3_PACE_coef)), c(1,4)]
res.model3 <- df[df$model3_PACE_coef != 0, ]

#----------------------------#

# Model 1

# Put the rows of beta.sub (cpgs) in the same order as the cpgs in res
beta.sub2 <- beta.sub[rownames(beta.sub) %in% res.model1$term, ]
beta.sub2 <- beta.sub2[match(res.model1$term, rownames(beta.sub2)), ]
identical(rownames(beta.sub2), res.model1$term)

# read in model 1 reference table for normalizaing betas
ref <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "predict_new/model1_EPIC_reference_table.txt"), header=T)

# take transpose to put cpgs as columns
mydat <- as.data.frame(t(beta.sub2))
# calculate psi on raw betas
my_psis <- calculate_psi(mydat, ref)

# normalize the betas
mydat.norm <- norm_betas(mydat, ref)
# calculate psi on normalized betas
norm_psis <- calculate_psi(mydat.norm, ref)

# clean and merge
colnames(my_psis) <- c("Sample_Name", "psi")
colnames(norm_psis) <- c("Sample_Name", "norm_psi")
myres <- merge(my_psis, norm_psis, by="Sample_Name")
myres2 <- merge(myres, pd.keep, by="Sample_Name")

# export
write.table(myres2,
  "model1_vcsip_psi_scores_and_norm_psi_scores.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#----------------------------#
# Model 2

# Put the rows of beta.sub (cpgs) in the same order as the cpgs in res
beta.sub2 <- beta.sub[rownames(beta.sub) %in% res.model2$term, ]
beta.sub2 <- beta.sub2[match(res.model2$term, rownames(beta.sub2)), ]
identical(rownames(beta.sub2), res.model2$term)

# read in model 2 reference table for normalizaing betas
ref <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "predict_new/model2_450k_reference_table.txt"), header=T)

# take transpose to put cpgs as columns
mydat <- as.data.frame(t(beta.sub2))
# calculate psi on raw betas
my_psis <- calculate_psi(mydat, ref)

# normalize the betas
mydat.norm <- norm_betas(mydat, ref)
# calculate psi on normalized betas
norm_psis <- calculate_psi(mydat.norm, ref)

# clean and merge
colnames(my_psis) <- c("Sample_Name", "psi")
colnames(norm_psis) <- c("Sample_Name", "norm_psi")
myres <- merge(my_psis, norm_psis, by="Sample_Name")
myres2 <- merge(myres, pd.keep, by="Sample_Name")

# export
write.table(myres2,
  "model2_vcsip_psi_scores_and_norm_psi_scores.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#----------------------------#
# Model 3

# Put the rows of beta.sub (cpgs) in the same order as the cpgs in res
beta.sub2 <- beta.sub[rownames(beta.sub) %in% res.model3$term, ]
beta.sub2 <- beta.sub2[match(res.model3$term, rownames(beta.sub2)), ]
identical(rownames(beta.sub2), res.model3$term)

# read in model 2 reference table for normalizaing betas
ref <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "predict_new/model3_PACE_reference_table.txt"), header=T)

# take transpose to put cpgs as columns
mydat <- as.data.frame(t(beta.sub2))
# calculate psi on raw betas
my_psis <- calculate_psi(mydat, ref)

# normalize the betas
mydat.norm <- norm_betas(mydat, ref)
# calculate psi on normalized betas
norm_psis <- calculate_psi(mydat.norm, ref)

# clean and merge
colnames(my_psis) <- c("Sample_Name", "psi")
colnames(norm_psis) <- c("Sample_Name", "norm_psi")
myres <- merge(my_psis, norm_psis, by="Sample_Name")
myres2 <- merge(myres, pd.keep, by="Sample_Name")

# export
write.table(myres2,
  "model3_vcsip_psi_scores_and_norm_psi_scores.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)
