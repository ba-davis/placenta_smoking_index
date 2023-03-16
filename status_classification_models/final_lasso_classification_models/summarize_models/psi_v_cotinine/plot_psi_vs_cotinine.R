
# USE NORM_PSI

# Plot PSI score against Cotinine values and Birth Weight

# conda activate tidymodels
library(tidyverse)
library(ggrepel)

# read in the pd file
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))

# read in the cotinine metadata
cotinine_dat <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
    "lyndsey_shorey_project/metadata/Cotinine_ng_ml_pd.txt"), header=T)
colnames(cotinine_dat)[1] <- "Sample_Name"
rownames(cotinine_dat) <- cotinine_dat$Sample_Name

# read in the summary table 
mypath <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "summarize_models/summary_tables")

#--------------------------------------------------------------------#

source("source_fxn.R")

plot_psi_v_cot(infile = "model1_vcsip_summary_table.txt",
  pd = pd,
  cotinine_dat = cotinine_dat,
  mypath = mypath,
  mod_num = "model1_EPIC",
  y_rho = -9,
  y_pval = -12)

plot_psi_v_cot(infile = "model2_vcsip_summary_table.txt",
  pd = pd,
  cotinine_dat = cotinine_dat,
  mypath = mypath,
  mod_num = "model2_450k",
  y_rho = -9,
  y_pval = -11)

plot_psi_v_cot(infile = "model3_vcsip_summary_table.txt",
  pd = pd,
  cotinine_dat = cotinine_dat,
  mypath = mypath,
  mod_num = "model3_PACE",
  y_rho = -17,
  y_pval = -21)
