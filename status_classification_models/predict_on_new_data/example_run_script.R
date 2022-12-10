


source("predict_functions.R")

#--------------------------------------#
# User Variables

# 1. input beta matrix ()
inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model1/")

# read in example cpg beta matrix
# (this is "df_clean", 96 samples as rows, cpgs returned from limma as columns
# and 2 other cols for sample name and status)
beta <- read.delim(paste0(inpath, "model1.clean_input_betas.txt"),
    header = TRUE)
# remove last two columns leaving only cpgs as columns
input_beta <- beta[, -c(103, 104)]

# 2. model reference table
mytab <- read.delim("model1_EPIC_reference_table.txt")

#-----------------------------------------#

res <- predict_on_new_data(input_beta, mytab)
