
# Create a reference table for each model
# Columns: <cpg> <estimate> <mean> <sd> <y_int>

#--------------------------------------------------------#

# Model 1 create reference table

inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model1/")

# read in the nonzero cpgs and their coefficients, and y-intercept value
ref <- read.delim(paste0(inpath,
    "nCpG_120_lambda_3.79269019073225e-09_nonzero_predictors.txt"),
    header = TRUE)
colnames(ref)[1] <- "cpg"

# read in the mean_sd table
mean_sd <- read.delim(paste0(inpath, "model1.mean_sd_table.txt"),
    header = TRUE)

# save y-intercept value
my_y <- ref$estimate[ref$cpg == "(Intercept)"]

# merge ref (removing last row which is intercept and last col which is penalty)
# and mean_sd tables on "cpg"
mytab <- merge(ref[-nrow(ref), -ncol(ref)], mean_sd, by = "cpg")

# add y-int column
mytab$y_int <- my_y

# export this model 1 reference table
write.table(mytab,
  "model1_EPIC_reference_table.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE)

#-------------------------------------------------------#

# Model 2 create reference table

inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model2/ncpg120/")

# read in the nonzero cpgs and their coefficients, and y-intercept value
ref <- read.delim(paste0(inpath,
    "nCpG_120_lambda_3.79269019073225e-09_nonzero_predictors.txt"),
    header = TRUE)
colnames(ref)[1] <- "cpg"

# read in the mean_sd table
mean_sd <- read.delim(paste0(inpath, "model2.mean_sd_table.txt"),
    header = TRUE)

# save y-intercept value
my_y <- ref$estimate[ref$cpg == "(Intercept)"]

# merge ref (removing last row which is intercept and last col which is penalty)
# and mean_sd tables on "cpg"
mytab <- merge(ref[-nrow(ref), -ncol(ref)], mean_sd, by = "cpg")

# add y-int column
mytab$y_int <- my_y

# export this model 1 reference table
write.table(mytab,
  "model2_450k_reference_table.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE)

#-------------------------------------------------------#

# Model 3 create reference table

inpath <- paste0("/home/groups/hoolock2/u0/bd/",
    "Projects/lyndsey_shorey_project/classification_models/",
    "fall_2022_final_models/model3/")

# read in the nonzero cpgs and their coefficients, and y-intercept value
ref <- read.delim(paste0(inpath,
    "nCpG_98_lambda_1.27427498570313e-08_nonzero_predictors.txt"),
    header = TRUE)
colnames(ref)[1] <- "cpg"

# read in the mean_sd table
mean_sd <- read.delim(paste0(inpath, "model3.mean_sd_table.txt"),
    header = TRUE)

# save y-intercept value
my_y <- ref$estimate[ref$cpg == "(Intercept)"]

# merge ref (removing last row which is intercept and last col which is penalty)
# and mean_sd tables on "cpg"
mytab <- merge(ref[-nrow(ref), -ncol(ref)], mean_sd, by = "cpg")

# add y-int column
mytab$y_int <- my_y

# export this model 1 reference table
write.table(mytab,
  "model3_PACE_reference_table.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE)
