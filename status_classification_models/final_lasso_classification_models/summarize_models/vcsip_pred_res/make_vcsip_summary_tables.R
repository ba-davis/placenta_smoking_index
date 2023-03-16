

# Model 1

# read in results summary table
df <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/model1/model1.prediction_table.txt"),
  header=T)

# read in the psi/norm_psi table
psi <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/summarize_models/psi_scores/",
  "normalized_psi_scores/model1_vcsip_psi_scores_and_norm_psi_scores.txt"),
  header=T)
psi2 <- psi[ ,-4]

# merge
res <- merge(df, psi2, by="Sample_Name")
# export
write.table(res,
  "model1_vcsip_summary_table.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#-----------------------------------------------------
# Model 2

# read in results summary table
df <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/model2/ncpg120/model2.prediction_table.txt"),
  header=T)

# read in the psi/norm_psi table
psi <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/summarize_models/psi_scores/",
  "normalized_psi_scores/model2_vcsip_psi_scores_and_norm_psi_scores.txt"),
  header=T)
psi2 <- psi[ ,-4]

# merge
res <- merge(df, psi2, by="Sample_Name")
# export
write.table(res,
  "model2_vcsip_summary_table.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#-----------------------------------------------------
# Model 3

# read in results summary table
df <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/model3/model3.prediction_table.txt"),
  header=T)

# read in the psi/norm_psi table
psi <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/summarize_models/psi_scores/",
  "normalized_psi_scores/model3_vcsip_psi_scores_and_norm_psi_scores.txt"),
  header=T)
psi2 <- psi[ ,-4]

# merge
res <- merge(df, psi2, by="Sample_Name")
# export
write.table(res,
  "model3_vcsip_summary_table.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)
