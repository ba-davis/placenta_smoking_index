


# read in cpgs file
df <- read.delim("nCpG_120_lambda_3.79269019073225e-09_nonzero_predictors.txt", header=T)
# remove the intercept column
df2 <- df[!(df$term=="(Intercept)"), ]
colnames(df2)[1] <- "IlmnID"

# read in cpg annotation file
annot <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/annot/EPIC_hg37_smoking_index_annot.txt"))

# merge the tables
newdf <- merge(df2, annot, by="IlmnID")
# remove penalty column
newdf$penalty <- NULL

# export merged table
write.table(newdf,
  "model1_nonzero_predictors_annot.txt",
  sep="\t", col.names=T, row.names=F, quote=F)
