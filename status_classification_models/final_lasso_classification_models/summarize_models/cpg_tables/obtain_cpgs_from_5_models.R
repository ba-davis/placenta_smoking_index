
# conda activate tidymodels

library(tidymodels)
library(writexl)

file1="../../model1/imp_vars.txt"
file2="../../model2/ncpg120/imp_vars.txt"
file3="../../model3/imp_vars.txt"

myfiles=c(file1, file2, file3)

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

for (i in 1:length(dfs)) {
  colnames(dfs[[i]])[2] <- paste0(names(dfs)[i], "_coef")
  dfs[[i]] <- dfs[[i]][ ,c(1,2)]
}

my_table <- Reduce(function(df1, df2) merge(df1, df2, by = "term", all = TRUE), dfs)
# remove intercept column
my_table <- my_table[-1,]
colnames(my_table)[1] <- "IlmnID"

#write_xlsx(my_table, "all_model_cpgs.xlsx")

# read in cpg annotation file
annot <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/annot/EPIC_hg37_smoking_index_annot.txt"))

# merge the tables
newdf <- merge(my_table, annot, by="IlmnID")
write_xlsx(newdf, "all_model_cpgs_annot.xlsx")

# now only keep rows with nonzero coef in at least 1 model
res2 <- newdf[rowSums(newdf[, c(2:4)], na.rm=TRUE) != 0, ]
write_xlsx(res2, "all_model_nonzero_cpgs_annot.xlsx")
write.table(res2,
  "all_model_nonzero_cpgs_annot.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

