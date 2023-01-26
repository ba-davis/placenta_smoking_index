

# conda activate tidymodels
library(tidyverse)
library(ggrepel)
library(readxl)

# read in the full beta matrix
beta <- read.delim(gzfile(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/new_beta.txt.gz")))

# read in the nonzero coef cpgs summary file for all 5 models
res <- read_xlsx(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
  "summarize_models/cpg_tables/all_model_nonzero_cpgs_annot.xlsx"))

colnames(res)[1] <- "term"
# subset the full beta matrix to the cpgs of interest
beta.sub <- beta[rownames(beta) %in% res$term, ]
# export this subsetted beta matrix for future use
write.table(beta.sub, "three_models_nonzero_cpgs_vcsip_beta_matrix.txt", sep="\t", col.names=T, row.names=F, quote=F)

# read in the metadata table
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))

# make new res where cpgs with NA in col of interest are removed
df <- res[!(is.na(res$model1_EPIC_coef)), c(1,2)]
res.model1 <- df[df$model1_EPIC_coef != 0, ]
df <- res[!(is.na(res$model2_450k_coef)), c(1,3)]
res.model2 <- df[df$model2_450k_coef != 0, ]
df <- res[!(is.na(res$model3_PACE_coef)), c(1,4)]
res.model3 <- df[df$model3_PACE_coef != 0, ]

#-----------------------------#

# Model 1

# Put the rows of beta.sub (cpgs) in the same order as the cpgs in res
beta.sub2 <- beta.sub[rownames(beta.sub) %in% res.model1$term, ]
beta.sub2 <- beta.sub2[match(res.model1$term, rownames(beta.sub2)), ]
identical(rownames(beta.sub2), res.model1$term)

# Make sure rows of pd (samples) are in same order as columns of beta.sub (samples)
identical(rownames(pd), colnames(beta.sub2))

## since I know the status of my own samples I can check to see if they cluster by score:

## y is a data frame multiplying each beta value by it's coefficient
y <- NULL
i <- NULL
A <- res.model1$model1_EPIC_coef
B <- beta.sub2
y <- matrix(nrow = nrow(B), ncol = ncol(B))
for (i in 1:nrow(B)){
    y[i,]<-rbind(as.numeric(B[i,])*A[i])
}

s <- colSums(y) #s is the linear combination of all beta values multiplied by their respective coefficient, giving a single score per sample

d2 <- data.frame(pd$Sample_Name, pd$status, s)

write.table(d2,
  "model1.psi_vcsip.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#---------------------#

# Model 2

# Put the rows of beta.sub (cpgs) in the same order as the cpgs in res
beta.sub2 <- beta.sub[rownames(beta.sub) %in% res.model2$term, ]
beta.sub2 <- beta.sub2[match(res.model2$term, rownames(beta.sub2)), ]
identical(rownames(beta.sub2), res.model2$term)

# Make sure rows of pd (samples) are in same order as columns of beta.sub (samples)
identical(rownames(pd), colnames(beta.sub2))

## since I know the status of my own samples I can check to see if they cluster by score:

## y is a data frame multiplying each beta value by it's coefficient
y <- NULL
i <- NULL
A <- res.model2$model2_450k_coef
B <- beta.sub2
y <- matrix(nrow = nrow(B), ncol = ncol(B))
for (i in 1:nrow(B)){
    y[i,]<-rbind(as.numeric(B[i,])*A[i])
}

s <- colSums(y) #s is the linear combination of all beta values multiplied by their respective coefficient, giving a single score per sample

d2 <- data.frame(pd$Sample_Name, pd$status, s)

write.table(d2,
  "model2.psi_vcsip.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#------------------------#

# Model 3

# Put the rows of beta.sub (cpgs) in the same order as the cpgs in res
beta.sub2 <- beta.sub[rownames(beta.sub) %in% res.model3$term, ]
beta.sub2 <- beta.sub2[match(res.model3$term, rownames(beta.sub2)), ]
identical(rownames(beta.sub2), res.model3$term)

# Make sure rows of pd (samples) are in same order as columns of beta.sub (samples)
identical(rownames(pd), colnames(beta.sub2))

## since I know the status of my own samples I can check to see if they cluster by score:

## y is a data frame multiplying each beta value by it's coefficient
y <- NULL
i <- NULL
A <- res.model3$model3_PACE_coef
B <- beta.sub2
y <- matrix(nrow = nrow(B), ncol = ncol(B))
for (i in 1:nrow(B)){
    y[i,]<-rbind(as.numeric(B[i,])*A[i])
}

s <- colSums(y) #s is the linear combination of all beta values multiplied by their respective coefficient, giving a single score per sample

d2 <- data.frame(pd$Sample_Name, pd$status, s)

write.table(d2,
  "model3.psi_vcsip.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#----------------------------------------#
