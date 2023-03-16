
# conda activate tidymodels

library(tidyverse)
library(readxl)
library(ggrepel)

myfiles <- list.files(getwd(), pattern="psi_scores.txt$", full=T)

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

#------------------#
# raw psi scores

# add model name to psi colname for each df
for (i in 1:length(dfs)) {
  colnames(dfs[[i]])[2] <- names(dfs)[i]
}

# combine dfs into one df
my_table <- Reduce(function(df1, df2) merge(df1, df2, by = "Sample_Name", all = TRUE), dfs)
my_table2 <- my_table[ ,c(1,4,2,5,8)]

# melt to long format
plotdf <- pivot_longer(my_table2, model1_EPIC:model3_PACE, names_to="model", values_to="psi")
colnames(plotdf)[2] <- "status"

# label outliers
# outlier function to label outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#plotdf <- plotdf %>%
#  group_by(model) %>%
#  group_by(status.x) %>%
#  mutate(outlier = ifelse(is_outlier(psi), Sample_Name, NA))
# Not working? doesn't label all outlier points, ignore for now

p <- ggplot(plotdf, aes(x=model, y=psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_PSIscore_by_status_vcsip_boxplot_outliers.png", p, dpi=600)

p <- ggplot(plotdf, aes(x=model, y=psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_PSIscore_by_status_vcsip_boxplot.png", p, dpi=600)

#--------------------------------------#

# normalized psi

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

# add model name to psi colname for each df
for (i in 1:length(dfs)) {
  colnames(dfs[[i]])[3] <- names(dfs)[i]
}

# combine dfs into one df
my_table <- Reduce(function(df1, df2) merge(df1, df2, by = "Sample_Name", all = TRUE), dfs)
my_table2 <- my_table[ ,c(1,4,3,6,9)]

# melt to long format
plotdf <- pivot_longer(my_table2, model1_EPIC:model3_PACE, names_to="model", values_to="norm_psi")
colnames(plotdf)[2] <- "status"

# label outliers
# outlier function to label outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))}

plotdf <- plotdf %>%
  group_by(model) %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))
# Not working? doesn't label all outlier points, ignore for now

p <- ggplot(plotdf, aes(x=model, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_normPSIscore_by_status_vcsip_boxplot_outliers.png", p, dpi=600)

p <- ggplot(plotdf, aes(x=model, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_normPSIscore_by_status_vcsip_boxplot.png", p, dpi=600)

















# this isn't adding outliers for model1 or model2
# manually add

plotdf$outlier[plotdf$pd.Sample_Name=="flora" & plotdf$pd.status.x=="nonsmoker" & plotdf$model=="model1_EPIC"] <- "flora"
plotdf$outlier[plotdf$pd.Sample_Name=="gaston" & plotdf$pd.status.x=="nonsmoker" & plotdf$model=="model1_EPIC"] <- "gaston"
plotdf$outlier[plotdf$pd.Sample_Name=="jockey" & plotdf$pd.status.x=="smoker" & plotdf$model=="model1_EPIC"] <- "jockey"

plotdf$outlier[plotdf$pd.Sample_Name=="advent" & plotdf$pd.status.x=="smoker" & plotdf$model=="model2_450k"] <- "advent"
plotdf$outlier[plotdf$pd.Sample_Name=="icebox" & plotdf$pd.status.x=="smoker" & plotdf$model=="model2_450k"] <- "icebox"
plotdf$outlier[plotdf$pd.Sample_Name=="bubble" & plotdf$pd.status.x=="smoker" & plotdf$model=="model2_450k"] <- "bubble"
plotdf$outlier[plotdf$pd.Sample_Name=="alcove" & plotdf$pd.status.x=="smoker" & plotdf$model=="model2_450k"] <- "alcove"
plotdf$outlier[plotdf$pd.Sample_Name=="jockey" & plotdf$pd.status.x=="smoker" & plotdf$model=="model2_450k"] <- "jockey"

plotdf$outlier[plotdf$pd.Sample_Name=="flora" & plotdf$model=="model3_PACE"] <- "flora"
plotdf$outlier[plotdf$pd.Sample_Name=="doc" & plotdf$model=="model3_PACE"] <- "doc"
plotdf$outlier[plotdf$pd.Sample_Name=="gaston" & plotdf$model=="model3_PACE"] <- "gaston"
plotdf$outlier[plotdf$pd.Sample_Name=="advent" & plotdf$model=="model3_PACE"] <- "advent"
plotdf$outlier[plotdf$pd.Sample_Name=="jockey" & plotdf$model=="model3_PACE"] <- "jockey"

p <- ggplot(plotdf, aes(x=model, y=psi, fill=pd.status.x)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_PSIscore_by_status_vcsip_boxplot_outliers.png", p, dpi=600)

p <- ggplot(plotdf, aes(x=model, y=psi, fill=pd.status.x)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_PSIscore_by_status_vcsip_boxplot.png", p, dpi=600)

#ggsave("combined_PSIscore_by_status_vcsip_boxplot_outliers.pdf", p)
