
# conda activate tidymodels

library(tidyverse)
library(readxl)
library(ggrepel)

myfiles <- list.files(getwd(), pattern="psi_scores.txt$", full=T)

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

#---------------------------------------#

# label outliers

# outlier function to label outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#-------------------------------#
# Model 1

# add outlier column to df
d2 <- dfs[[1]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model1 EPIC PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model1_PSIscore_by_status_vcsip_boxplot_outliers.pdf", p)

p <- ggplot(d2, aes(x=status, y=psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model1 EPIC PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model1_PSIscore_by_status_vcsip_boxplot.pdf", p)

# now do normalized psi scores
d2 <- dfs[[1]]
colnames(d2)[1] <- "Sample_Name"
d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model1 EPIC Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model1_normPSIscore_by_status_vcsip_boxplot_outliers.pdf", p)

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model1 EPIC Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model1_normPSIscore_by_status_vcsip_boxplot.pdf", p)


#-----------------------------#
# Model 2

# add outlier column
d2 <- dfs[[2]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model2 450k PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model2_PSIscore_by_status_vcsip_boxplot_outliers.pdf", p)

p <- ggplot(d2, aes(x=status, y=psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model2 450k PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model2_PSIscore_by_status_vcsip_boxplot.pdf", p)

# now do for norm_psi

d2 <- dfs[[2]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model2 450k Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model2_normPSIscore_by_status_vcsip_boxplot_outliers.pdf", p)

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model2 450k Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model2_normPSIscore_by_status_vcsip_boxplot.pdf", p)

#-----------------------------#
# Model 3

# add outlier
d2 <- dfs[[3]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model3 PACE PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model3_PSIscore_by_status_vcsip_boxplot_outliers.pdf", p)

p <- ggplot(d2, aes(x=status, y=psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model3 PACE PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model3_PSIscore_by_status_vcsip_boxplot.pdf", p)

# now do for norm_psi

d2 <- dfs[[3]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model3 PACE Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model3_normPSIscore_by_status_vcsip_boxplot_outliers.pdf", p)

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data Model3 PACE Norm PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("model3_normPSIscore_by_status_vcsip_boxplot.pdf", p)
