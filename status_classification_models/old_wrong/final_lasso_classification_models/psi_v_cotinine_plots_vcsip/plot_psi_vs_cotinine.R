
# Plot PSI score against Cotinine values and Birth Weight

# conda activate tidymodels
library(tidyverse)
library(ggrepel)

infile <- "../psi_scores/model1.psi_vcsip.txt"
mod_num <- "model1"
y_rho <- 0
y_pval <- 0.05

infile <- "../psi_scores/model2.psi_vcsip.txt"
mod_num <- "model2"
y_rho <- -0.2
y_pval <- -0.21

infile <- "../psi_scores/model3.psi_vcsip.txt"
mod_num <- "model3"
y_rho <- -0.9
y_pval <- -1.0

infile <- "../psi_scores/model4.psi_vcsip.txt"
mod_num <- "model4"
y_rho <- -2
y_pval <- -2.1

infile <- "../psi_scores/model5.psi_vcsip.txt"
mod_num <- "model5"
y_rho <- -0.6
y_pval <- -0.63

#----------------------------#

# read in the pd file
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))

# read in the cotinine metadata
cotinine_dat <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
    "lyndsey_shorey_project/metadata/Cotinine_ng_ml_pd.txt"), header=T)
colnames(cotinine_dat)[1] <- "Sample_Name"


# read in the d2 file showing psi score per sample
mypath <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/spring_2022_final_models/",
  "summarize_models/psi_scores")

d2 <- read.delim(paste0(mypath, "/", infile), header=T)
colnames(d2)[1] <- "Sample_Name"
colnames(d2)[2] <- "smoking_status"
colnames(d2)[3] <- "PSI"

identical(d2$Sample_Name, pd$Sample_Name)
identical(cotinine_dat$Sample_Name, d2$Sample_Name)

#d2$infant_birth_weight <- pd$infant_birth_weight

# add cotinine data to d2
d2$Rand_Cotinine_ng_ml <- cotinine_dat$Rand_Cotinine_ng_ml
d2$Mid_Cotinine_ng_ml <- cotinine_dat$Mid_Cotinine_ng_ml
d2$Late_Cotinine_ng_ml <- cotinine_dat$Late_Cotinine_ng_ml

# plot PSI vs birth weight
#p <- ggplot(d2, aes(x=infant_birth_weight, y=PSI, color=smoking_status)) +
#  geom_point() +
#  theme_classic() +
#  labs(x="Infant Birth Weight", y="PSI Score", title="VCSIP Data PSI Score by Infant Birth Weight") +
#  theme(plot.title = element_text(hjust = 0.5))
#ggsave("model1_PSIscore_by_birthweight_vcsip_scatterplot.pdf", p)

#-------------------------------------------------------------------#

# get pearson correlation coefficient and pval
rho <- paste0("r=", round(cor.test(d2$PSI, d2$Rand_Cotinine_ng_ml, method="pearson")$estimate, 4))
pval <- paste0("p=", round(cor.test(d2$PSI, d2$Rand_Cotinine_ng_ml, method="pearson")$p.value, 7))

# plot PSI vs Rand cotinine
p <- ggplot(d2, aes(x=Rand_Cotinine_ng_ml, y=PSI, color=smoking_status)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate("text", x=50, y=y_rho, label=rho) +
  annotate("text", x=50, y=y_pval, label=pval) +
  theme_classic() +
  labs(x="Rand Cotinine ng/ml", y="PSI Score", title=paste0(mod_num, " VCSIP PSI Score by Rand Cotinine")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(mod_num, "_PSIscore_by_randCot_vcsip_scatterplot.pdf"), p)

#-----------#

rho <- paste0("r=", round(cor.test(d2$PSI, d2$Mid_Cotinine_ng_ml, method="pearson")$estimate, 4))
pval <- paste0("p=", round(cor.test(d2$PSI, d2$Mid_Cotinine_ng_ml, method="pearson")$p.value, 7))

# plot PSI vs Rand cotinine
p <- ggplot(d2, aes(x=Mid_Cotinine_ng_ml, y=PSI, color=smoking_status)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate("text", x=50, y=y_rho, label=rho) +
  annotate("text", x=50, y=y_pval, label=pval) +
  theme_classic() +
  labs(x="Mid Cotinine ng/ml", y="PSI Score", title=paste0(mod_num, " VCSIP PSI Score by Mid Cotinine")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(mod_num, "_PSIscore_by_midCot_vcsip_scatterplot.pdf"), p)

#-----------#

rho <- paste0("r=", round(cor.test(d2$PSI, d2$Late_Cotinine_ng_ml, method="pearson")$estimate, 4))
pval <- paste0("p=", round(cor.test(d2$PSI, d2$Late_Cotinine_ng_ml, method="pearson")$p.value, 7))

# plot PSI vs Rand cotinine
p <- ggplot(d2, aes(x=Late_Cotinine_ng_ml, y=PSI, color=smoking_status)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate("text", x=50, y=y_rho, label=rho) +
  annotate("text", x=50, y=y_pval, label=pval) +
  theme_classic() +
  labs(x="Late Cotinine ng/ml", y="PSI Score", title=paste0(mod_num, " VCSIP PSI Score by Late Cotinine")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(mod_num, "_PSIscore_by_lateCot_vcsip_scatterplot.pdf"), p)
