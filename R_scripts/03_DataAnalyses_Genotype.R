# Colour Polymorphism Manuscript
# 3. Data analyses - genotypic data


rm(list=ls())

# load packages
library(openxlsx)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)
library(stringr)
library(ggrepel)
library(HardyWeinberg)
library(tinytable)
library(car)
library(ggpubr)

# set working directory
#setwd("C:/Users/Tobias/Dropbox/Tobias Documents/Lund/Papers/WallLizardPapers/ColourPolymorphism_LossIntrogression/DataAnalysis")
setwd("C:/Users/feiner/Dropbox/DataAnalysis")

#Read in the data
#Do yellow first
GT_results_yellow_pheno <- read.csv("DataS4_Genotyping_BCO2_yellow.csv")

GT_results_yellow_pheno <- GT_results_yellow_pheno %>% 
  mutate(Match = case_when(Genotype == "y/y" & grepl("yellow", morph) ~ "match", #condition 1
                           Genotype == "Y/y" & !grepl("yellow", morph) ~ "match", # Condition 2
                           Genotype == "Y/Y" & !grepl("yellow", morph) ~ "match", # Condition 3
                                  TRUE ~ "mismatch"))

# For this purpose, treat IT as only the pure IT and hybrid as SA (there are only 2 populations from
# the hybrid zone proper)

GT_results_yellow_pheno <- GT_results_yellow_pheno %>% 
  mutate(lineage2 = case_when(lineage == "SA" ~ "SA",
                              lineage == "Hybrid" ~ "SA", 
                              lineage == "IT" ~ "IT"))

GT_results_yellow_pheno$lineage <- as.factor(GT_results_yellow_pheno$lineage)
GT_results_yellow_pheno$lineage2 <- as.factor(GT_results_yellow_pheno$lineage2)
GT_results_yellow_pheno$Match <- as.factor(GT_results_yellow_pheno$Match)
GT_results_yellow_pheno$abbpop <- as.factor(GT_results_yellow_pheno$abbpop)

#explore the relationship between mismatches and predictors Greenness and Lineage:

Green_mismatch_yellow <- GT_results_yellow_pheno %>%  group_by(GreenResolved) %>%   # Group by GreenResolved
  summarize(total = n(),                                   # Total count per GreenResolved
    mismatch_count = sum(Match == "mismatch"),     # Count of "mismatch"
    mismatch_percentage = (mismatch_count / total) * 100)  # Percentage of "mismatch"

# generate tables of mismatches across relevant categories
GT_results_yellow_pheno %>% 
  group_by(lineage2, abbpop) %>% # remove for overall %
  summarize(total = n(),                             
            match_count = sum(Match == "match"),     # Count of "mismatch"
            match_percentage = (match_count / total) * 100) %>% tt()

GT_results_yellow_pheno %>% 
  group_by(Genotype, morph) %>% 
  summarise(Individuals=n()) %>% tt()

GT_results_yellow_pheno %>% 
  group_by(Genotype) %>% 
  summarise(Individuals=n()) %>% tt()

GT_results_yellow_pheno %>% 
  group_by(Genotype, morph, lineage2) %>% 
  summarise(Individuals=n(),
            SVL=mean(svl),
            Green=mean(GreenResolved, na.rm = TRUE)) %>% tt()

# same for orange
GT_results_orange_pheno <- read.csv("DataS5_Genotyping_SPR_orange.csv")

GT_results_orange_pheno <- GT_results_orange_pheno %>% 
  mutate(Match = case_when(
    Genotype == "o/o" & grepl("orange", morph) ~ "match",  # Condition 1
    Genotype == "O/o" & !grepl("orange", morph) ~ "match", # Condition 2
    Genotype == "O/O" & !grepl("orange", morph) ~ "match", # Condition 3
    TRUE ~ "mismatch"))

# For this purpose, treat IT as only the pure IT and SA+hybrid as the alternative category
GT_results_orange_pheno <- GT_results_orange_pheno %>% 
  mutate(lineage2 = case_when(lineage == "SA" ~ "SA",
                              lineage == "Hybrid" ~ "SA", 
                              lineage == "IT" ~ "IT"))

GT_results_orange_pheno$lineage <- as.factor(GT_results_orange_pheno$lineage)
GT_results_orange_pheno$lineage2 <- as.factor(GT_results_orange_pheno$lineage2)
GT_results_orange_pheno$Match <- as.factor(GT_results_orange_pheno$Match)
GT_results_orange_pheno$abbpop <- as.factor(GT_results_orange_pheno$abbpop)

# explore the relationship between mismatches and predictors Greenness and Lineage

Green_mismatch_orange <- GT_results_orange_pheno %>%  group_by(GreenResolved) %>%   # Group by GreenResolved
  summarize(total = n(),                                   # Total count per GreenResolved
            mismatch_count = sum(Match == "mismatch"),     # Count of "mismatch"
            mismatch_percentage = (mismatch_count / total) * 100)  # Percentage of "mismatch"


# generate tables of mismatches across relevant categories
GT_results_orange_pheno %>% 
  group_by(lineage2) %>% # remove for overall %
  summarize(total = n(),                             
            match_count = sum(Match == "match"),     # Count of "mismatch"
            match_percentage = (match_count / total) * 100) %>% tt()

GT_results_orange_pheno %>% 
  group_by(Genotype, morph) %>% 
  summarise(Individuals=n()) %>% tt()

GT_results_orange_pheno %>% 
  group_by(Genotype) %>% 
  summarise(Individuals=n()) %>% tt()

GT_results_orange_pheno %>% 
  group_by(Genotype, morph) %>% 
  summarise(Individuals=n(),
            SVL=mean(svl),
            Green=mean(GreenResolved, na.rm = TRUE)) %>% tt()

Mismatch_Yellow <- subset(GT_results_yellow_pheno, Match=="mismatch")
Mismatch_Orange <- subset(GT_results_orange_pheno, Match=="mismatch")
Mismatch_Yellow <- Mismatch_Yellow[,c(1,2,10,5,7,8,6,3,9)]
Mismatch_Orange <- Mismatch_Orange[,c(1,2,10,5,7,8,6,3,9)]
colnames(Mismatch_Yellow)[colnames(Mismatch_Yellow) == 'Genotype'] <- "Genotype_yellow"
colnames(Mismatch_Yellow)[colnames(Mismatch_Yellow) == 'Match'] <- "Match_yellow"
colnames(Mismatch_Orange)[colnames(Mismatch_Orange) == 'Genotype'] <- "Genotype_orange"
colnames(Mismatch_Orange)[colnames(Mismatch_Orange) == 'Match'] <- "Match_orange"

Mismatch_Yellow_Orange <- merge(Mismatch_Yellow, GT_results_orange_pheno, by="ID", all.x=T)
Mismatch_Orange_Yellow <- merge(Mismatch_Orange, GT_results_yellow_pheno, by="ID", all.x=T)

Mismatch_Yellow_Orange <- Mismatch_Yellow_Orange[,c(1:9,11,17)]
Mismatch_Orange_Yellow <- Mismatch_Orange_Yellow[,c(1:7,11,17,8,9)]

colnames(Mismatch_Yellow_Orange)[colnames(Mismatch_Yellow_Orange) == 'Genotype'] <- "Genotype_orange"
colnames(Mismatch_Yellow_Orange)[colnames(Mismatch_Yellow_Orange) == 'Match'] <- "Match_orange"
colnames(Mismatch_Orange_Yellow)[colnames(Mismatch_Orange_Yellow) == 'Genotype'] <- "Genotype_yellow"
colnames(Mismatch_Orange_Yellow)[colnames(Mismatch_Orange_Yellow) == 'Match'] <- "Match_yellow"

Mismatch_All <- rbind(Mismatch_Yellow_Orange,Mismatch_Orange_Yellow)

write.csv(Mismatch_All, "TableS2_Genotype_Phenotype_Mismatches.csv", row.names = F)

#############
#split the genotypes into alleles
GT_results_yellow_pheno[c('Allele1', 'Allele2')] <- str_split_fixed(GT_results_yellow_pheno$Genotype, '/', 2)
GT_results_orange_pheno[c('Allele1', 'Allele2')] <- str_split_fixed(GT_results_orange_pheno$Genotype, '/', 2)

#convert to long format
GT_results_yellow_long <- gather(GT_results_yellow_pheno, Allele_id, Alleles, Allele1:Allele2)
GT_results_orange_long <- gather(GT_results_orange_pheno, Allele_id, Alleles, Allele1:Allele2)

#summarize allele frequencies per population
GT_results_yellow_sum <- as.data.frame(GT_results_yellow_long %>%
                                         group_by(abbpop, Alleles) %>%
                                         dplyr::summarize(Freq=n()) %>%
                                         mutate(countT= sum(Freq)) %>%
                                         mutate(per=paste0(round(Freq/countT,2))))
GT_results_orange_sum <- as.data.frame(GT_results_orange_long %>%
                                         group_by(abbpop, Alleles) %>%
                                         dplyr::summarize(Freq=n()) %>%
                                         mutate(countT= sum(Freq)) %>%
                                         mutate(per=paste0(round(Freq/countT,2))))

GT_results_Y_sum <- subset(GT_results_yellow_sum, Alleles == "Y")
GT_results_Y_sum <- GT_results_Y_sum %>% rename(Y = Freq)

GT_results_O_sum <- subset(GT_results_orange_sum, Alleles == "O")
GT_results_O_sum <- GT_results_O_sum %>% rename(O = Freq)

Pheno_data <- read.csv("DataS3_141_Populations_ForPolymorphism.csv")

# We only need a few variables
Pheno_data <- Pheno_data[,c(1,14,18)]

#merge datasets
Yellow <- merge(GT_results_Y_sum, Pheno_data, by="abbpop")
Orange <- merge(GT_results_O_sum, Pheno_data, by="abbpop")

Yellow$lineage <- as.factor(Yellow$lineage)
Yellow$per <- as.numeric(Yellow$per)
Yellow$Green <- as.numeric(Yellow$Green)
str(Yellow)

Orange$lineage <- as.factor(Orange$lineage)
Orange$per <- as.numeric(Orange$per)
Orange$Green <- as.numeric(Orange$Green)
str(Orange)

levels(Yellow$lineage)[levels(Yellow$lineage) == "Hybrid"] <- "SA"
levels(Orange$lineage)[levels(Orange$lineage) == "Hybrid"] <- "SA"

Yellow_IT <- subset(Yellow, lineage == "IT")
Yellow_SA <- subset(Yellow, lineage == "SA")

Orange_IT <- subset(Orange, lineage == "IT")
Orange_SA <- subset(Orange, lineage == "SA")

##################### YELLOW #####################################
# Fit the binomial glm
model_allele_Y <- glm(cbind(Y, countT - Y) ~ Green*lineage, # testing for interaction first
                 family = binomial(link = "logit"), 
                 data = Yellow)

Anova(model_allele_Y, type="III")
summary(model_allele_Y)
# note that the interaction is non-significant
# if one wishes to remove the interaction to fit a single slope the analysis can be run without it

# model_allele_Y_2 <- glm(cbind(Y, countT - Y) ~ Green + lineage, 
#                      family = binomial(link = "logit"), 
#                      data = Yellow)

# Anova(model_allele_Y_2, type="II")
# summary(model_allele_Y_2)

# some model diagnostics
with(summary(model_allele_Y), c(Deviance = deviance, Null_Deviance = null.deviance, df = df.residual))

par(mfrow = c(2, 2))  # Set up a 2x2 plotting grid
plot(model_allele_Y)
overdispersion <- sum(residuals(model_allele_Y, type = "pearson")^2) / df.residual(model_allele_Y)
overdispersion


# Create a data frame for all combinations of Green and lineage
new_data_y <- expand.grid(
  Green = seq(min(Yellow$Green), max(Yellow$Green), length.out = 1000),
  lineage = levels(Yellow$lineage)
)

# Predict on the link scale
predictions_Y <- predict(model_allele_Y, newdata = new_data_y, type = "link", se.fit = TRUE)

# Add predictions and standard errors to the data frame
new_data_y$fit <- predictions_Y$fit  # Predicted logit values
new_data_y$se.fit <- predictions_Y$se.fit  # Standard errors

# Define the inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Calculate confidence intervals on the logit scale
new_data_y$lower <- inv_logit(new_data_y$fit - 1.96 * new_data_y$se.fit)
new_data_y$upper <- inv_logit(new_data_y$fit + 1.96 * new_data_y$se.fit)

# Transform fitted values to probabilities
new_data_y$predicted_Y <- inv_logit(new_data_y$fit)

Geno_freq_y <- ggplot() +
  geom_line(data = new_data_y, aes(x = Green, y = predicted_Y, color = lineage), linewidth = 1) +
  geom_ribbon(data = new_data_y, aes(x = Green, ymin = lower, ymax = upper, fill = lineage), alpha = 0.2, color = NA) +
  geom_point(data = Yellow, aes(x = Green, y = per, color = lineage, size = countT/2, fill = lineage), shape = 16, alpha=0.8) +
  scale_color_manual(values = c("#77AADD", "#EE8866")) +
  scale_fill_manual(values = c("#77AADD", "#EE8866")) +
  labs(x = "Greenness", y = "Frequency of dominant Y allele") + ylim(0,1) +
  scale_x_continuous(limits=c(1, 10), breaks=c(1,5,10)) +
  theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 8))

Geno_freq_y 

min(Yellow$countT/2) #6
mean(Yellow$countT/2) #18.86
sd(Yellow$countT/2) #4.87

############################# ORANGE ##################

# Fit the binomial GLM
model_allele_O <- glm(cbind(O, countT - O) ~ Green*lineage,  # testing for interaction
                      family = binomial(link = "logit"), 
                      data = Orange)
Anova(model_allele_O, type = "III")
summary(model_allele_O)

# note that the interaction is non-significant
# if one wishes to remove the interaction to fit a single slope the analysis can be run without it

# model_allele_O_2 <- glm(cbind(O, countT - O) ~ Green + lineage,  # the interaction is non-significant
#                      family = binomial(link = "logit"), 
#                      data = Orange)

#Anova(model_allele_O_2, type = "II")
#summary(model_allele_O_2)

# some model diagnostics
with(summary(model_allele_O), c(Deviance = deviance, Null_Deviance = null.deviance, df = df.residual))

par(mfrow = c(2, 2))  # Set up a 2x2 plotting grid
plot(model_allele_Y)
overdispersion <- sum(residuals(model_allele_O, type = "pearson")^2) / df.residual(model_allele_O)
overdispersion

# Create a data frame for all combinations of Green and lineage
new_data_o <- expand.grid(
  Green = seq(min(Orange$Green), max(Orange$Green), length.out = 1000),
  lineage = levels(Orange$lineage)
)

# Predict on the link scale
predictions_O <- predict(model_allele_O, newdata = new_data_o, type = "link", se.fit = TRUE)

# Add predictions and standard errors to the data frame
new_data_o$fit <- predictions_O$fit  # Predicted logit values
new_data_o$se.fit <- predictions_O$se.fit  # Standard errors


# Define the inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Calculate confidence intervals on the logit scale
new_data_o$lower <- inv_logit(new_data_o$fit - 1.96 * new_data_o$se.fit)
new_data_o$upper <- inv_logit(new_data_o$fit + 1.96 * new_data_o$se.fit)

# Transform fitted values to probabilities
new_data_o$predicted_O <- inv_logit(new_data_o$fit)

Geno_freq_o <- ggplot() +
  geom_line(data = new_data_o, aes(x = Green, y = predicted_O, color = lineage), linewidth = 1) +
  geom_ribbon(data = new_data_o, aes(x = Green, ymin = lower, ymax = upper, fill = lineage), alpha = 0.2, color = NA) +
  geom_point(data = Orange, aes(x = Green, y = per, color = lineage, size = countT/2, fill = lineage), shape = 16, alpha=0.8) +
  scale_color_manual(values = c("#77AADD", "#EE8866")) +
  scale_fill_manual(values = c("#77AADD", "#EE8866")) +
  labs(x = "Greenness", y = "Frequency of dominant O allele") + ylim(0,1) +
  scale_x_continuous(limits=c(1, 10), breaks=c(1,5,10)) +
  theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 8))

Geno_freq_o

min(Orange$countT/2) #11
mean(Orange$countT/2) #20.22
sd(Orange$countT/2) #4.34

############### Fig. 3C and D - Genotype frequencies for Yellow and Orange

pdf("./Plots/GenotypeFreq_Greenness.pdf", useDingbats=FALSE)
ggarrange(Geno_freq_y, Geno_freq_o, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()
