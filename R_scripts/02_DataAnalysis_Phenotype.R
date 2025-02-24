# Colour Polymorphism Manuscript
# 2. Data analyses - phenotypic data

# Aim: to test if the frequency of W morph within populations is explained by the population-average 
# expression of the nigriventris phenotype

# Note that both males and females are included for morph frequencies while the estimate of nigriventris 
# phenotype (greenness) is calculated from data on males

rm(list=ls())

# load libraries
library(openxlsx)
library(dplyr)
library(tibble)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggtern)
library(tinytable)
library(car)

# set working directory (pick one)
#setwd("C:/Users/Tobias/Dropbox/Tobias Documents/Lund/Papers/WallLizardPapers/ColourPolymorphism_LossIntrogression/DataAnalysis")
setwd("C:/Users/feiner/Dropbox/DataAnalysis")

# Read in data set
data <- read.csv("DataS3_141_Populations_ForPolymorphism.csv") 

str(data)
data$abbpop <- as.factor(data$abbpop)
data$lineage <- as.factor(data$lineage)

# we will do a single analysis that include both IT and SA

# For this purpose, treat IT as only the pure IT and SA+hybrid as the alternative category

data <- data %>% 
  mutate(lineage2 = case_when(lineage == "SA" ~ "SA",
                           lineage == "Hybrid" ~ "SA", 
                           lineage == "IT" ~ "IT"))

data$lineage2 <- as.factor(data$lineage2)

data %>% 
  group_by(lineage) %>% 
  summarise(Individuals=n()) %>% tt()

# Fit the binomial GLMs using lme4

model_1 <- glm(cbind(white, totalN - white) ~ Green*lineage2, 
               family = binomial(link = "logit"), 
               data = data)

Anova(model_1, type="III")

# note that the interaction is non-significant
# if one wishes to remove the interaction to fit a single slope the analysis can be run without it
# model_2 <- glm(cbind(white, totalN - white) ~ Green+lineage2, 
#               family = binomial(link = "logit"), 
#               data = data)

# Anova(model_2, type="II")
# summary(model_2)

# some model diagnostics
with(summary(model_1), c(Deviance = deviance, Null_Deviance = null.deviance, df = df.residual))

par(mfrow = c(2, 2))  # Set up a 2x2 plotting grid
plot(model_1)
overdispersion <- sum(residuals(model_1, type = "pearson")^2) / df.residual(model_1)
overdispersion

# Generate predicted values

# Create a data frame for all combinations of Green and lineage2
new_data <- expand.grid(
  Green = seq(min(data$Green), max(data$Green), length.out = 1000),
  lineage2 = levels(data$lineage2)
)

# Predict on the link scale
predictions <- predict(model_1, newdata = new_data, type = "link", se.fit = TRUE)

# Add predictions and standard errors to the data frame
new_data$fit <- predictions$fit  # Predicted logit values
new_data$se.fit <- predictions$se.fit  # Standard errors

# Define the inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Calculate confidence intervals on the logit scale
new_data$lower <- inv_logit(new_data$fit - 1.96 * new_data$se.fit)
new_data$upper <- inv_logit(new_data$fit + 1.96 * new_data$se.fit)

# Transform fitted values to probabilities
new_data$predicted <- inv_logit(new_data$fit)

############### Fig. 2B - Frequencies of white morphs for IT and SA lineages by greenness

#plot the fit and the observed data in one graph
pdf("./Plots/MorphFreq_Lineage_Greenness.pdf", height=4, width=7, useDingbats=T)
ggplot() +
  geom_line(data = new_data, aes(x = Green, y = predicted, color = lineage2), linewidth = 1) +
  geom_ribbon(data = new_data, aes(x = Green, ymin = lower, ymax = upper, fill = lineage2), alpha = 0.2, color = NA) +
  geom_point(data = data, aes(x = Green, y = freqW, color = lineage2, size = totalN, fill = lineage2), shape = 16, alpha=0.8) +
  scale_color_manual(values = c("#EE8866", "#77AADD")) +
  scale_fill_manual(values = c("#EE8866", "#77AADD")) +
  labs(x = "Greenness", y = "Frequency of white morph") + ylim(0,1) +
  scale_x_continuous(limits=c(1, 10), breaks=c(1,5,10)) +
  theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 8))
dev.off()

############### Fig. 2C and D - Ternary plots for IT and SA lineage

data$freqO_new <- data$freqO + data$freqOY/2
data$freqY_new <- data$freqY + data$freqOY/2

data_IT <- subset(data, lineage2 =="IT")
data_SA <- subset(data, lineage2 =="SA")

base <- ggtern(data = data_IT, aes(freqY_new, freqW, freqO_new)) +
  geom_mask() + 
  geom_point(position= position_jitter_tern(x=0.005, y=0.005, z=0.005), alpha = 0.7, size = 4, shape=16, aes(color=Green)) + 
  scale_color_gradientn(colors = c("#996633","#996633","#666600","#669900","#33CC00","#339900","#336600")) +
  xlab("Yellow")+ylab("White")+zlab("Orange") +
  theme_rgbw() + theme_light() + theme(legend.position="none")  + theme_showgrid()
zoom <- base + theme_zoom_T(0.2) + theme_hidegrid_major()
pdf("./Plots/Ternary_byLineage_IT.pdf", height=5, useDingbats=F)
grid.arrange(base,zoom, nrow=2, 
             ncol=1, heights=c(5,3))
dev.off()

base <- ggtern(data = data_SA, aes(freqY_new, freqW, freqO_new)) +
  geom_mask() + 
  geom_point(position= position_jitter_tern(x=0.005, y=0.005, z=0.005), alpha = 0.7, size = 4, shape=16, aes(color=Green)) + 
  scale_color_gradientn(colors = c("#996633","#996633","#666600","#669900","#33CC00","#339900","#336600")) +
  xlab("Yellow")+ylab("White")+zlab("Orange") +
  theme_rgbw() + theme_light() + theme(legend.position="none")
zoom <- base + theme_zoom_T(0.2) + theme_hidegrid_major()
base_grob <- ggplotGrob(base)
zoom_grob <- ggplotGrob(zoom)
pdf("./Plots/Ternary_byLineage_SA.pdf", height=5, useDingbats=T)
grid.arrange(base_grob,zoom_grob, nrow=2, 
             heights=c(5,3))
dev.off()
