# Colour Polymorphism Manuscript
# 4. Data analyses - genomics (LD decay and permutations)

rm(list=ls())

library(openxlsx)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)
library(stringr)
library(ggrepel)
library(HardyWeinberg)

setwd("C:/Users/feiner/Dropbox/DataAnalysis")

LD_BCO2 <- read.csv("DataS6_LD_from_BCO2_100kb.csv", header=T)

gene_start <- 26165779  
gene_end <- 26188197
exons <- data.frame(
  exon_start = c(26165779,26175083,26176676,26179731,26181132,26184487,26185863,26187111,26187895),
  exon_end = c(26165863,26175287,26176887,26179846,26181234,26184618,26186023,26187278,26188197))

############### Fig. S1A - Patterns of LD decay - BCO2

pdf("./Plots/LD_BCO2.pdf", height = 4, width = 8, useDingbats=FALSE)
ggplot(data=LD_BCO2, aes(x=BP_B, y=R2)) + 
  geom_point(alpha=0.3, shape=16, size=5) + theme_bw() + 
  xlim(26161882-27000,26161882+27000) + 
  geom_vline(xintercept=26161882, linetype="solid", color = "red") +
  geom_vline(xintercept=26161882+2500, linetype="dashed", color = "red") +
  geom_vline(xintercept=26161882-2500, linetype="dashed", color = "red") +
  geom_vline(xintercept=26164143, linetype="dashed", color = "blue") +
  #geom_vline(xintercept=26169143, linetype="dashed", color = "blue") +
  # Add gene as a horizontal line
  geom_segment(aes(x=gene_start, xend=gene_end, y=1.05, yend=1.05), 
               color="black", lwd=1) +
  # Add exons as thicker segments (or rectangles)
  geom_segment(data=exons, aes(x=exon_start, xend=exon_end, y=1.05, yend=1.05), 
               color="black", lwd=10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

LD_SPR <- read.csv("DataS7_LD_from_SPR_100kb.csv", header=T)

gene_start <- 77997110  
gene_end <- 77994329
exons <- data.frame(
  exon_start = c(77997110,77995004,77994519),
  exon_end = c(77996786,77994714,77994329))

############### Fig. S1B - Patterns of LD decay - SPR

pdf("./Plots/LD_SPR.pdf", height = 4, width = 8, useDingbats=FALSE)
ggplot(data=LD_SPR, aes(x=BP_B, y=R2)) + 
  geom_point(alpha=0.3, shape=16, size=5) + theme_bw() + 
  xlim(77999982-27000,77999982+27000) + 
  geom_vline(xintercept=77999982, linetype="solid", color = "red") +
  geom_vline(xintercept=77999982+2500, linetype="dashed", color = "red") +
  geom_vline(xintercept=77999982-2500, linetype="dashed", color = "red") +
  # Add gene as a horizontal line
  geom_segment(aes(x=gene_start, xend=gene_end, y=1.05, yend=1.05), 
               color="black", lwd=1) +
  # Add exons as thicker segments (or rectangles)
  geom_segment(data=exons, aes(x=exon_start, xend=exon_end, y=1.05, yend=1.05), 
               color="black", lwd=10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


#Calculate likelihood of getting overlap with one of the morph loci
  sharedExpected <- array() 
for (j in 1:10000) {
SET1 <- sample(135441, 1222, replace=F) #FST IT
SET2 <- sample(121615, 1021, replace=F) #FST SA
SET3 <- sample(142465, 6434, replace=F) #RAD IT
SET4 <- sample(144629, 3244, replace=F) #RAD SA
SET5 <- sample(269548, 647, replace=F) #Tajima IT
SET6 <- sample(261491, 625, replace=F) #Tajima SA
sharedExpected[j] <- sum(intersect(SET1,1),intersect(SET2,1),intersect(SET3,1),intersect(SET4,1),intersect(SET5,1),intersect(SET6,1))
}

obs_higher1 <- length(sharedExpected[sharedExpected>=1])
obs_higher1/100 #in percent

hist(sharedExpected) #distribution of shared transcripts expected by chance
quantile(sharedExpected, c(0.025, 0.975)) #this gives the upper and lower 95% confidence interval
