#!/usr/bin/env Rscript

##############################################################################
#LOAD ALL LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES
library(fossil)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)

#ASSIGN USER INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################

#chao1 - species richness estimator will accounting for library size
#shannon - gives weight to evenness
#simpson - gives weight to abundance/species richness

##############################################################################
#DETERMINE ALPHA DIVERSITY

#CREATE EMPTY LISTS TO APPEND TO
chao1_results <- NULL
shannon_results <- NULL
shannon_equitability <- NULL
simpson_results <- NULL
for(i in 1:length(community_merge)){
  #FOR EACH SAMPLE IN COMMUNITY ONE DETERMINE ITS CHAO1, SHANON, AND SIMPSON RESULT THEN APPEND TO RESPECTIVE LIST
  chao1_results <- append(chao1_results, chao1(community_merge[,i], taxa.row=TRUE))
  shannon_results <- append(shannon_results, vegan::diversity(community_merge[,i], index="shannon"))
  sample_size <- sum(community_merge[,i] != 0)
  shannon_equitability <- append(shannon_equitability, shannon_results[i] / log(sample_size))
  simpson_results <- append(simpson_results, vegan::diversity(community_merge[,i], index="simpson"))
}

#ASSIGN NEW LISTS INTO METADATA DATAFRAME
metadata$chao1 <- chao1_results
metadata$shannon <- shannon_results
metadata$shannon_equitability <- shannon_equitability
metadata$simpson <- simpson_results
metadata$inverse_simpson <- 1 - metadata$simpson

##############################################################################


##############################################################################
#WRTE ALPHA DIVERSITY SIGNIFICANCE TO FILES

metadata$sample_group <- as.factor(metadata$sample_group)

#OPEN ALPHA DIVERSITY FILE
sink(file=paste(output_dir, "/alpha_diversity/alpha_diversity_significance.txt", sep=""))

#CHAO1 INDEX
print("Chao1")
chao1_significance <- t.test(metadata$chao1 ~ metadata$sample_group)
chao1_significance

#SHANNON INDEX
print("Shannon")
shannon_significance <- t.test(metadata$shannon ~ metadata$sample_group)
shannon_significance

print("Shannon Equitability Index")
shannon_equitability_significance <- t.test(metadata$shannon_equitability ~ metadata$sample_group)
shannon_equitability_significance

#SIMPSON INDEX
print("Simpson")
simpson_significance <- t.test(metadata$simpson ~ metadata$sample_group)
simpson_significance

print("Inverse Simpson")
inverse_simpson_significance <- t.test(metadata$inverse_simpson ~ metadata$sample_group)
inverse_simpson_significance

#CLOSE ALPHA DIVERSITY FILE
sink()

##############################################################################

##############################################################################
#PLOT ALPHA DIVERSITY

#SHANNON BOTH COMMUNITIES NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Shannon Diversity Index / Community Comparison / Notched") +
  geom_boxplot(aes(sample_group, shannon), notch=TRUE) +
  geom_jitter(aes(sample_group, shannon), width=0.1, height=0)

ggsave(filename="shannon_diversity_index_notched.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#sHANNON BOTH COMMUNITIES NON-NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Shannon Diversity Index / Community Comparison") +
  geom_boxplot(aes(sample_group, shannon), notch=FALSE) +
  geom_jitter(aes(sample_group, shannon), width=0.1, height=0)

ggsave(filename="shannon_diversity_index.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SHANNON EQUITABILITY NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Shannon Equitability Index / Community Comparison") +
  geom_boxplot(aes(sample_group, shannon_equitability), notch=TRUE) +
  geom_jitter(aes(sample_group, shannon_equitability), width=0.1, height=0)

ggsave(filename="shannon_equitability_index_notched.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SHANNON EQUITABILITY NON-NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Shannon Equitability Index / Community Comparison / Notched") +
  geom_boxplot(aes(sample_group, shannon_equitability), notch=FALSE) +
  geom_jitter(aes(sample_group, shannon_equitability), width=0.1, height=0)

ggsave(filename="shannon_equitability_index.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#CHAO1 BOTH COMMUNITIES NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Chao1 Diversity Index / Community Comparison / notched") +
  geom_boxplot(aes(sample_group, chao1), notch=TRUE) +
  geom_jitter(aes(sample_group, chao1), width=0.1, height=0)

ggsave(filename="chao1_diversity_index_notched.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#CHAO1 BOTH COMMUNITIES NON-NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Chao1 Diversity Index / Community Comparison") +
  geom_boxplot(aes(sample_group, chao1), notch=FALSE) +
  geom_jitter(aes(sample_group, chao1), width=0.1, height=0)

ggsave(filename="chao1_diversity_index.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SIMPSON BOTH COMMUNITIES NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Simpson Diversity Index / Community Comparison / notched") +
  geom_boxplot(aes(sample_group, simpson), notch=TRUE) +
  geom_jitter(aes(sample_group, simpson), width=0.1, height=0)

ggsave(filename="simpson_diversity_index_notched.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SIMPSON BOTH COMMUNITIES NON-NOTHCED
plot <- ggplot(data=metadata) +
  ggtitle("Simpson Diversity Index / Community Comparison") +
  geom_boxplot(aes(sample_group, simpson)) +
  geom_jitter(aes(sample_group, simpson), width=0.1, height=0)

ggsave(filename="simpson_diversity_index.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#INVERSE SIMPSON BOTH COMMUNITIES NOTCHED
plot <- ggplot(data=metadata) +
  ggtitle("Inverse Simpson Diversity Index / Community Comparison / notched") +
  geom_boxplot(aes(sample_group, inverse_simpson), notch=TRUE) +
  geom_jitter(aes(sample_group, inverse_simpson), width=0.1, height=0)

ggsave(filename="inverse_simpson_diversity_index_notched.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#INVERSE SIMPSON BOTH COMMUNITIES NON-NOTHCED
plot <- ggplot(data=metadata) +
  ggtitle("Inverse Simpson Diversity Index / Community Comparison") +
  geom_boxplot(aes(sample_group, inverse_simpson)) +
  geom_jitter(aes(sample_group, inverse_simpson), width=0.1, height=0)

ggsave(filename="inverse_simpson_diversity_index.png", plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

metadata$sample_group <- as.character(metadata$sample_group)

#SHANNON COMMUNITY ONE
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_one_name),]) +
  ggtitle(paste("Shannon Diversity Index / ", community_one_name, sep="")) +
  geom_boxplot(aes(0, shannon), width=10) + 
  geom_jitter(aes(0, shannon), height=0, width=3) +
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("shannon_diversity_index_",community_one_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SHANNON COMMUNITY TWO
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_two_name),]) +
  ggtitle(paste("Shannon Diversity Index / ", community_two_name, sep="")) +
  geom_boxplot(aes(0, shannon), width=10) + 
  geom_jitter(aes(0, shannon), height=0, width=3) + 
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("shannon_diversity_index_",community_two_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SHANNON EQUITABILITY COMMUNITY ONE
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_one_name),]) +
  ggtitle(paste("Shannon Equitability Index / ", community_one_name, sep="")) +
  geom_boxplot(aes(0, shannon_equitability), width=10) + 
  geom_jitter(aes(0, shannon_equitability), height=0, width=3) +
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("shannon_equitability_index_",community_one_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SHANNON EQUITABILITY COMMUNITY TWO
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_two_name),]) +
  ggtitle(paste("Shannon Equitability Index / ", community_two_name, sep="")) +
  geom_boxplot(aes(0, shannon_equitability), width=10) + 
  geom_jitter(aes(0, shannon_equitability), height=0, width=3) + 
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("shannon_equitability_index_",community_two_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#CHAO1 COMMUNITY ONE
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_one_name),]) +
  ggtitle(paste("Chao1 Diversity Index / ", community_one_name, sep="")) +
  geom_boxplot(aes(0, chao1), width=10) + 
  geom_jitter(aes(0, chao1), height=0, width=3) +
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("Chao1_diversity_index_",community_one_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#CHAO1 COMMUNITY TWO
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_two_name),]) +
  ggtitle(paste("Chao1 Diversity Index / ", community_two_name, sep="")) +
  geom_boxplot(aes(0, chao1), width=10) + 
  geom_jitter(aes(0, chao1), height=0, width=3) +
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("Chao1_diversity_index_",community_two_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SIMPSON COMMUNITY ONE
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_one_name),]) +
  ggtitle(paste("Simpson Diversity Index / ", community_one_name, sep="")) +
  geom_boxplot(aes(0, simpson), width=10) + 
  geom_jitter(aes(0, simpson), height=0, width=3) + 
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("simpson_diversity_index_",community_one_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#SIMPSON COMMUNITY TWO
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_two_name),]) +
  ggtitle(paste("Simpson Diversity Index / ", community_two_name, sep="")) +
  geom_boxplot(aes(0, simpson), width=10) + 
  geom_jitter(aes(0, simpson), height=0, width=3) +
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("simpson_diversity_index_", community_two_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#iNVERSE SIMPSON COMMUNITY ONE
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_one_name),]) +
  ggtitle(paste("Inverse Simpson Diversity Index / ", community_one_name, sep="")) +
  geom_boxplot(aes(0, inverse_simpson), width=10) + 
  geom_jitter(aes(0, inverse_simpson), height=0, width=3) + 
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("inverse_simpson_diversity_index_",community_one_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

#INVERSE SIMPSON COMMUNITY TWO
plot <- ggplot(data=metadata[ which(metadata$sample_group == community_two_name),]) +
  ggtitle(paste("Inverse Simpson Diversity Index / ", community_two_name, sep="")) +
  geom_boxplot(aes(0, inverse_simpson), width=10) + 
  geom_jitter(aes(0, inverse_simpson), height=0, width=3) +
  xlim(c(-10,10)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste("inverse_simpson_diversity_index_", community_two_name, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/alpha_diversity/", sep=""))

##############################################################################

#WRITE METADATA TO FILE
write.table(metadata, file=paste(output_dir, "/alpha_diversity/", "alpha_diversity_metadata.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#SAVE R WORKSPACE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT TRACKER
print("alpha diversity done")