#!/usr/bin/env Rscript


##############################################################################
#LOAD ALL LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES

#LOAD LIBRARIES
library(edgeR)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(vegan)

#ASSIGN USER INPUT
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
norm_technique <- args[2]
libsize <- args[3]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################

##############################################################################
#NORMALIZATION PROCESS USES COMMUNITY MERGE THEN RESPLITS THIS DATAFRAME

#IF NO NORMALIZATION IS WANTED
if(norm_technique == "None"){
  community_merge_norm <- community_merge
  community_one_norm <- community_one
  community_two_norm <- community_two
} else {
  #IF NORMALIZATION IS WANTED

    
    if(norm_technique == "rarefaction"){
        #RAREFY 
        community_merge_norm <- as.data.frame(rrarefy(community_merge, libsize))
        print("rarefaction")
    }

    if(norm_technique == "DESeq2"){
        #DESEQ NORMALIZATION
        dds <- DESeqDataSetFromMatrix(countData = community_merge, colData = metadata, design=~sample_group)
        dds <- estimateSizeFactors(dds)
        community_merge_norm <- counts(dds, normalized=TRUE)
        community_merge_norm <- as.data.frame(community_merge_norm, col.names=1)
        colnames(community_merge_norm) <- colnames(community_merge)
        print("DESeq2 used")
    }

    if(norm_technique == "RPM"){
      #READS PER MILLION NORMALIZATION
      community_merge_norm <- community_merge
      colsums <- colSums(community_merge)

      for(column in 1:ncol(community_merge_norm)){
        community_merge_norm[,column] <- community_merge_norm[,column] / colsums[column] * 1000000
      }

      print("RPM Used")
    }

    #SPLIT MERGED DATAFRAME INTO TWO COMMUNITIES
    community_one_norm <- dplyr::select(community_merge_norm, sample_one_names)
    community_two_norm <- dplyr::select(community_merge_norm, sample_two_names)

    #REMOVE ANY ROWS CONTAINING ALL ZEROES
    community_one_norm <- community_one_norm[rowSums(community_one_norm[,-1])>0,]
    community_two_norm <- community_two_norm[rowSums(community_two_norm[,-1])>0,]
}

##############################################################################

##############################################################################
#PLOTTING SEQUENCE DEPTH BEFORE AND AFTER NORMALIZATION

#CALCULATE TOTAL NUMBER OF ASVS IN EACH COLUMN (SAMPLE)
sample_sums_non_norm <- colSums(community_merge)
sample_sums_norm <- colSums(community_merge_norm)

#ADD NEWELY CALCULATED SUMS TO METADATA
metadata$non_norm_sums <- sample_sums_non_norm
metadata$normalized_sums <- sample_sums_norm

#PLOT USING GGPLOT
non_normalized <- ggplot(data=metadata) +
  geom_bar(stat="identity", mapping =aes(x=sample_name,y=non_norm_sums)) +
  ggtitle("Non-Normalized")

normalized <- ggplot(data=metadata)+
  geom_bar(stat="identity", mapping=aes(x=sample_name, y=normalized_sums)) +
  ggtitle("Normalized")

#SAVE PLOTS TO FILES
ggsave(filename="non_normalized_seq_counts_plot.png", plot=non_normalized, device="png", path=paste(output_dir, "/normalization_results/", sep=""))
ggsave(filename="normalized_seq_counts_plot.png", plot=normalized, device="png", path=paste(output_dir, "/normalization_results/", sep=""))

##############################################################################

#WRITE NORMALIZED TABLES TO FILES FOR POTENTIAL USER REFERENCE
write.table(community_merge_norm, file=paste(output_dir, "/normalization_results/total_asv_counts_normalized.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(community_one_norm, file=paste(output_dir, "/normalization_results/", community_one_name, "_asv_counts_normalized.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(community_two_norm, file=paste(output_dir, "/normalization_results/", community_two_name, "_asv_counts_normalized.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#SAVE R WORKSPACE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT TRACKER
print("normalization done")