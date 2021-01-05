#!/usr/bin/env Rscript

##############################################################################
#LIBRARY LOAD AND USER INPUT ASSIGNED TO SCRIPT VARIABLES

#LOAD LIBRARIES
library("phyloseq")
library("tidyverse")
library("DESeq2")

#ASSIGN INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################

#CREATE A COPY OF METADATA TABLE
metadata_rearranged <- metadata
#ASSIGN FILENAMES A NEW COLUMN INSTEAD OF BEING THE ROWNAMES
metadata_rearranged$file_name <- rownames(metadata)
#CHANGE THE ROWNAMES TO BE SAMPLE NAMES
rownames(metadata_rearranged) <- metadata$sample_name

##############################################################################
#IDENTIFICATION FUNCTION DETERMINES HOW MANY ASVS OF EACH CLASSIFICATION THERE ARE AT EACH TAXONOMIC LEVEL


identification <- function(phy_obj, level){
    #CONGLOMERATE A PHYLOSEQ OBJECT TO GIVEN TAXONOMY (LEVEL)
    temp <- tax_glom(phy_obj, level)
    #EXTRACT THE COUNTS FROM NEW PHYLOSEQ OBJECT
    counts <- data.frame(otu_table(temp))
    #CORRECT THE ROWNAMES
    rownames(counts) <- data.frame(tax_table(temp))[,level]
    #WRITE NORMALIZED COUNTS FOR EVERY CLASSFICIATION AT EACH TAXONOMIC LEVEL
    write.table(counts, file=paste(output_dir, "/taxonomic_summary/", level, "_counts_normalized.txt", sep=""), sep="\t", row.names=TRUE, col.names=TRUE, append=FALSE)

    #SPLIT COUNTS BY SAMPLE GROUP
    community_one_id <- counts[,sample_one_names]
    community_two_id <- counts[,sample_two_names]
    
    #CREATE NEW DATAFRAME CONTAINING ROWSUMS OF COMMUNITY ONE (SUM OF EACH GIVEN CLASSIFIER AT SPECIFIED TAXONOMIC LEVEL)
    summary_stats <- data.frame(rowSums(community_one_id))
    #ADD THE MEAN OF COMMUNITY ONE
    summary_stats[,2] <- rowMeans(community_one_id)
    #ADD THE STANDARD DEVIATION OF COMMUNITY ONE
    summary_stats[,3] <- apply(community_one_id, 1, sd)
    #ADD THE ROWSUMS OF COMMUNITY TWO
    summary_stats[,4] <- rowSums(community_two_id)
    #ADD ROWMEANS OF COMMUNITY TWO
    summary_stats[,5] <- rowMeans(community_two_id)
    #ADD STANDARD DEVIATION OF COMMUNITY TWO
    summary_stats[,6] <- apply(community_two_id, 1, sd)
    t_test <- NULL
    
    #CONDUCT A T-TEST COMPARING EACH CLASSIFICATION BETWEEN COMMUNITY ONE AND COMMUNITY TWO
    for(i in c(1:length(rownames(community_one_id)))){
        p_value <- t.test(community_one_id[i,], community_two_id[i,])$p.value
        t_test <- append(t_test, p_value)
    }

    #ADD T-TEST TO DATAFRAME
    summary_stats[,7] <- t_test
    
    #ASSIGN COLUMN NAMES TO NEW DATAFRAME
    colnames(summary_stats) <- c(paste(community_one_name, "_sum"), paste(community_one_name, "_mean"), paste(community_one_name, "_standard_deviation"), paste(community_two_name, "_sum"), paste(community_two_name, "_mean"), paste(community_two_name, "_standard_deviation"), "T_Test_significance")
    #WRITE THIS DATAFRAME TO NEW FILE
    write.table(summary_stats, file=paste(output_dir, "/taxonomic_summary/", level, "_summary_stats.txt", sep=""), sep="\t", row.names=TRUE, col.names=TRUE, append=FALSE)
}


#CALL ABOVE FUNCTION FOR EACH TAXONOMIC LEVEL
identification(community_merge_norm_phy_obj, "Kingdom")
identification(community_merge_norm_phy_obj, "Phylum")
identification(community_merge_norm_phy_obj, "Class")
identification(community_merge_norm_phy_obj, "Order")
identification(community_merge_norm_phy_obj, "Family")
identification(community_merge_norm_phy_obj, "Genus")

#PRINT TRACKER
print("Identification Summary Done")