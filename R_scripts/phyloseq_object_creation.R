#!/usr/bin/env Rscript

##############################################################################
#LOAD LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES

#LOAD LIBRARIES
library(phyloseq)
library(microbiome)
library(tidyverse)

#ASSIGN USER INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################


##############################################################################
#PHYLOSEQ OBJECT CREATION

#LOAD TAXONOMY FILE OUTPUTTED BY DADA2
tax_file <- paste(output_dir, "/taxmat_samples.txt", sep="")
tax_assignment <- read.table(tax_file)

#CONVERT COUNT TABLE TO OTU TABLE
community_merge_otu <- otu_table(community_merge, taxa_are_rows=TRUE)

#CONVERT TAXONOMY TABLE TO TAX TABLE
tax <- tax_table(as.matrix(tax_assignment))

#CREATE NON-NORMALIZED PHYLOSEQ OBJECT CONTAINING BOTH COMMUNITIES
community_merge_phy_obj <- phyloseq(community_merge_otu, tax)

#CONVERT COMMUNITY ONE COUNT DATA TO OTU TABLE
community_one_otu <- otu_table(community_one, taxa_are_rows=TRUE)

#CREATE COMMUNITY ONE NON-NORMALIZED PHYLOSEQ OBJECT
community_one_phy_obj <- phyloseq(community_one_otu, tax)

#CONVERT COMMUNITY TWO COUNT DATA TO OTU TABLE
community_two_otu <- otu_table(community_two, taxa_are_rows=TRUE)

#CREATE COMMUNITY TWO NON-NORMALIZED PHYLOSEQ OBJECT
community_two_phy_obj <- phyloseq(community_two_otu, tax)

#CONVERT NORMALIZED COUNT DATA TO OTU TABLE
community_merge_otu_norm <- otu_table(community_merge_norm,taxa_are_rows=TRUE)

#CREATE NORMALIZED PHYLOSEQ OBJECT CONTAINING BOTH COMMUNITIES
community_merge_norm_phy_obj <- phyloseq(community_merge_otu_norm, tax)

#CONVERT NORMALIZED COMMUNITY ONE COUNT DATA TO OTU TABLE
community_one_otu_norm <- otu_table(community_one_norm, taxa_are_rows=TRUE)

#CREATE NORMALIZED COMMUNITY ONE PHYLOSEQ OBJECT
community_one_norm_phy_obj <- phyloseq(community_one_otu_norm, tax)

#CONVERT NORMALIZED COMMUNITY TWO COUNT DATA TO OTU TABLE
community_two_otu_norm <- otu_table(community_two_norm, taxa_are_rows=TRUE)

#CREATE NORMALIZED COMMUNITY TWO PHYLOSEQ OBJECT
community_two_norm_phy_obj <- phyloseq(community_two_otu_norm, tax)


#TRANSFORM TO PERCENT ABUNDANCE USING NORMALIZED PHYLOSEQ OBJECTS
community_one_phy_obj_perc_abundance <- transform(community_one_norm_phy_obj, 'compositional')
community_two_phy_obj_perc_abundance <- transform(community_two_norm_phy_obj, 'compositional')
community_merge_phy_obj_perc_abundance <- transform(community_merge_norm_phy_obj, 'compositional')

##############################################################################
#RETURN PERCENT OF ASVS IDENTIFIED TO EACH TAXONOMIC LEVEL
plot_identification_levels <- function(tax_table){
  tax_table %>%
    as("matrix") %>%
    as_tibble(rownames="OTU") %>%
    gather("Rank", "Name", rank_names(tax_table)) %>%
    group_by(Rank) %>%
    summarize(OTUs_classified = sum(!is.na(Name))) %>%
    mutate(Frac_classified = OTUs_classified / ntaxa(tax_table)) %>%
    mutate(Rank = factor(Rank, rank_names(tax_table))) %>%
    arrange(Rank)
}

#SHOW ALL ASSIGNMENTS OF UNFILTERED DATA
all_asv <- plot_identification_levels(tax)
#WRITE TO FILE
write.table(all_asv, file=paste(output_dir, "total_asv_taxonomic_assignments_summary.txt", sep=""), row.names=FALSE, col.names=TRUE)


#SHOW ASSIGNMENTS OF FILTERED DATA
tax_table_used <- tax[row.names(community_merge)]
used_asvs <- plot_identification_levels(tax_table_used)
#WRITE TO FILE
write.table(used_asvs, file=paste(output_dir, "used_asv_taxonomic_assignments.txt", sep=""), row.names=FALSE, col.names=TRUE)
##############################################################################

#SAVE R WORKSPACE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT TRACKER
print("Phyloseq Object Creation Done")