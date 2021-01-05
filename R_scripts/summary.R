#!/usr/bin/env Rscript

##############################################################################
#LOAD LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES

#LOAD LIBRARIES
library(phyloseq)

#ASSIGN USER INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################


##############################################################################
#SUMMARY FILE

#OPEN FILE
sink(paste(output_dir, "summary_stats.txt", sep=""))
print(paste("Name of Community One: ", community_one_name, sep=""))
print(paste("Name of Community Two: ", community_two_name, sep=""))
#PRINT DESCRIPTIVE STATISTICS
print("Summary Statistics")
print(paste("Total Number of ASVs Prefiltered: ", merge_unfiltered_length, sep=""))
print(paste("Total Number of ASVs Postfiltered: ", merge_filtered_all, sep=""))
print(paste("Total Number of ASVs Removed: ", merge_unfiltered_length - merge_filtered_all, sep=""))
print(paste("Total Number of ASVs Removed for < minimum sum: ", merge_unfiltered_length - merge_filtered_min_sum_length, sep=""))
print(paste("Total Number of ASVs Removed for < minimum sample: ", merge_filtered_min_sum_length - merge_filtered_all))
print("")
print("")
print(paste("Number of ASVs in Community One Prefiltering: ", community_one_pre_length, sep=""))
print(paste("Number of ASVs in Community One Postfiltering: ", community_one_post_length, sep=""))
print(paste("Number of ASVs in Community Two Prefiltering: ", community_two_pre_length, sep=""))
print(paste("Number of ASVs in Community Two Postfiltering: ", community_two_post_length, sep=""))

#DETERMINE NUMBER OF ASVS UNIQUE TO EACH COMMUNITY OR SHARED BETWEEN COMMUNITIES
unique_one <- 0
shared <- 0
for(item in rownames(community_one)){
    #FOR EACH ASV IN COMMUNITY ONE
    if(item %in% rownames(community_two)){
        #IS THIS ASV IN COMMUNITY TWO?
        shared <- shared + 1
    }
    else{
        unique_one <- unique_one + 1
    }
}
unique_two <- 0
for(item in rownames(community_two)){
    if(!(item %in% rownames(community_one))){
        unique_two <- unique_two + 1
    }
}
print(paste("Number of ASVs unique to Community One: ", unique_one, sep=""))
print(paste("Number of ASVs unique to Community Two: ", unique_two, sep=""))
print(paste("Number of Shared ASVs: ", shared, sep=""))
print("")
print("")


print(paste("Number of ASVs Differentially Expressd Higher in Community One: ", community_one_differential_all, sep=""))
print(paste("Number of ASVs Differentially Expressd Higher in Community One (significant): ", community_one_differential_significant, sep=""))
print(paste("Number of ASVs Differentially Expressed Higher in Community Two: ", community_two_differential_all, sep=""))
print(paste("Number of ASVs Differentially Expressed Higher in Community Two (significant): ", community_two_differential_significant, sep=""))
print(paste("Number of ASVs Differentially Expressed Higher in Community One also found in Community Two: ", community_one_differential_in_two, sep=""))
print(paste("Number of ASVs Differentially Expressed Higher in Community Two also found in Community One: ", community_two_differential_in_one, sep=""))

print("")
print("")

print(paste("ASV with largest Log 2 Fold Change in Community One: ", community_one_log_fold_max, sep=""))
print(paste("ASV with largest Log 2 Fold Change in Community Two: ", community_two_log_fold_max, sep=""))
print(paste("ASV with smallest Differentially Abundant adjusted p-value in Community One: ", community_one_p_min, sep=""))
print(paste("ASV with smallest Differentially Abundant adjusted p-value in Community Two: ", community_two_p_min, sep=""))


#FUNCTION THAT COUNTS THE NUMBER OF UNIQUE CLASSIFICATIONS AT EACH TAXONOMIC LEVEL
count_tax_units <- function(phy_obj, level){
    temp <- tax_table(phy_obj)
    rownames(temp) <- NULL
    levels <- unique(na.omit(temp[,level]))
    return(length(levels))
}


print(paste("Total Number of Unique Kingdoms: ", count_tax_units(community_merge_norm_phy_obj, "Kingdom"), sep=""))
print(paste("Total Number of Unique Kingdoms in Community One: ", count_tax_units(community_one_norm_phy_obj, "Kingdom"), sep=""))
print(paste("Total Number of Unique Kingdoms in Community Two: ", count_tax_units(community_two_norm_phy_obj, "Kingdom"), sep=""))


print(paste("Total Number of Unique Phylum: ", count_tax_units(community_merge_norm_phy_obj, "Phylum"), sep=""))
print(paste("Total Number of Unique Phylum in Community One: ", count_tax_units(community_one_norm_phy_obj, "Phylum"), sep=""))
print(paste("Total Number of Unique Phylum in Community Two: ", count_tax_units(community_two_norm_phy_obj, "Phylum"), sep=""))

print(paste("Total Number of Unique Classes: ", count_tax_units(community_merge_norm_phy_obj, "Class"), sep=""))
print(paste("Total Number of Unique Classes in Community One: ", count_tax_units(community_one_norm_phy_obj, "Class"), sep=""))
print(paste("Total Number of Unique Classes in Community Two: ", count_tax_units(community_two_norm_phy_obj, "Class"), sep=""))

print(paste("Total Number of Unique Orders: ", count_tax_units(community_merge_norm_phy_obj, "Order"), sep=""))
print(paste("Total Number of Unique Orders in Community One: ", count_tax_units(community_one_norm_phy_obj, "Order"), sep=""))
print(paste("Total Number of Unique Orders in Community Two: ", count_tax_units(community_two_norm_phy_obj, "Order"), sep=""))

print(paste("Total Number of Unique Families: ", count_tax_units(community_merge_norm_phy_obj, "Family"), sep=""))
print(paste("Total Number of Unique Families in Community One: ", count_tax_units(community_one_norm_phy_obj, "Family"), sep=""))
print(paste("Total Number of Unique Families in Community Two: ", count_tax_units(community_two_norm_phy_obj, "Family"), sep=""))

print(paste("Total Number of Unique Genus: ", count_tax_units(community_merge_norm_phy_obj, "Genus"), sep=""))
print(paste("Total Number of Unique Genus in Community One: ", count_tax_units(community_one_norm_phy_obj, "Genus"), sep=""))
print(paste("Total Number of Unique Genus in Community Two: ", count_tax_units(community_two_norm_phy_obj, "Genus"), sep=""))

#CLOSE FILE
sink()

#PRINT TRACKER
print("Summary File Created")