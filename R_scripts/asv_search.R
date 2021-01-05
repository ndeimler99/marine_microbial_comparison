#!/usr/bin/env Rscript

##############################################################################
#ASSIGN INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
asv <- args[2]

#LOAD WORKSPACE
load(paste(output_dir, "/r_save", sep=""))
##############################################################################

#CONVERT USER INPUT INTO ALL CAPITAL LETTERS AS DADA2 OUTPUTS
asv <- toupper(asv)

#OPEN FILE BASED ON ASV DESIRED
sink(file = paste(output_dir, asv,"_summary.txt", sep=""))

#PRINT ASV NAME
print(asv)

#SAVE ASV TAXONOMY
taxonomy <- tax_assignment[asv,]

#PRINT ASV SUMMARY RESULTS
print(paste("Kingdom:", taxonomy[asv, "Kingdom"]))
print(paste("Phylum:", taxonomy[asv, "Phylum"]))
print(paste("Class:", taxonomy[asv, "Class"]))
print(paste("Order:", taxonomy[asv, "Order"]))
print(paste("Family:", taxonomy[asv, "Family"]))
print(paste("Genus:", taxonomy[asv, "Genus"]))
print(paste("Species:", taxonomy[asv, "Species"]))
print(community_merge[asv,])
sink()