#!/usr/bin/env Rscript

##############################################################################
#ASSIGN USER INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]

#LOAD R WORKSPACE AND ALL PAKCAGES
load(paste(output_dir, "r_save", sep=""))
library(vegan)
##############################################################################

##############################################################################
#OPEN PNG FILE
png(file=paste(output_dir, "rarefaction_curve.png", sep=""))

#CREATE RAREFACTION PLOT
rarecurve(t(community_merge), col=metadata$color, step=1, ylab="ASVs", label=T)
abline(v=min(rowSums(t(community_merge))))

#CLOSE PNG FILE
dev.off()
##############################################################################

#SAVE R WORKSPACE
save.image(paste(output_dir, "r_save"))

#PRINT TRACKER
print("rarefaction done")