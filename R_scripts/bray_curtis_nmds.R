#!/usr/bin/env Rscript

##############################################################################
#LOAD ALL LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES

#LOAD LIBRARIES
library(vegan)
library(ggplot2)

#ASSIGN INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################


##############################################################################
#CREATE AND WRITE BRAY CURTIS DISSIMILARITY MATRICES TO FILES

#COMMUNITY ONE NON-NORMALIZED MATRIX
bray_curtis_community_one <- vegdist(t(community_one), method="bray")
write.table(as.matrix(bray_curtis_community_one), file=paste(output_dir, "/bray_nmds/", "bray_curtis_disimilarity_matrix_", community_one_name, ".txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#COMMUNITY ONE NORMALIZED MATRIX
bray_curtis_community_one_norm <- vegdist(t(community_one_norm), method="bray")
write.table(as.matrix(bray_curtis_community_one_norm), file=paste(output_dir, "/bray_nmds/", "bray_curtis_disimilarity_matrix_", community_one_name, "_norm.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#COMMUNITY ONE PRESENCE ABSENCE MATRIX
bray_curtis_community_one_presence_absence <- vegdist(t(community_one_presence_absence), method="bray")
write.table(as.matrix(bray_curtis_community_one_presence_absence), file=paste(output_dir, "/bray_nmds/", "bray_curtis_disimilarity_matrix_", community_one_name, "_presence_absence.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#COMMUNITY TWO NON-NORMALIZED MATRIX
bray_curtis_community_two <- vegdist(t(community_two), method="bray")
write.table(as.matrix(bray_curtis_community_two), file=paste(output_dir, "/bray_nmds/", "bray_curtis_disimilarity_matrix_", community_two_name, ".txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#COMMUNITY TWO NORMALIZED MATRIX
bray_curtis_community_two_norm <- vegdist(t(community_two_norm), method="bray")
write.table(as.matrix(bray_curtis_community_two_norm), file=paste(output_dir, "/bray_nmds/", "bray_curtis_disimilarity_matrix_", community_two_name, "_norm.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#COMMUNITY TWO PRESENCE ABSENCE MATRIX
bray_curtis_community_two_presence_absence <- vegdist(t(community_two_presence_absence), method="bray")
write.table(as.matrix(bray_curtis_community_two_presence_absence), file=paste(output_dir, "/bray_nmds/", "bray_curtis_disimilarity_matrix_", community_two_name, "_presence_absence.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#BOTH COMMUNITIES NON-NORMALIZED MATRIX
bray_curtis_community_merge <- vegdist(t(community_merge), method="bray")
write.table(as.matrix(bray_curtis_community_merge), file=paste(output_dir, "/bray_nmds/",  "bray_curtis_disimilarity_matrix_all_samples.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#BOTH COMMUNITIES NORMALIZED MATRIX
bray_curtis_community_merge_norm <- vegdist(t(community_merge_norm), method="bray")
write.table(as.matrix(bray_curtis_community_merge_norm), file=paste(output_dir, "/bray_nmds/",  "bray_curtis_disimilarity_matrix_all_samples_norm.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#BOTH COMMUNITIES PRESENCE ABSENCE MATRIX
bray_curtis_community_merge_presence_absence <- vegdist(t(community_merge_presence_absence), method="bray")
write.table(as.matrix(bray_curtis_community_merge_presence_absence), file=paste(output_dir, "/bray_nmds/",  "bray_curtis_disimilarity_matrix_all_samples_presence_absence.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

##############################################################################

##############################################################################
#CONVERSION FROM BRAY CURTIS MATRIX TO NMDS COORDINATES
NMDS_community_merge <- metaMDS(bray_curtis_community_merge)
NMDS_community_merge_norm <- metaMDS(bray_curtis_community_merge_norm)
NMDS_community_merge_pres_absence <- metaMDS(bray_curtis_community_merge_presence_absence)

NMDS_community_one <- metaMDS(bray_curtis_community_one)
NMDS_community_one_norm <- metaMDS(bray_curtis_community_one_norm)
NMDS_community_one_presence_absence <- metaMDS(bray_curtis_community_one_presence_absence)

NMDS_community_two <- metaMDS(bray_curtis_community_two)
NMDS_community_two_norm <- metaMDS(bray_curtis_community_two_norm)
NMDS_community_two_presence_absence <- metaMDS(bray_curtis_community_two_presence_absence)
##############################################################################

##############################################################################
#NMDS PLOTTER
NMDS_plotter <- function(NMDS, title){
  #EXTRACT COORDINATES FROM metaMDS OBJECT
  nmds1 <- NMDS$points[,1]
  nmds2 <- NMDS$points[,2]
  colors <- NULL
  
  #CREATE LIST OF COLORS
  for(sample in rownames(NMDS$points)){
    colors <- append(colors, metadata$color[sample==metadata$sample_name])
  }
  
  #CREATE TEMPORARY DATAFRAME CONTAINING COORDINATES
  temp_frame <- data.frame(cbind(nmds1, nmds2))

  #EXTRACT STRESS FROM metaMDS OBJECT AND ASSIGN TO VARIABLE
  if(NMDS$stress <= 0.001){
    stress <- "<= 0.001"
  }
  else{
    stress <- paste("=", round(NMDS$stress, 3), sep= " ")
  }

  #PLOT NMDS
  plot <- ggplot(data=temp_frame, aes(nmds1, nmds2)) +
    #stat_ellipse(type="t", size=0.1) +
    ggtitle(paste(title, " / stress ", stress, sep="")) + xlab("NMDS1") + ylab("NMDS2") +
    geom_text(color=colors, label=rownames(NMDS$points))

  return(plot)
}

#BOTH COMMUNITIES NORMALIZED
plot <- NMDS_plotter(NMDS_community_merge_norm, "All Samples Normalized NMDS Plot")
ggsave(filename="all_samples_normalized_nmds.png", plot=plot, device="png", path=paste(output_dir, "/bray_nmds/", sep=""))

#COMMUNITY ONE NORMALIZED
plot <- NMDS_plotter(NMDS_community_one_norm, paste(community_one_name, " Normalized NMDS Plot", sep=""))
ggsave(filename=paste(community_one_name, "_normalized_nmds.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/bray_nmds/", sep=""))

#COMMUNITY TWO NORMALIZED
plot <- NMDS_plotter(NMDS_community_two_norm, paste(community_two_name, " Normalized NMDS Plot", sep=""))
ggsave(filename=paste(community_two_name, "_normalized_nmds.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/bray_nmds/", sep=""))

#BOTH COMMUNITIES NON-NORMALIZED
plot <- NMDS_plotter(NMDS_community_merge, "All Samples NMDS Plot")
ggsave(filename="all_samples_nmds.png", plot=plot, device="png", path=paste(output_dir, "/bray_nmds/", sep=""))

#COMMUNITY ONE NON-NORMALIZED
plot <- NMDS_plotter(NMDS_community_one, paste(community_one_name, " NMDS Plot", sep=""))
ggsave(filename=paste(community_one_name, "_nmds.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/bray_nmds/", sep=""))

#COMMUNITY TWO NON-NORMALIZED
plot <- NMDS_plotter(NMDS_community_two, paste(community_two_name, " NMDS Plot", sep=""))
ggsave(filename=paste(community_two_name, "_nmds.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/bray_nmds/", sep=""))

##############################################################################


#SAVE R WORKSPACE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT TRACKER
print("Bray Curtis/NMDS Done")