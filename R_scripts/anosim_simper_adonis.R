#!/usr/bin/env Rscript

##############################################################################
#LOAD LIBRARIES
library(vegan)

#ASSIGN INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
low_perc <- args[2]
cum_sum <- args[3]

#LOAD WORKSPACE
load(paste(output_dir, "/r_save", sep=""))
##############################################################################


##############################################################################
#ANOSIM

#MAKE SURE SAMPLE_GROUP IS A FACTOR
metadata$sample_group <- as.factor(metadata$sample_group)

#ANOSIM BETWEEN COMMUNITIES BY SAMPLE_GROUP
anosim_between_communities <- anosim(t(community_merge_norm), metadata$sample_group)

#OPEN ANOSIM FILE
sink(file=paste(output_dir, "/anosim_simper/anosim_results.txt", sep=""))
#PRINT ANOSIM RESULTS
summary(anosim_between_communities)
#CLOSE ANOSIM FILES
sink()

#OPEN ANOSIM PLOT FILE
png(paste(output_dir, "/anosim_simper/anosim_between_communities.png", sep=""))
#AUTO PLOT
plot(anosim_between_communities)
#CLOSE ANOSIM PLOT FILE
while (!is.null(dev.list()))  dev.off()
##############################################################################



## The SIMPER Test and Understanding it's Results
#The Simper test compares the two communitites and tells you which ASV's (taxonomy) contribute the most to the bray-curtis dissimilarity index.  
#Cumsum will always add up to one (bray-curtis index) and average is average contribution of rowname-taxonomy. 
#However, in order to determine if the communities are significantly different the above anosim must be referenced.

#We see a number of OTUs that may differ 
#However, these are just the OTUs that most contribute to Bray-Curtis measures between our age groups.
#They are not necessarily significantly different, Kruskal test required.

#kruskal transform community data

##############################################################################
#SIMPER

#BOTH COMMUNITES SIMPER ON NORMALIZED DATA - CREATES UNIQUE SIMPER OBJECT
community_merge_normalized_simper <- simper(t(community_merge_norm), metadata$sample_group)

#WRITE OUT SUMMARY RESULTS
simp_table_name <- paste(community_two_name, "_", community_one_name, sep="")
community_merge_summary <- summary(community_merge_normalized_simper)
write.table(data.frame(community_merge_summary[1]), file=paste(output_dir, "/anosim_simper/", "simper_summary_all_ASVs_", simp_table_name, ".txt", sep=""), append=FALSE, sep="\t")

#BOTH COMMUNITIES SIMPER PRESENCE ABSENCE
community_merge_presence_absence_simper <- simper(t(community_merge_presence_absence), metadata$sample_group)

#WIRTE OUT SUMMARY RESULTS
community_pres_abs <- summary(community_merge_presence_absence_simper)
write.table(data.frame(community_pres_abs[1]), file=paste(output_dir, "/anosim_simper/", "simper_summary_all_ASVs_presence_absence_", simp_table_name, ".txt", sep=""), append=FALSE, sep="\t")
##############################################################################

##############################################################################
#SIMPER FILTER FUNCTION
#TAKES SIMPER DATAFRAME AND FILTERS ASVS WITH LOW CONTRIBUTION OUT BASED ON LOW_PERCENT AND HIGH_CUMULATIVE
simper_filter <- function(simper, low_perc, high_cum, rm){
  #SIMPER MUST BE DATA FRAME, LOW_PERCENT IS THE MINIMUM PERCENTAGE CONTRIBUTION AN asv MUST HAVE TO BE KEPT
  #SIMPER IS SORTED FROM LARGEST TO LEAST CONTRIBUTION, HIGH_CUM MARKS THE CUMULATIVE CONTRIBUTION CUTOFF
  #RM SHOULD BE FALSE IF SIMPER DATAFRAME IS FROM BOTH COMMUNITIES, TRUE IF IT IS FROM ONLY ONE COMMUNITY
  if(rm==TRUE){
    #REMOVE BETWEEN COMMUNITY DIFFERENCES (AVA AND AVB)
    colnames(simper) <- column_names
    simper <- dplyr::select(simper, -c('sd','ratio','ava', 'avb'))
  }
  #FILTER SIMPER RESULTS AND REMOVE ALL ASV'S BELOW LOW PERCENTAGE
  simper <- simper[(simper$average > low_perc),]

  #FILTER SIMPER RESULTS AND REMOVE ALL ASV'S ABOVE HIGH_CUM
  simper <- simper[(simper$cumsum < high_cum),]

  #CREATE A LIST OF ALL ASV'S
  species <- rownames(simper)
  sig_list <- NULL
  for(i in species){
    #FOR EACH ASV
    #CHECK FOR SIGNIFICANT DIFFERENE BETWEEN COMMUNITIES AND APPEND TO SIG_LIST
    sig_list <- append(sig_list, kruskal.test(data.frame(t(community_merge))[,i]~metadata$sample_group)$p.value)
  }
  #hOCHBERG CORRECT SIGNIFICANCE VALUES
  sig_list_adjusted <- p.adjust(sig_list, method="hochberg")

  #ADD SIGNIFICANCE VALUES TO SIMPER TABLE
  simper$Significance <- sig_list
  simper$Adjusted_Significance <- sig_list_adjusted

  #RETURN MODIFIED SIMPER LIST
  return(simper)
}
#############################################################################

#############################################################################
#DETERMINE WHAT COLUMN NAMES SHOULD BE
column_names <- c("average", "sd","ratio","ava","avb", "cumsum")
#EXTRACT DATA FROM SIMPER OBJECT SUMMARY
community_merge_summary_data <- data.frame(community_merge_summary[1])
#ASSIGN NEW COLUMN NAMES
colnames(community_merge_summary_data) <- column_names
#CALL SIMPER FILTER FUNCTION AND WRITE TO TABLE
community_merge_normalized_simper_analysis <- simper_filter(community_merge_summary_data, low_perc, cum_sum, FALSE)
write.table(data.frame(community_merge_normalized_simper_analysis), file=paste(output_dir, "/anosim_simper/", "simper_summary_filtered_ASVs_both_communities.txt", sep=""), append=FALSE, sep="\t")

#REPEAT PREVIOUS FOR PRESENCE ABSENCE DATA
community_merge_summary_presence_absence_data <- data.frame(community_pres_abs[1])
colnames(community_merge_summary_presence_absence_data) <- column_names
community_merge_presence_absence_simper_analysis <- simper_filter(community_merge_summary_presence_absence_data, low_perc, cum_sum, FALSE)
write.table(data.frame(community_merge_presence_absence_simper_analysis), file=paste(output_dir, "/anosim_simper/", "simper_summary_filtered_ASVs_presence_absence_both_communities.txt", sep=""), append=FALSE, sep="\t")

#average = species contribution to between-group disimilarity
#sd = standard deviation of contribution
#ratio = average to standard deviation
#ava,avb average abundance per group
#cumsum = cumulative sum of differences
#significance = kruskal test significance

#SIMPER FOR COMMUNITY ONE - CREATES SIMPER OBJECT FOR EVERY COMPARISON BETWEEN SAMPLES
community_one_normalized_simper <- simper(t(community_one_norm), colnames(community_one))

#EXTRACT SUMMARY FROM SIMPER OBJECT
community_one_normalized_simper_summary <- summary(community_one_normalized_simper)


#FOR EACH COMPARISON IN SIMPER SUMMARY RUN THROUGH SIMPER SUMMARY AND WRITE TO FILES
for(i in 1:length(community_one_normalized_simper_summary)){
  name <- names(community_one_normalized_simper_summary[i])
  simper <- simper_filter(data.frame(community_one_normalized_simper_summary[1]), low_perc, cum_sum, TRUE)
  write.table(simper, file=paste(output_dir, "/anosim_simper/community_one/", name, ".txt", sep=""), append=FALSE, sep="\t")
}


#REPEAT FOR COMMUNITY TWO
community_two_normalized_simper <- simper(t(community_two_norm), colnames(community_two))

community_two_normalized_simper_summary <- summary(community_two_normalized_simper)

for(i in 1:length(community_two_normalized_simper_summary)){
  name <- names(community_two_normalized_simper_summary[i])
  simper <- simper_filter(data.frame(community_two_normalized_simper_summary[1]), low_perc, cum_sum, TRUE)
  write.table(simper, file=paste(output_dir, "/anosim_simper/community_two/", name, ".txt", sep=""), append=FALSE, sep="\t")
}

##############################################################################

##############################################################################
#ADONIS

#BOTH COMMUNITIES NORMALIZED
sink(file=paste(output_dir, "/adonis/", "community_merge.txt",sep=""))
adonis(bray_curtis_community_merge_norm ~ metadata$sample_group)
sink()

#BOTH COMMUNITIES PRESENCE ABSENCE
sink(file=paste(output_dir, "/adonis/", "community_merge_presence_absence.txt",sep=""))
adonis(bray_curtis_community_merge_presence_absence ~ metadata$sample_group)
sink()

#COMMUNITY ONE NORMALIZED
sink(file=paste(output_dir, "/adonis/", "community_one.txt",sep=""))
adonis(bray_curtis_community_one_norm ~ colnames(community_one))
sink()

#COMMUNITY ONE PRESENCE ABSENCE
sink(file=paste(output_dir, "/adonis/", "community_one_presence_absence.txt",sep=""))
adonis(bray_curtis_community_one_presence_absence ~ colnames(community_one))
sink()

#COMMUNITY TWO NORMALIZED
sink(file=paste(output_dir, "/adonis/", "community_two.txt", sep=""))
adonis(bray_curtis_community_two_norm ~ colnames(community_two))
sink()

#COMMUNITY TWO PRESENCE ABSENCE
sink(file=paste(output_dir, "/adonis/", "community_two_presence_absence.txt", sep=""))
adonis(bray_curtis_community_two_presence_absence ~ colnames(community_two))
sink()
##############################################################################

#SAVE R WORKSAPCE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT TRACKER
print("Anosim/Simper/Adonis Done")