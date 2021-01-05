#!/usr/bin/env Rscript

##############################################################################
#ASSIGN INPUTS TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)

#minimum sum = lowest sum of asvs in all samples to keep asv
minimum_sum <- args[1]
#minimum sample count = lowest number of samples that an asv must appear in to keep asv
minimum_sample_count <- args[2]
#sample save = if minimum sample count not met, sample save (sum of asv across all samples) must be greater than sample save
sample_save <- args[3]
#metadata file containing filenames and sample groups
metadata_file <- args[4]
#output directory specified in original pipeline call
output_dir <- args[5]
##############################################################################

#LOAD IN DADA2 ASV COUNTS
community_merge_unfiltered <- read.table(paste(output_dir, "/seqcounts.txt", sep=""), header=TRUE)


#DETERMINE NUMBER OF ASVS IN BOTH COMMUNITIES
merge_unfiltered_length <- length(rownames(community_merge_unfiltered))

##############################################################################
#METADATA TABLE MANIPULATION


#LOAD IN METADATA FILE
metadata <- read.table(metadata_file, header=FALSE, row.names=1)

#ADD COLUMN NAMES TO METADATA DATAFRAME
colnames(metadata) <- c("sample_group", "sample_name")

#CHANGE SAMPLE GROUP TO FACTORS SO YOU CAN MANIPULATE IT EASIER
metadata$sample_group <- as.factor(metadata$sample_group)
##############################################################################


##############################################################################
#FILTER RAW ASV COUNTS

#BY MINIMUM SUM
community_merge <- community_merge_unfiltered[rowSums(community_merge_unfiltered) > as.integer(minimum_sum),]
#SAVE NUMBER OF ASVS AFTER MINIMUM SUM FILTRATION
merge_filtered_min_sum_length <- length(rownames(community_merge))

#BY NUMBER OF SAMPLES THAT OCCUR IN EACH SAMPLE
community_merge <- community_merge[rowSums(community_merge!=0) > as.integer(minimum_sample_count) | rowSums(community_merge) > as.integer(sample_save),]
#SAVE NUMBER OF ASVS AFTER FILTRATION
merge_filtered_all <- length(rownames(community_merge))

##############################################################################



##############################################################################
#SEPERATING MERGED DATA BY COMMUNITY

#SPLIT METADATA BY COMMUNITY
split_list <- split(metadata, metadata$sample_group)
#DETERMINE WHICH SAMPLES BELONGS TO WHICH COMMUNITY
sample_one_names <- rownames(data.frame(split_list[1]))
sample_two_names <- rownames(data.frame(split_list[2]))

#ADD COLORS TO METDATA DATAFRAME BY COMMUNITY
metadata[sample_one_names, "color"] <- "blue"
metadata[sample_two_names, "color"] <- "red"

#DETERMINE COMMUNITY NAMES BY SAMPLE_GROUP IN METADATA
community_one_name <- levels(metadata$sample_group)[1]
community_two_name <- levels(metadata$sample_group)[2]

#CHANGE HEADERS TO REFLECT SAMPLE NAMES, NOT FILE NAMES
new_col_names <- NULL
for(column_name in colnames(community_merge)){
  if(column_name %in% sample_one_names){
    column_name <- metadata[column_name,]$sample_name
    new_col_names <- append(new_col_names, column_name)
  }
  else{
    column_name <- metadata[column_name,]$sample_name
    new_col_names <- append(new_col_names, column_name)
  }
}

colnames(community_merge) <- new_col_names
colnames(community_merge_unfiltered) <- new_col_names

#SORT COMMUNITY MERGE DATAFRAME IN ORDER BY SAMPLE NAMES
community_merge <- community_merge[, order(names(community_merge))]
community_merge_unfiltered <- community_merge_unfiltered[, order(names(community_merge_unfiltered))]

#ORDER METADATA BY SAMPLE NAMES
metadata <- metadata[order(metadata$sample_name),]

#SPLIT MERGED DATAFRAME INTO COMMUNITY DATAFRAMES BY SAMPLE NAME
split_list <- split(metadata, metadata$sample_group)
sample_one_names <- data.frame(split_list[1])
colnames(sample_one_names) <- c("sample_group", "sample_name", "color")
sample_one_names <- sample_one_names$sample_name
sample_two_names <- data.frame(split_list[2])
colnames(sample_two_names) <- c("sample_group", "sample_name", "color")
sample_two_names <- sample_two_names$sample_name

community_one <- dplyr::select(community_merge, all_of(sample_one_names))
community_two <- dplyr::select(community_merge, all_of(sample_two_names))
community_one_unfiltered <- dplyr::select(community_merge_unfiltered, sample_one_names)
community_two_unfiltered <- dplyr::select(community_merge_unfiltered, sample_two_names)

#IF ANY ROWS CONTAIN ALL 0'S REMOVE IT FROM DATAFRAME
community_one <- community_one[rowSums(community_one) > 0,]
community_two <- community_two[rowSums(community_two) > 0,]
community_one_unfiltered <- community_one_unfiltered[rowSums(community_one_unfiltered) > 0,]
community_two_unfiltered <- community_two_unfiltered[rowSums(community_two_unfiltered) > 0,]

#SAVE LENGTHS OF COMMUNITIES PRE AND POST FILTERING
community_one_pre_length <- length(rownames(community_one_unfiltered))
community_two_pre_length <- length(rownames(community_two_unfiltered))
community_one_post_length <- length(rownames(community_one))
community_two_post_length <- length(rownames(community_two))


#CREATE PRESENCE ABSENCE TRANSFORMED COMMUNITY AND MERGED DATAFRAMES
community_one_presence_absence <- (community_one > 0) * 1
community_two_presence_absence <- (community_two > 0) * 1
community_merge_presence_absence <- (community_merge > 0) * 1
##############################################################################


#WRITE OUT FILTERED ASV COUNTS TO NEW FILES
write.table(community_merge, file=paste(output_dir, "asv_counts_filtered.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(community_one, file=paste(output_dir, community_one_name, "_asv_counts_filtered.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(community_two, file=paste(output_dir, community_two_name, "_asv_counts_filtered.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#SAVE ANY VARIABLES CREATED IN THIS SCRIPT TO AN R WORKSPACE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT PIPELINE TRACKER
print("Data Manipulation and Filtering Done")