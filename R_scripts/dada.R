#!/usr/bin/env Rscript

library(dada2)
library(sys)
library(ggplot2)

#### parsing user arguments from master bash script

args = commandArgs(trailingOnly=TRUE)

output_dir <- args[2]

fh.path <- args[1]

setwd(fh.path)

# capturing some useful output to a summary file
summ_file <- paste(output_dir, "/dada_summary.txt", sep = "")
sink(file = summ_file, append = TRUE)

print("script start")
Sys.time()

#### initial file management 

# checking whether input files are in sequencer format or SRA format
# sequencer: SAMPLENAME_R1_001.fastq, SAMPLENAME_R2_001.fastq
# SRA: SAMPLENAME_1.fastq, SAMPLENAME_2.fastq
# then splitting into forward and reverse reads

file_check <- list.files(fh.path, pattern = "_R1_001.fastq", full.names = FALSE)

if(length(file_check) != 0){

  fw.reads <- sort(list.files(fh.path, pattern = "_R1_001.fastq", full.names = TRUE))

  rv.reads <- sort(list.files(fh.path, pattern = "_R2_001.fastq", full.names = TRUE))

} else {

  fw.reads <- sort(list.files(fh.path, pattern = "_1.fastq", full.names = TRUE))

  rv.reads <- sort(list.files(fh.path, pattern = "_2.fastq", full.names = TRUE))

}

# extracting sample names: works for both sequencer and SRA format
# basename ditches the path, then separates the name by underscore, then selects first object
# which is the sample name, iterated for each file in fw.reads

sample.names <- sapply(strsplit(basename(fw.reads), "_"), `[`, 1)

#### plotting the quality profile for all unfiltered/untrimmed files

# saving each plot independently as a .png in quality plots subdirectory initialized in bash script

quality_plots_dir <- paste(output_dir, "/quality_plots", sep = "")

for(i in 1:length(fw.reads)){
  forward_quality <- plotQualityProfile(fw.reads[i])
  current_name <- paste(sample.names[i], "_fw_quality.png", sep = "")
  ggsave(filename=current_name, plot=forward_quality, device="png", path=quality_plots_dir)
}

for(i in 1:length(rv.reads)){
  reverse_quality <- plotQualityProfile(rv.reads[i])
  current_name <- paste(sample.names[i], "_rv_quality.png", sep = "")
  ggsave(filename=current_name, plot=reverse_quality, device="png", path=quality_plots_dir)
}

#### Quality filtering and trimming

# creating a new subdirectory to store filtered files
# keeping naming consistent with SRA fw/rv convention
# named to be compressed upon creation

fw.filtered <- file.path(output_dir, "filtered", paste(sample.names, "_1_filt.fastq.gz", sep = ""))

rv.filtered <- file.path(output_dir, "filtered", paste(sample.names, "_2_filt.fastq.gz", sep = ""))

# naming the objects here allows the sample names to be carried over into the sequence table

names(fw.filtered) <- sample.names

names(rv.filtered) <- sample.names

print("Filtering")
Sys.time()

# paramters for the filtering step should be tweaked depending on the data, but this is a good starting point

filter.out <- filterAndTrim(fw.reads, fw.filtered, rv.reads, rv.filtered, trimLeft = c(19,20),
                            truncQ = 11, maxN = 0, maxEE = 2, rm.phix = TRUE,
                            compress = TRUE, multithread = TRUE)

# printing the filter object to the summary file
filter.out

#### Learning error rates 

# this will often be the rate-limiting step so don't be surprised if the program hangs up here for a while

# Plots of the error rates for each sample can be easily generated if desired, see dada tutorial

print("Error rates")
Sys.time()

fw.error <- learnErrors(fw.filtered, multithread = TRUE)

rv.error <- learnErrors(rv.filtered, multithread = TRUE)

#### Sample inference... THE COOL PART!

print("Sample inference")
Sys.time()

fw.dada <- dada(fw.filtered, err = fw.error, multithread = TRUE)

rv.dada <- dada(rv.filtered, err = rv.error, multithread = TRUE)

# saving some useful summaries to be outputted later

fw.dada.summary <- fw.dada[[1]]

rv.dada.summary <- rv.dada[[1]]

fw.dada.summary

rv.dada.summary

#### Merge read pairs 

print("Merging pairs")
Sys.time()

merged.reads <- mergePairs(fw.dada, fw.filtered, rv.dada, rv.filtered)

#### Creating an ASV table

seq.table <- makeSequenceTable(merged.reads)

seq.dims <- dim(seq.table)

# checking the output here will tell us if the seqs fall within
# a reasonable length for the V4 region

table.check <- table(nchar(getSequences(seq.table)))

#### Removing chimeras

seq.table.nochimeras <- removeBimeraDenovo(seq.table, method = "consensus", 
                                           multithread = TRUE)

noch.dims <- dim(seq.table.nochimeras)

# dividing chimeras by seq abundance can tell us percent of non chimeric seqs

perc.nonchim <- sum(seq.table.nochimeras)/sum(seq.table)

print("Sequence table dimensions")
seq.dims

print("Sequence length check")
table.check

print("sequence table dimensions, no chimeras")
noch.dims

print("proportion of non-chimeric sequences")
perc.nonchim

Sys.time()

#### Checking reads that made it through dada processing

# a function to count unique reads in an object

count.fun <- function(reads_file){
  sum(getUniques(reads_file))
}

# summarizing the number of reads that passed each stage

track.reads <- cbind(filter.out, sapply(fw.dada, count.fun),
                     sapply(rv.dada, count.fun), sapply(merged.reads, count.fun), rowSums(seq.table.nochimeras))

colnames(track.reads) <- c("input", "filtered", "denoised.fw", "denoised.rv", "merged", "nonchimera")

rownames(track.reads) <- sample.names

print("Tracking reads through Dada2 processing")

track.reads

#### Taxonomic assignment 

silva_ref <- "./silva_nr99_v138_train_set.fa.gz"

species_ref <- "./silva_species_assignment_v138.fa.gz"

taxa <- assignTaxonomy(seq.table.nochimeras, silva_ref, multithread = TRUE)

taxa <- addSpecies(taxa, species_ref)

# transposing the sequence table so ASVs are rows and samples are columns

seq.table.nochimeras <- t(seq.table.nochimeras)

#### Naming output dataframes with shorthand ASV names rather than the sequences

make_ASVnames <- function(n){
  # n is the number of rows in seqtable
  # function takes n, and creates n-number asv labels
  asv_list <- vector()
  for(num in 1:n){
    current_asv <- paste("ASV", toString(num), sep = "")
    asv_list <- append(asv_list, current_asv)
  }
  return(asv_list)
}

# creating a df with the shorthand ASV name and the full sequence

asv.seqs <- data.frame(row.names = make_ASVnames(nrow(seq.table.nochimeras)),
                       "Sequence" = row.names(seq.table.nochimeras))

# reassigning the rownames of the seqtable and taxa table with shorthand ASV names

row.names(seq.table.nochimeras) <- row.names(asv.seqs)

row.names(taxa) <- row.names(asv.seqs)

#### Output files

# taxonomic assignment table
write.table(taxa, file = paste(output_dir, "/taxmat_samples.txt", sep=""), sep = "\t", row.names = TRUE, col.names = TRUE)

# sequence table
write.table(seq.table.nochimeras, file = paste(output_dir, "/seqcounts.txt", sep=""), sep = "\t", row.names = TRUE, col.names = TRUE)

# ASV names and sequences in fasta format
for(i in 1:nrow(asv.seqs)){
  cat(paste(">", row.names(asv.seqs)[i], "\n", sep = ""), file = paste(output_dir, "asv_sequences.fasta", sep=""), append = TRUE)
  cat(paste(asv.seqs$Sequence[i], "\n", sep = ""), file = paste(output_dir, "/asv_sequences.fasta", sep=""), append = TRUE)
}

# ASV names and sequences in tabular format
write.table(asv.seqs, file = paste(output_dir, "/asv_sequences.txt", sep=""), sep = "\t", row.names = TRUE, col.names = TRUE)

print("Files written out, end script")
Sys.time()

# closing the summary file (rerouting output back to standard out)
sink()

















