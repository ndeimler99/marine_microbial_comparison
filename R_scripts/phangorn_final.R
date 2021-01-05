#!/usr/bin/env Rscript

# LOADING PACKAGES

library("Biostrings")
library("msa")
library(ape)
library(phangorn)
library(phyloseq)


#####################################################
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

load(paste(output_dir, "r_save", sep=""))

# LOADING R ENVIRONMENT


#########################################
seqs <- read.table(paste(output_dir, "asv_sequences.txt", sep=""))

# <<<<<<<<<<<<<<<<<<<<< For taxonomic table, work from this <<<<<<<<<<<<<<<<<<<

# Does this command make a tax table?
#####

present_one <- ceiling(ncol(community_one) / 4)


#########################################################################
# PREPARING FILTERED DATA SETS TO WORK WITH
# 2 CRITERIA TO INCLUDE:
# 100 COUNTS, AVERAGE PER COLUMN
# MUST BE PRESENT IN AT LEAST 25% OF THE SAMPLES (ROUNDED UP)

# DEFINING 25% PRESENT & INITIALIZING LIST
present_one <- ceiling(ncol(community_one) / 4)
community_one_list <- c()

# SELECTING ASVs THAT MEET CRITERIA FOR COMMUNITY 1 WITH LOOP THROUGH PREVIOUSLY FILTERED SET
for (row in 1:nrow(community_one)){
  if ((sum(community_one[row,]) / ncol(community_one)) >=100) {
    if (rowSums(community_one[row,] != 0) >= present_one){
          community_one_list <- c(community_one_list, rownames(community_one[row,]))
    }
  
  }
}

# SETTING VARIABLE FOR FILTERED ASVs MEETING CRITERIA
comm_one_filter <- seqs[seqs$ASV %in% community_one_list,]


#########################################################################
# REPEATING ABOVE FOR COMMUNITY 2
present_two <- ceiling(ncol(community_two) / 4)
community_two_list <- c()

for (row in 1:nrow(community_two)){
  if ((sum(community_two[row,]) / ncol(community_two)) >=100) {
    if (rowSums(community_two[row,] != 0) >= present_two){
      community_two_list <- c(community_two_list, rownames(community_two[row,]))
    }
    
  }
}


comm_two_filter <- seqs[seqs$ASV %in% community_two_list,]

#########################################################################
# MAKING DATAFRAME INCLUDING ALL SEQS MEETING CRITERIA
# FROM BOTH COMMUNITIES, COMBINED

merged_filter <- seqs[seqs$ASV %in% c(community_one_list, community_two_list),]

#########################################################################
# RETRIEVING TAXONOMY NAMES

# TAXONOMY TABLE
tax <- tax_table(community_merge_phy_obj)

# REMOVING ASVs NOT IN COMMUNITY MERGED_FILTER TABLE
# EVERY ASV IN ALL TREES ARE IN THIS TABLE
tax_selected <- tax[rownames(community_merge),]


# ADDS NEW COLUMN TO tax_selected, "Final", WITH PLACEHOLDER NAs
tax_selected$Final <- NA



# IN LOOP:
# SETS COLUMN "Final" EQUAL TO MOST SPECIFIC TAXONOMIC CLASSIFICATION (LAST ONE NOT "NA")
# IF NO TAXONOMY LEVELS ARE KNOWN (ALL "NA") - SETS "Final" COLUMN TO "?"

tax_col_num <- ncol(tax_selected)

for (row in 1:nrow(tax_selected) ) {
  current <- tax_selected[row,]
  na_num <- apply(tax_selected[row,], MARGIN = 1, function(x) sum(is.na(x)))
  
  if ( (tax_col_num - na_num) != 0) {
  
    tax_selected[row,]$Final <- (current[,(tax_col_num - na_num)])
  } else {
    tax_selected[row,]$Final <- "?"
  }
  
}

#########################################################################
# OUTPUT FILTERED ASVs FOR EACH COMMUNITY (AND SHARED) TO 3 FASTA FILES

# <<<<<<<<<<<<<< MAY NEED TO CHANGE THE OUTPUT DIR - GET NATHAN FOR HELP >>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<< The output directory for these FASTA files needs to be accessible

output_dir <- getwd()


for(i in 1:nrow(comm_one_filter)){
  cat(paste(">", row.names(comm_one_filter)[i], "\n", sep = ""), file = paste(output_dir, "/comm_one_filtered_sequences.fasta", sep=""), append = TRUE)
  cat(paste(comm_one_filter$Sequence[i], "\n", sep = ""), file = paste(output_dir, "/comm_one_filtered_sequences.fasta", sep=""), append = TRUE)
}

for(i in 1:nrow(comm_two_filter)){
  cat(paste(">", row.names(comm_two_filter)[i], "\n", sep = ""), file = paste(output_dir, "/comm_two_filtered_sequences.fasta", sep=""), append = TRUE)
  cat(paste(comm_two_filter$Sequence[i], "\n", sep = ""), file = paste(output_dir, "/comm_two_filtered_sequences.fasta", sep=""), append = TRUE)
}

for(i in 1:nrow(merged_filter)){
  cat(paste(">", row.names(merged_filter)[i], "\n", sep = ""), file = paste(output_dir, "/comm_merged_filtered_sequences.fasta", sep=""), append = TRUE)
  cat(paste(merged_filter$Sequence[i], "\n", sep = ""), file = paste(output_dir, "/comm_merged_filtered_sequences.fasta", sep=""), append = TRUE)
}

#########################################################################

# AT THIS POINT:
# FILTERED SEQS ARE TURNED INTO FASTA FILES - ALL OUTPUT TO WORKING DIRECTORY
# IN tax_selected: ALL POST-FILTERED ASVs HAVE ASSOCIATED "FINAL" TAXONOMIES, FOR LATER USE IN TREES

#########################################################################

#########################################################################
# ALIGNING SEQUENCES, FOR USE IN DISTANCE MATRICES; TREE BUILDING

# READING IN FASTA FILES AS DNAStringSets
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# <<<< NEED TO CHECK THAT THIS OUTPUT_DIR PORTION WILL WORK WITH THE MAIN SCRIPT and current setup for directories
orf_comm_one <- readDNAStringSet( (paste(output_dir, "/comm_one_filtered_sequences.fasta", sep="")), format = "fasta" )
orf_comm_two <- readDNAStringSet( (paste(output_dir, "/comm_two_filtered_sequences.fasta", sep="")), format = "fasta" )
orf_comm_merged <- readDNAStringSet( (paste(output_dir, "/comm_merged_filtered_sequences.fasta", sep="")), format = "fasta" )

# MAKING ALIGNMENT OBJECTS WITH MSA - MULTIPLE SEQUENCE ALIGNMENT
align_c_one <- msa(orf_comm_one)
align_c_two <- msa(orf_comm_two)
align_c_merged <- msa(orf_comm_merged)

#########################################################################
# Tree building


# CONVERTING ALIGNMENTS INTO phyDat OBJECTS
c_one_phydat <- msaConvert(align_c_one, "phangorn::phyDat")
c_two_phydat <- msaConvert(align_c_two, "phangorn::phyDat")
c_merged_phydat <- msaConvert(align_c_merged, "phangorn::phyDat")


#########################################################################
# BUILDING COMMUNITY 1 TREE (PYROSOMES)
dmat_pyr <- dist.ml(c_one_phydat)
pyr_nj <- nj(dmat_pyr)


# RE-ROOT OUTPUT OF NEIGHBOR-JOINING TREE METHOD
# WITH LAST TIP IN $tip.label VECTOR (OUTLIER OF GROUP)
pyr_last_tip <- pyr_nj$tip.label[length(pyr_nj$tip.label)]

# SETTING NUMBER OF TIPS VARIABLE
pyr_tip_num <- length(pyr_nj$tip.label)

reroot_pyr <- root(pyr_nj, pyr_last_tip, resolve.root=TRUE)


# <<<<<<<<<<<<< After making sure this works hard-coded
# <<<<<<<<<<<<<<<< change all the variables to these

#########################################################################
# SETTING VARIABLES FOR BOOTSTRAPPING SUPPORT & OPTIMIZATION

# NUMBER OF ITERATIONS
num_bootstrap <- 500

# SELECTS CUT-OFF PERCENTAGE FOR SHOWING JUNCTION VALUES
bs_support <- 70

# TREE REARRANGEMENT CHOICE OPTIONS: TBR (DEFAULT), NNI, SPR
rearrangement_choice <- 'tbr'

# USING NUMBER OF CORES DETECTED WHERE MULTICORE OPTION IS ENABLED
cores <- detectCores()

#########################################################################
# BOOTSTRAPPING AND OPTIMIZATION


# OPTIMIZING WITH SELECTED METHOD
pml_pyr <- pml(reroot_pyr, c_one_phydat)
fit_pyr_tbr <- optim.pml(pml_pyr, model='JC', rearrangement='tbr')

# BOOTSTRAPPING
bs_pyr_tbr <- bootstrap.pml(fit_pyr_tbr, bs=100, optNni=TRUE, control=pml.control(trace=0),
        multicore=TRUE, mc.cores=cores)


# PLOTTING TREE WITH BS SUPPORT VALUES AT JUNCTIONS
bs_pyr_title <- (paste("Bootstrap Support Values - 
Pyrosome Community - Root: ", pyr_last_tip))
plotBS(root((fit_pyr_tbr$tree), pyr_last_tip), bs_pyr_tbr, p=bs_support, type='p', 
       main=bs_pyr_title, cex=0.4, show.tip.label=FALSE)


# MAKING TBR OPTIMIZED TREE OUTGROUP ROOTED
pyr_optim_tree <- root((fit_pyr_tbr$tree), pyr_last_tip, resolve.root=TRUE)


######## ORIGINAL TREE (NON-AGGLOMERATED) WITH LABELS
# PARSIMONY MEASURE IS INCLUDED WITH ALL TREES
# AS TREES ARE TIP-AGGLOMERATED, THE LEVEL OF COMPLEXITY DROPS
# PARSIMONY QUANTIFIES THIS DROP IN COMPLEXITY WITH INCREMENTAL AGGLOMERATION

# PARSIMONY AS MEASURE OF COMPLEXITY
pyr_orig_parsimony <- parsimony(pyr_optim_tree, c_one_phydat)

# SUB-TITLE TO INCLUDE PARSIMONY AND LABEL
sub_title_orig_pyr <- (paste("Parsimony: ", pyr_orig_parsimony))

# PLOTTING ORIGINAL TREE
plot(pyr_optim_tree, cex=0.4, main='Pyrosome - Community 1
Full Tree', sub=sub_title_orig_pyr, show.tip.label=TRUE)

######## ORIGINAL TREE WITHOUT LABELS - FOR EASY VISUALIZATION OF UNDERLYING STRUCTURE
plot(pyr_optim_tree, main='Pyrosomes - Community 1
Full Tree', sub=sub_title_orig_pyr, show.tip.label=FALSE)


#########################################################
######## PLOTTING TREE AFTER TIP AGGLOMERATION (h = 0.15)

# TIP AGGLOMERATION SET TO 0.15
pyr_glom_15 <- tip_glom(pyr_optim_tree, h=0.15)

# GETTING PARSIMONY AS MEASURE OF COMPLEXITY, POST-AGGLOMERATION
pyr_g15_parsimony <- parsimony(pyr_glom_15, c_one_phydat)

# ADDING PARSIMONY TO SUB-TITLE
sub_title_15_pyr <- (paste("Parsimony: ", pyr_g15_parsimony))

# PLOTTING TREE (FIRST LEVEL AGGLOMERATION)
plot(pyr_glom_15, cex=0.4, main='Pyrosomes - Community 1
Tip Agglomerated
(Partial)', sub=sub_title_15_pyr)

#########################################################
######## PLOTTING TREE AFTER TIP AGGLOMERATION (h = 0.30)

# TIP AGGLOMERATION SET TO 0.30
pyr_glom_30 <- tip_glom(pyr_optim_tree, h=0.3)

# GETTING PARSIMONY AS MEASURE OF COMPLEXITY, POST-AGGLOMERATION
pyr_g30_parsimony <- parsimony(pyr_glom_30, c_one_phydat)
sub_title_30_pyr <- paste("Parsimony: ", pyr_g30_parsimony)
plot(pyr_glom_30, cex=0.7, main='Pyrosomes - Community 1
Tip Agglomerated
(Complete)', sub=sub_title_30_pyr)

#######################
# PLOTTING TREE WITH TAXON LABELS INSTEAD OF ASVs, FOR FINAL AGGLOMERATION (h = 0.30)

# COLLECTING ORIGINAL TIP LABELS (ASVs) IN VECTOR - ORDERED
pyr_tip_names <- pyr_glom_30$tip.label

# CONVERTING TIP LABELS TO THE MOST SPECIFIC TAXONOMIC LEVEL (FROM tax_selected)
pyr_glom_30$tip.label <- tax_selected[pyr_tip_names,]$Final

# PLOTTING TREE WITH TAXONOMIC LABELS
plot(pyr_glom_30, cex=0.7, main='Pyrosomes - Community 1
Tip Agglomerated with Taxonomies Included
(Complete)', sub=sub_title_30_pyr)

#########################################################################

#########################################################################
# SAME PROCESS FOR COMMUNITY 2 - SEAWATER


# BUILDING TREE
dmat_sea <- dist.ml(c_two_phydat)
sea_nj <- nj(dmat_sea)

# RE-ROOTING
sea_last_tip <- sea_nj$tip.label[length(sea_nj$tip.label)]
sea_tip_num <- length(sea_nj$tip.label)

reroot_sea <- root(sea_nj, sea_last_tip, resolve.root=TRUE)

# BS AND PARSIMONY
pml_sea <- pml(reroot_sea, c_two_phydat)

# OPTIMIZING
fit_sea_tbr <- optim.pml(pml_sea, model='JC', rearrangement='tbr')

# BOOTSTRAPPING
bs_sea_tbr <- bootstrap.pml(fit_sea_tbr, bs=100, optNni=TRUE, control=pml.control(trace=0),
                            multicore=TRUE, mc.cores=cores)

# PLOTTING BS TREE
bs_sea_title <- (paste("Bootstrap Support Values - 
Seawater Community - Root: ", sea_last_tip))
plotBS(root(fit_sea_tbr$tree, sea_last_tip), bs_sea_tbr, p=bs_support, type='p', 
       main=bs_sea_title, cex=0.4, show.tip.label=FALSE)

# RE-ROOTING
sea_optim_tree <- root(fit_sea_tbr$tree, sea_last_tip, resolve.root=TRUE)

# PLOTTING ORIGINAL WITH LABELS (AND PARSIMONY)
sea_orig_parsimony <- parsimony(sea_optim_tree, c_two_phydat)
sub_title_orig_sea <- (paste("Parsimony: ", sea_orig_parsimony))
plot(sea_optim_tree, cex=0.4, main='Seawater - Community 2
Full Tree', sub=sub_title_orig_sea)
# WITHOUT LABELS
plot(sea_optim_tree, cex=0.4, main='Seawater - Community 2
Full Tree', sub=sub_title_orig_sea, show.tip.label=FALSE)

# TIP AGGLOMERATION TO h = 0.15 & PARSIMONY
sea_glom_15 <- tip_glom(sea_optim_tree, h=0.15)
sea_g15_parsimony <- parsimony(sea_glom_15, c_two_phydat)
sub_title_15_sea <- (paste("Parsimony: ", sea_g15_parsimony))
plot(sea_glom_15, cex=0.4, main='Seawater - Community 2
Tip Agglomerated
(Partial)', sub=sub_title_15_sea)

# TIP AGGLOMERATION TO h = 0.30 & PARSIMONY
sea_glom_30 <- tip_glom(sea_optim_tree, h=0.3)
sea_g30_parsimony <- parsimony(sea_glom_30, c_two_phydat)
sub_title_30_sea <- paste("Parsimony: ", sea_g30_parsimony)
plot(sea_glom_30, cex=0.7, main='Seawater - Community 2
Tip Agglomerated
(Complete)', sub=sub_title_30_sea)


# PLOTTING TREE WITH TAXON LABELS INSTEAD OF ASVs
sea_tip_names <- sea_glom_30$tip.label
sea_glom_30$tip.label <- tax_selected[sea_tip_names,]$Final

plot(sea_glom_30, cex=0.7, main='Seawater - Community 2
Tip Agglomerated with Taxonomies Included
(Complete)', sub=sub_title_30_sea)

#########################################################################


#########################################################################
# SHARED (MERGED) COMMUNITY TREE - SAME PROCESS

# BUILDING TREE
dmat_shared <- dist.ml(c_merged_phydat)
shared_nj <- nj(dmat_shared)

# RE-ROOTING
shared_last_tip <- shared_nj$tip.label[length(shared_nj$tip.label)]
shared_tip_num <- length(shared_nj$tip.label)
reroot_shared <- root(shared_nj, shared_last_tip, resolve.root=TRUE)


# BS & PARSIMONY
pml_shared <- pml(reroot_shared, c_merged_phydat)

# OPTIMIZATION
fit_shared_tbr <- optim.pml(pml_shared, model='JC', rearrangement='tbr')

# BOOTSTRAPPING
bs_shared_tbr <- bootstrap.pml(fit_shared_tbr, bs=100, optNni=TRUE, control=pml.control(trace=0),
                            multicore=TRUE, mc.cores=cores)

# PLOTTING BS TREE
bs_shared_title <- (paste("Bootstrap Support Values - 
Shared Community - Root: ", shared_last_tip))
plotBS(root(fit_shared_tbr$tree, shared_last_tip), bs_shared_tbr, p=70, type='p', 
       main=bs_shared_title, cex=0.4, show.tip.label=FALSE)


# RE-ROOTING
shared_optim_tree <- root(fit_shared_tbr$tree, shared_last_tip, resolve.root=TRUE)


################### SHARED (MERGED) COMMUNITY TREES

# ORIGINAL TREE WITH LABELS
shared_orig_parsimony <- parsimony(shared_optim_tree, c_merged_phydat)
sub_title_orig_shared <- (paste("Parsimony: ", shared_orig_parsimony))
plot(shared_optim_tree, cex=0.4, main='Shared - Both Communities
Full Tree', sub=sub_title_orig_shared)

# ORIGINAL TREE WIHTOUT LABELS
plot(shared_optim_tree, cex=0.4, main='Shared - Both Communities
Full Tree', sub=sub_title_orig_shared, show.tip.label=FALSE)

# # TIP AGGLOMERATION TO h = 0.15 & PARSIMONY
shared_glom_15 <- tip_glom(shared_optim_tree, h=0.15)
shared_g15_parsimony <- parsimony(shared_glom_15, c_merged_phydat)
sub_title_15_shared <- (paste("Parsimony: ", shared_g15_parsimony))
plot(shared_glom_15, cex=0.4, main='Shared - Both Communities 
tip Agglomerated
(Partial)', sub=sub_title_15_shared)

# # TIP AGGLOMERATION TO h = 0.30 & PARSIMONY
shared_glom_30 <- tip_glom(shared_optim_tree, h=0.3)
shared_g30_parsimony <- parsimony(shared_glom_30, c_merged_phydat)
sub_title_30_shared <- paste("Parsimony: ", shared_g30_parsimony)
plot(shared_glom_30, cex=0.7, main='Shared - Both Communities
Tip Agglomerated
(Complete)', sub=sub_title_30_shared)


# PLOTTING TREE WITH TAXON LABELS INSTEAD OF ASVs
shared_tip_names <- shared_glom_30$tip.label
shared_glom_30$tip.label <- tax_selected[shared_tip_names,]$Final

plot(sea_glom_30, cex=0.7, main='Seawater - Community 2
Tip Agglomerated with Taxonomies Included
(Complete)', sub=sub_title_30_sea)

#########################################################################



#########################################################################

# <<<<<<<<<<<<<<<<<<< NEXT: OUTPUT ONE PDF WITH ALL THE PLOTS ON THEM LATER >>>>>>>>>>>


'''pdf(file="Pyrosomes_agglomeration2.pdf", width=10, height=8)
plot(pyr_glom_30, cex=1, main=```Pyrosomes - Tip Agglomerated```, 
     sub=sub_title_30_pyr, cex.sub=2, cex.main=2)
dev.off()'''


# UNIFRAC PORTION IS IN phangorn_6_outgroup_root_attempt.R