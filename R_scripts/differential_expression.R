#!/usr/bin/env Rscript

##############################################################################
#LOAD LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES

#LOAD LIBRARIES
library("DESeq2")
library("phyloseq")
library("ggplot2")
library("ashr")

#ASSIGN INPUT TO SCRIPT VARIABLES
args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
shrinkage_method <- args[2]
transform <- args[3]

#LOAD R WORKSPACE
load(paste(output_dir, "r_save", sep=""))

##############################################################################



##############################################################################
#PREPARING FOR DESEQ2 DIFFERENTIAL ABUNDANCE ANALYSIS

#CREATE NEW TABLE CONTAINING METADATA
metadata_rearranged <- metadata
#SAVING CURRENT ROWNAMES IN NEW COLUMN
metadata_rearranged$file_name <- rownames(metadata)
#SETTING ROWNAMES TO SAMPLE NAMES INSTEAD OF FILE NAMES
rownames(metadata_rearranged) <- metadata$sample_name

#ASSIGNING MODIFIED METADATA TABLE TO SAMPLE_DATA OF PHYLOSEQ OBJECT
sample_data(community_merge_phy_obj) <- metadata_rearranged

##############################################################################


##############################################################################
#DESEQ2 DIFFERENTIAL ABUNDANCE ANALYSIS

#CONVERT PHYLOSEQ OBJECT TO DESEQ OBJECT, SEPERATING BY SAMPLE GROUP
deseq_merge <- phyloseq_to_deseq2(community_merge_phy_obj, ~sample_group)

#ACTUAL DIFFERENTIAL ANALYSIS
deseq_merge <- DESeq(deseq_merge)
deseq_merge <- deseq_merge[rowSums(counts(deseq_merge)) > 1,]

#SHRINK DATA (REDUCES LOG FOLD CHANGE AND SIGNIFICANCE OF ASVS WITH A LARGE DEGREE OF VARIANCE)
if(shrinkage_method=="normal"){
    shrunk <- lfcShrink(deseq_merge, coef=2, type="normal")
}

if(shrinkage_method=="ashr"){
    shrunk <- lfcShrink(deseq_merge, coef=2, type="ashr")
}

if(shrinkage_method=="apeglm"){
    shrunk <- lfcShrink(deseq_merge, coef=2, type="apeglm")
}

#GET RESULTS OF SHRUNK DIFFERENTIAL ABUNDANCE
differential_expression_results <- shrunk
#SORT RESULTS WITH MOST SIGNIFICANT FIRST
differential_expression_results <- differential_expression_results[order(differential_expression_results$padj),]

#WRITE ALL RESULTS TO TABLE
write.table(data.frame(differential_expression_results), file=paste(output_dir, "/differential_expression/results_table.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#REMOVE ALL ASVS WHERE PADJ == null
expression_results <- differential_expression_results[is.na(differential_expression_results$padj) == FALSE,]
#SAVED ASVS IN WHICH PADJ IS LESS THAN 0.01 SIGNIFICANCE
differentially_expressed <- expression_results[expression_results$padj < 0.01,]
#WRITE SIGNIFICANT ASVS TO NEW TABLE
write.table(data.frame(differentially_expressed), file=paste(output_dir, "/differential_expression/differentially_expressed.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#OPEN FILE FOR SUMMARY OF ASV DIFFERENTIAL ABUNDANCE
sink(file=paste(output_dir, "/differential_expression/all_asv_summary.txt", sep=""))
summary(differential_expression_results)
sink()

#OPEN FILE FOR LOG2FOLD CHANGE OF ASV DIFFERENTIAL ABUNDANCE
png(paste(output_dir, "/differential_expression/all_asvs_MA.png", sep=""))
plotMA(shrunk, alpha=0.01, main="Differential Expression of all ASVs")
while (!is.null(dev.list()))  dev.off()


#PCA PLOTS USING VARIOUS METHODS OF TRANSFORMATION
if(transform=="rld"){   
    transformed <- rlog(deseq_merge, blind=FALSE)
}
if(transform=="vsd"){
    transformed <- varianceStabilizingTransformation(deseq_merge, blind=FALSE)
}
if(transform=="vsd.fast"){
    transformed <- vst(deseq_merge, blind=FALSE)
}

#PLOT PCA
plot <- plotPCA(transformed, intgroup=c("sample_group"))
#SAVE PCA
ggsave(filename="all_asvs_pca.png", plot=plot, device="png", path=paste(output_dir, "/differential_expression/", sep=""))

##############################################################################


##############################################################################
#FURTHER ANALYSIS OF DIFFERENTIAL ABUNDANCE RESULTS

#IF LOG FOLD CHANGE IS LESS THAN 0, IT IS DIFFERENTIALLY ABUNDANT IN COMMUNITY ONE
community_one_differential <- expression_results[expression_results$log2FoldChange < 0,]
#DETERMINE HOW MANY ASVS ARE DIFFERENTIALLY ABUNDANT IN COMMUNITY ONE
community_one_differential_all <- length(rownames(community_one_differential))
#SIGNIFICANTLY DIFFERENTIAL ABUNDANT IN COMMUNITY ONE
community_one_differential <- community_one_differential[community_one_differential$padj < 0.01,]
#DETERMINE HOW MANY ASVS ARE SIGNIFICANT
community_one_differential_significant <- length(rownames(community_one_differential))
#REMOVE ANY ASVS IN WHICH PADJ IS NULL
community_one_differential <- community_one_differential[is.na(community_one_differential$padj) == FALSE,]

#IF LOG FOLD CHANGE IS GREATER THAN 0, IT IS DIFFERENTIALLY ABUNDANT IN COMMUNITY TWO
community_two_differential <- expression_results[expression_results$log2FoldChange > 0,]
#DETERMINE HOW MANY ASVS ARE DIFFERENTIALLY ABUNDANT IN COMMUNITY TWO
community_two_differential_all <- length(rownames(community_two_differential))
#sIGNIFICANTLY DIFFERENTIAL ABUNDANT IN COMMUNITY TWO
community_two_differential <- community_two_differential[community_two_differential$padj < 0.01,]
#DETERMINE HOW MANY ASVS ARE SIGNIFICANT
community_two_differential_significant <- length(rownames(community_two_differential))
#REMOVE ANY ASVS IN WHICH PADJ IS NULL
community_two_differential <- community_two_differential[is.na(community_two_differential$padj) == FALSE,]

#DETERMINE IF THOSE THAT ARE DIFFERENTIALLY ABUNDANT ALSO APPEAR IN THE OTHER COMMUNITY AS WELL
#CREATE EMPTY LIST
cross <- NULL
community_one_differential_in_two <- 0
#FOR ASV IN COMMUNITY ONE DIFFERENTIAL
for(name in rownames(community_one_differential)){
    #IF THIS ASV IS IN COMMUNITY TWO COUNT DATAFRAME APPEND A YES TO THE LIST, ELSE APPEND A NO TO THE LIST
    if(name %in% rownames(community_two)){
        cross <- append(cross, "Yes")
        community_one_differential_in_two <- community_one_differential_in_two + 1
    }
    else{
        cross <- append(cross, "No")
    }
}

#ASSIGN AS NEW COLUMN IN DIFFERENTIAL ABUNDANCE RESULTS
community_one_differential$`community_two_name` <- cross


#REPEAT FOR COMMUNITY TWO SIGNIFICANCE WHILE COMPARING TO COMMUNITY ONE COUNT DATAFRAME
community_two_differential_in_one <- 0
cross <- NULL
for(name in rownames(community_two_differential)){
    if(name %in% rownames(community_one)){
        cross <- append(cross, "Yes")
        community_two_differential_in_one <- community_two_differential_in_one + 1
    }
    else{
        cross <- append(cross, "No")
    }
}

community_two_differential$`community_one_name` <- cross

#DETERMINE LIST OF ASVS THAT ARE EXTREMELY DIFFERENTIALLY ABUNDANT IN EACH COMMUNITY
extremely_different_community_one <- community_one_differential[community_one_differential$padj < 0.00001,]
extremely_different_community_two <- community_two_differential[community_two_differential$padj < 0.00001,]

#WRITE TO TABLES
write.table(community_one_differential, file=paste(output_dir, "/differential_expression/", community_one_name, "_differential_0.01.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(community_two_differential, file=paste(output_dir, "/differential_expression/", community_two_name, "_differential_0.01.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(extremely_different_community_one, file=paste(output_dir, "/differential_expression/", community_one_name, "_differential_0.00001.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")
write.table(extremely_different_community_two, file=paste(output_dir, "/differential_expression/", community_two_name, "_differential_0.00001.txt", sep=""), row.names=TRUE, col.names=TRUE, append=FALSE, sep="\t")

#PLOT DIFFERENTIAL ABUNDANCE RESULTS (Y=LOG2FOLDCHANGE)
plot <- ggplot() +
    geom_point(mapping=aes(x=rownames(differential_expression_results), y=differential_expression_results$log2FoldChange)) +
    geom_point(mapping=aes(x=rownames(extremely_different_community_one), y=extremely_different_community_one$log2FoldChange), color="red") +
    geom_point(mapping=aes(x=rownames(extremely_different_community_two), y=extremely_different_community_two$log2FoldChange), col="blue")

#SAVE PLOT
ggsave(filename="differential_expression_plot.png", plot=plot, device="png", path=paste(output_dir, "/differential_expression/", sep=""))
##############################################################################


##############################################################################
#DIFFERENTIAL ABUNDANCY BY TAXONOMIC LEVEL

taxonomic_differentiation <- function(level){
    #CONGLOMERATE ALL ASVS THAT BELONG TO EACH UNIQUE CLASSIFICATION AT DESIRED TAXONOMIC LEVEL
    temp <- tax_glom(community_merge_phy_obj, level)
    
    #ASSIGN SAMPLE DATA TO NEW CONGLOMERATED PHYLOSEQ OBJECT
    sample_data(temp) <- metadata_rearranged
    #CONVERT PHYLOSEQ OBJECT TO DESEQ OBJECT
    deseq_merge <- phyloseq_to_deseq2(temp, ~sample_group)
    #DESEQ DIFFERENTIAL ABUNDANCE
    deseq_merge <- DESeq(deseq_merge, fitType="mean")

    #SHRINK USING DESIRED METHOD
    if(shrinkage_method=="normal"){
        shrunk <- lfcShrink(deseq_merge, coef=2, type="normal")
    }

    if(shrinkage_method=="ashr"){
        shrunk <- lfcShrink(deseq_merge, coef=2, type="ashr")
    }

    if(shrinkage_method=="apeglm"){
        shrunk <- lfcShrink(deseq_merge, coef=2, type="apeglm")
    }

    #EXTRACT RESULTS FROM DESEQ OBJECT
    deseq <- shrunk
    #CHANGE TO DATAFRAME
    expression <- data.frame(deseq)
    #ASSIGN APPROPRIATE ROWNAMES TO DATAFRAME
    rownames(expression) <- data.frame(tax_table(temp))[,level]
    #WRITE DIFFERENTIAL ABUNDANCE RESULTS FOR EACH UNIQUE CLASSIFICATION TO TABLE
    write.table(expression, file=paste(output_dir, "/differential_expression/", level, "_differential.txt", sep=""), sep="\t", row.names=TRUE, col.names=TRUE, append=FALSE)
    
    #REMOVE PADJ==NULL
    deseq <- deseq[is.na(deseq$padj)==FALSE,]
    #EXTRACT SIGNIFICANT DIFFERENTIAL ABUNDANCE
    deseq_merge_sig <- deseq[deseq$padj < 0.01,]

    #PLOT DIFFERENTIAL ABUNDANCE
    plot <- ggplot() +  
        geom_point(mapping=aes(x = rownames(deseq), y=deseq$log2FoldChange)) +
        geom_point(mapping=aes(x=rownames(deseq_merge_sig), y= deseq_merge_sig$log2FoldChange), col="Red")
    #SAVE PLOT
    ggsave(filename=paste(level, "_differential_expression.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/differential_expression/", sep=""))

    #OPEN FILE FOR SUMMARY OF DIFFERENTIAL ABUNDANCE
    sink(file=paste(output_dir, "/differential_expression/", level, "_differential_expression_summary.txt", sep=""))
    print(summary(deseq))
    sink()

    #OPEN FILE TO PLOT DIFFERENTIAL ABUNDANCE
    png(paste(output_dir, "/differential_expression/", level, "_MA.png", sep=""))
    plotMA(deseq, alpha=0.1, main="Differential Expression of all ASVs")
    while (!is.null(dev.list()))  dev.off()

    #TRANSFORM DATA
    if(transform == "rld"){   
        transformed <- rlog(deseq_merge, blind=FALSE)
    }
    if(transform == "vsd"){
        transformed <- varianceStabilizingTransformation(deseq_merge, blind=FALSE)
    }
    if(transform == "vsd.fast"){
        transformed <- vst(deseq_merge, blind=FALSE)
    }

    #PLOT AND SAVE DATA
    plot <- plotPCA(transformed, intgroup=c("sample_group"))
    ggsave(filename=paste(level, "_pca.png", sep=""), plot=plot, device=png, path=paste(output_dir, "/differential_expression/", sep=""))
}

taxonomic_differentiation("Kingdom")
taxonomic_differentiation("Phylum")
taxonomic_differentiation("Class")
taxonomic_differentiation("Order")
taxonomic_differentiation("Family")
taxonomic_differentiation("Genus")

##############################################################################


#DETERMINE WHICH ASV HAD LARGEST LOG_FOLD CHANGE AND LOWEST PADJ FROM COMMUNITY ONE
community_one_log_fold_max <- rownames(community_one_differential[which.min(community_one_differential$log2FoldChange),])
community_one_p_min <- rownames(community_one_differential[which.min(community_one_differential$padj),])

#DETERMINE WHICH ASV HAD LARGEST LOG FOLD CHANGE AND LOWEST PADJ FROM COMMUNITY TWO
community_two_log_fold_max <- rownames(community_two_differential[which.max(community_two_differential$log2FoldChange),])
community_two_p_min <- rownames(community_two_differential[which.min(community_two_differential$padj),])

#SAVE R WORKSPACE
save.image(paste(output_dir, "r_save", sep=""))

#PRINT TRACKER
print("differential expression done")