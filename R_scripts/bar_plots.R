#!/usr/bin/env Rscript

##############################################################################
#LOAD LIBRARIES AND ASSIGN USER INPUT TO SCRIPT VARIABLES

library(phyloseq)
library(ggplot2)
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
number_of_clades <- as.integer(args[2])

load(paste(output_dir, "r_save", sep=""))
##############################################################################


##############################################################################
#Function that takes a phyloseq object and a taxonomy level that will return a stacked bar plot containing number_of_clades clades
#Will then save two files a stacked bar plot of counts and percent abundance with inputted title

plotter <- function(phy_object, taxonomy_level, title){
    
    temp <- tax_glom(phy_object, taxonomy_level)
    counts <- otu_table(temp)
    rownames(counts) <- tax_table(temp)[,taxonomy_level]
    counts <- data.frame(counts)
    counts$sum <- rowSums(counts)
    counts <- counts[order(counts$sum, decreasing=TRUE),]
    if(length(rownames(counts)) > number_of_clades){
        
        new_counts <- counts[0:number_of_clades,]
        new_counts <- select(new_counts, -sum)
        other_start = number_of_clades + 1
        row_names <- rownames(counts)[0:number_of_clades]
        row_names <- append(row_names, "Other")
        counts <- counts[other_start:length(rownames(counts)),]
        other <- colSums(counts)
        new_counts <- rbind(new_counts, other)
    }
    else{
        new_counts <- select(counts, -sum)
        row_names <- rownames(counts)
    }
    rownames(new_counts) <- row_names

    df <- NULL

    new_counts <- as.data.frame(t(new_counts))
    for(name in colnames(new_counts)){
        clade <- rep(name, length(rownames(new_counts)))
        sample <- rownames(new_counts)
        count <- new_counts[,name]
        temp_df <- data.frame(sample, clade, count)
        df <- rbind(df, temp_df)
    }

    df <- as.data.frame(df)
    df$clade <- factor(df$clade)
    df$sample <- factor(df$sample, levels=unique(df$sample))


    plot <- ggplot(data=df) +
        geom_bar(mapping=aes(x=sample, y=count, fill=clade), stat="identity", position="fill") +
        ggtitle(paste(title, "Percent Abundance")) +
        theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color="black", angle=30, hjust=0.5)) +
        xlab("Sample") +
        ylab("Percent Abundance") +
        labs(fill=taxonomy_level)

    ggsave(filename=paste(title, "_percent_abundance.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

    plot <- ggplot(data=df) +
        geom_bar(mapping=aes(x=sample, y=count, fill=clade), stat="identity", position="stack") +
        ggtitle(paste(title, "Counts")) +
        theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(color="black", angle=30, hjust=0.5)) +
        xlab("Sample") +
        ylab("Count") +
        labs(fill=taxonomy_level)

    ggsave(filename=paste(title, "_counts.png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
}
##############################################################################


##############################################################################
#Call above function for each community individually as well as both communities normalized and not

taxonomic_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
for(level in taxonomic_levels){
    plotter(community_merge_norm_phy_obj, level, paste(level, " - Both Communities Normalized"))
    plotter(community_merge_phy_obj, level, paste(level, " - Both Community Non-normalized"))
    plotter(community_one_norm_phy_obj, level, paste(level, " - ", community_one_name, " Normalized"))
    plotter(community_one_phy_obj, level, paste(level, " - ", community_one_name, " Non-normalized"))
    plotter(community_two_norm_phy_obj, level, paste(level, " - ", community_two_name, " Normalized"))
    plotter(community_two_phy_obj, level, paste(level, " - ", community_two_name, " Non-normalized"))
}

##############################################################################
#Save R Work Space
save.image(paste(output_dir, "r_save", sep=""))

#Print Tracker
print("Stacked Bar Plots Done")
##############################################################################
