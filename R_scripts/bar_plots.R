#!/usr/bin/env Rscript

library(phyloseq)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
load(paste(output_dir, "r_save", sep=""))

taxonomic_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

plotter <- function(phy_object, taxonomy_level){

    temp_tax_table <- data.frame(tax_table(phy_object), header=TRUE)
    temp_otu_table <- data.frame(otu_table(phy_object), header=TRUE)
    truths <- temp_tax_table[is.na(temp_tax_table[taxonomy_level])==FALSE,]
    otus <- temp_otu_table[rownames(truths),]
    
    otus <- dplyr::select(otus, -header)
    truths <- dplyr::select(truths, -header)

    otus <- otus[order(rowSums(otus), decreasing=TRUE),]
    truths <- truths[order(rowSums(otus), decreasing=TRUE),]

    otus <- otu_table(otus, taxa_are_rows=TRUE)
    truths <- tax_table(as.matrix(truths))
    
    object <- phyloseq(otus, truths)
    phy_object <- tax_glom(phy_object, taxonomy_level)
    object <- tax_glom(object, taxonomy_level)
    plot_NA <- plot_bar(phy_object, fill=taxonomy_level)+
        geom_bar(stat="identity")
    plot_noNA <- plot_bar(object, fill=taxonomy_level) + 
        geom_bar(stat="identity")

    plots <- list(plot_NA, plot_noNA)
    return(plots)
}

for(level in taxonomic_levels){

    plots <- plotter(community_one_norm_phy_obj, level)
    plot <- plots[[1]] +
        ggtitle(paste(community_one_name, "normalized", level, sep=" "))
    ggsave(filename=paste(community_one_name, "_normalized_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
     
    plot <- plots[[2]] +
        ggtitle(paste(community_one_name, "normalized", level, sep=" "))
    ggsave(filename=paste(community_one_name, "_normalized_NA_removed", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

}

for(level in taxonomic_levels){
    plots <- plotter(community_two_norm_phy_obj, level)

    plot <- plots[[1]] +
        ggtitle(paste(community_two_name, "normalized", level, sep=" "))
    ggsave(filename=paste(community_two_name, "_normalized_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
     
    plot <- plots[[2]] +
        ggtitle(paste(community_two_name, "normalized", level, sep=" "))
    ggsave(filename=paste(community_two_name, "_normalized_NA_removed", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

}

for(level in taxonomic_levels){
    plots <- plotter(community_merge_norm_phy_obj, level)

    plot <- plots[[1]] +
        ggtitle(paste("all samples normalized", level, sep=" "))
    ggsave(filename=paste("all_samples_normalized_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
     
    plot <- plots[[2]] +
        ggtitle(paste("all samples normalized", level, sep=" "))
    ggsave(filename=paste("all_samples_normalized_NA_removed", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

}

#below will be for perc abundance
for(level in taxonomic_levels){
    plots <- plotter(community_one_phy_obj_perc_abundance, level)

    plot <- plots[[1]] +
        ggtitle(paste(community_one_name, "normalized", level, sep=" "))
    ggsave(filename=paste(community_one_name, "_percent_abundance_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
     
    plot <- plots[[2]] +
        ggtitle(paste(community_one_name, "Percent Abundance", level, sep=" "))
    ggsave(filename=paste(community_one_name, "_percent_abundance_NA_removed", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

}

for(level in taxonomic_levels){
    plots <- plotter(community_two_phy_obj_perc_abundance, level)

    plot <- plots[[1]] +
        ggtitle(paste(community_two_name, "percent abundance", level, sep=" "))
    ggsave(filename=paste(community_two_name, "_percent_abundance_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
     
    plot <- plots[[2]] +
        ggtitle(paste(community_two_name, "percent abundance NA Removed", level, sep=" "))
    ggsave(filename=paste(community_two_name, "_percent_abundance_NA_removed_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

}

for(level in taxonomic_levels){
    plots <- plotter(community_merge_phy_obj_perc_abundance, level)

    plot <- plots[[1]] +
        ggtitle(paste("all samples percent abundace", level, sep=" "))
    ggsave(filename=paste("all_samples_percent_abundance_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))
     
    plot <- plots[[2]] +
        ggtitle(paste("all samples percent_abundance NA removed", level, sep=" "))
    ggsave(filename=paste("all_samples_percent_abundance_NA_removed_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/stacked_bar_plots/", sep=""))

}
save.image(paste(output_dir, "r_save", sep=""))
print("Stacked Bar Plots Done")