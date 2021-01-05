#!/usr/bin/env Rscript

library(phyloseq)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
output_dir <- args[1]

load(paste(output_dir, "r_save", sep=""))

taxonomic_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

for(level in taxonomic_levels){
    plot <- plot_heatmap(community_one_norm_phy_obj, method="NMDS", distance="bray", taxa.label=level, taxa.order=level) +
        ggtitle(paste(community_one_name, "normalized", level, "Bray Heatmap", sep=""))

    ggsave(filename=paste(community_one_name, "_normalized_bray_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/heat_maps/", sep=""))
}

for(level in taxonomic_levels){
    plot <- plot_heatmap(community_one_norm_phy_obj, sample.order=colnames(community_one), taxa.label=level, taxa.order=level) +
        ggtitle(paste(community_one_name, "normalized ", level, "heatmap", sep=""))

    ggsave(filename=paste(community_one_name, "_normalized_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/heat_maps/", sep=""))
}

for(level in taxonomic_levels){
    plot <- plot_heatmap(community_two_norm_phy_obj, method="NMDS", distance="bray", taxa.label=level, taxa.order=level) +
        ggtitle(paste(community_two_name, "normalized", level, "Bray Heatmap", sep=""))

    ggsave(filename=paste(community_two_name, "_normalized_bray", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/heat_maps/", sep=""))
}

for(level in taxonomic_levels){
    plot <- plot_heatmap(community_two_norm_phy_obj, sample.order=colnames(community_two), taxa.label=level, taxa.order=level) +
        ggtitle(paste(community_two_name, "normalized ", level, "heatmap", sep=""))

    ggsave(filename=paste(community_two_name, "_normalized_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/heat_maps/", sep=""))
}

for(level in taxonomic_levels){
    plot <- plot_heatmap(community_merge_norm_phy_obj, method="NMDS", distance="bray", taxa.label=level, taxa.order=level) +
        ggtitle(paste("All Samples normalized", level, "Bray Heatmap", sep=""))

    ggsave(filename=paste("all_samples_normalized_bray_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/heat_maps/", sep=""))
}

for(level in taxonomic_levels){
    plot <- plot_heatmap(community_merge_norm_phy_obj, sample.order=colnames(community_merge), taxa.label=level, taxa.order=level) +
        ggtitle(paste("All Samples normalized ", level, "heatmap", sep=""))

    ggsave(filename=paste("all_samples_normalized_", level, ".png", sep=""), plot=plot, device="png", path=paste(output_dir, "/heat_maps/", sep=""))
}
save.image(paste(output_dir, "r_save", sep=""))
print("Heat maps Done")