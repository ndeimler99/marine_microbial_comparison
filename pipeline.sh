#!/bin/bash

while getopts d:m:o: flag
do
    case "${flag}" in
        d) dataset_dir=${OPTARG};;
        m) metadata=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done


#######################################################################################################################################################
#HYPERPARAMETERS BELOW

#filtering (must be integer values)
minimum_sum_ASV=100
min_sample_save=1500
minimum_sample_count=1

#normalization technique can either be "DESeq2", "rarefaction", or "RPM".  rare_faction_lib_size only needs to be adjusted if choosing rarefaction normalization technique
normalization_technique="DESeq2"
rarefaction_lib_size=10000

#bar plots
number_of_clades=10

#anosim (must be less than 1, greater than 0)
min_percentage=0.01
high_cumulative=0.8

#differential abundance
#Shrinkage method can either be "normal", "ashr", "apeglm"
#Transform can either be "rld", "vsd", or "vsd.fast"
shrinkage_method="normal"
transform="rld"
#######################################################################################################################################################

current_dir=$(pwd)


#check for absolute or relative paths and reassign if necesarry
if [[ "$dataset_dir" != /* ]]; then
    dataset_dir=${dataset_dir:1}
    dataset_dir=$current_dir$dataset_dir
fi

if [[ "$output_dir" != /* ]]; then
    output_dir=${output_dir:1}
    output_dir=$current_dir$output_dir
fi

output_dir=$output_dir"/"

#make output directory
mkdir $output_dir

#make output directory for quality plots
mkdir $output_dir/quality_plots

#pass dataset and output to dada script
conda activate R_Environment

###THE BELOW SCRIPT CALLS WILL HAVE TO CHANGE DEPENDING ON HOW FINAL PACKAGE IS SENT/LAYED OUT
./R_scripts/dada.R $dataset_dir $output_dir


#THE BELOW SCRIPT IS REQUIRED AND SHOULD NOT BE REMOVED
#TAKES RAW ASV COUNTS, FILTERS THEM, AND SPLITS ASV COUNT DATA BY COMMUNITY THROUGH THE USE OF METADATA FILE
./R_scripts/data_manipulation_and_filtering.R $minimum_sum_ASV $minimum_sample_count $min_sample_save $metadata $output_dir

#Rarefaction Curve, this was not working in R studio, maybe it will on talapas
#YOU CANNOT DRAW CONCLUSIONS ABOUT SPECIES RICHNESS FROM RAREFACTION CURVES.  THIS SIMPLY TELLS YOU IF YOUR SAMPLING IS SUFFICIENT
./R_scripts/rarefaction.R $output_dir


#Normalization
mkdir $output_dir/normalization_results
#rarefaction_lib_size must only be modified if normalization_technique="rarefaction"
./R_scripts/normalization.R $output_dir $normalization_technique $rarefaction_lib_size

#Phyloseq Object Creation
./R_scripts/phyloseq_object_creation.R $output_dir

mkdir $output_dir/stacked_bar_plots/
#stacked bar plots
./R_scripts/bar_plots.R $output_dir $number_of_clades

#heatmaps
mkdir $output_dir/heat_maps/
./R_scripts/heat_maps.R $output_dir
#NOT COMMENTED

#alpha diversity
mkdir $output_dir/alpha_diversity/
./R_scripts/alpha_diversity.R $output_dir

#bray curtis matrices and nmds plots
mkdir $output_dir/bray_nmds/
./R_scripts/bray_curtis_nmds.R $output_dir

#anosim, simper, permanova/adonis
mkdir $output_dir/anosim_simper/
mkdir $output_dir/anosim_simper/community_one/
mkdir $output_dir/anosim_simper/community_two/
mkdir $output_dir/adonis/
./R_scripts/anosim_simper_adonis.R $output_dir $low_percentage $cum_sum

#differential_abundance
mkdir $output_dir/differential_expression/
./R_scripts/differential_expression.R $output_dir $shrinkage_method $transform

#taxonomic summary
mkdir $output_dir/taxonomic_summary
./R_scripts/identification_summary.R $output_dir

#major summary statistics
./R_scripts/summary.R $output_dir

#phylogenetic analysis
mkdir $output_dir/phylogenetics
mkdir $output_dir/phylogenetics/community_one
mkdir $output_dir/phylogenetics/community_two
mkdir $output_dir/phylogenetics/community_merged
./R_scripts/phy_tree_1.7.21.R $output_dir
