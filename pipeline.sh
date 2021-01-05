#!/bin/bash

#This is the script that will be called to start the program

#it will take the options of where the data sets are stored, a file containing list of file names/corresponding to sample names and any other important information, as well as an output directory

while getopts d:m:o: flag
do
    case "${flag}" in
        d) dataset_dir=${OPTARG};;
        m) metadata=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

#echo "dataset=$dataset_dir";
#echo "metdata: $metadata";
#echo "output: $output_dir";

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


minimum_sum_ASV=100
min_sample_save=1500
minimum_sample_count=1

#THE BELOW SCRIPT IS REQUIRED AND SHOULD NOT BE REMOVED
#TAKES RAW ASV COUNTS, FILTERS THEM, AND SPLITS ASV COUNT DATA BY COMMUNITY THROUGH THE USE OF METADATA FILE
./R_scripts/data_manipulation_and_filtering.R $minimum_sum_ASV $minimum_sample_count $min_sample_save $metadata $output_dir

#Rarefaction Curve, this was not working in R studio, maybe it will on talapas
#YOU CANNOT DRAW CONCLUSIONS ABOUT SPECIES RICHNESS FROM RAREFACTION CURVES.  THIS SIMPLY TELLS YOU IF YOUR SAMPLING IS SUFFICIENT
./R_scripts/rarefaction.R $output_dir


#Normalization
#normalization technique can either be DESeq2, rarefaction, or RPM
mkdir $output_dir/normalization_results
normalization_technique="DESeq2"
#rarefaction_lib_size must only be modified if normalization_technique="rarefaction"
rarefaction_lib_size=10000
./R_scripts/normalization.R $output_dir $normalization_technique $rarefaction_lib_size

#Phyloseq Object Creation
./R_scripts/phyloseq_object_creation.R $output_dir

mkdir $output_dir/stacked_bar_plots/
#stacked bar plots
./R_scripts/bar_plots.R $output_dir
#NOT COMMENTED

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
min_percentage=0.01
high_cumulative=0.8
./R_scripts/anosim_simper_adonis.R $output_dir $low_percentage $cum_sum

#Shrinkage method can either be "normal", "ashr", "apeglm"
#Transform can either be "rld"
shrinkage_method="ashr"
transform="rld"
mkdir $output_dir/differential_expression/
./R_scripts/differential_expression.R $output_dir $shrinkage_method $transform

mkdir $output_dir/taxonomic_summary
./R_scripts/identification_summary.R $output_dir

./R_scripts/summary.R $output_dir

mkdir $output_dir/phylogenetics
#./R_scripts/phangorn_final.R $output_dir/phylogenetics
