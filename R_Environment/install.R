#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
cran <- args[1]

install.packages("versions", repos=cran, quiet=TRUE)
library(versions)

#install package "versions"

#install packages
install.versions("phangorn", "2.5.5", quiet=TRUE)
install.versions("MASS", "7.3-53", quiet=TRUE)
install.versions("ape", "5.4-1", quiet=TRUE)
install.versions("ashr", "2.2-47", quiet=TRUE)
install.versions("tidyverse", "1.3.0", quiet=TRUE)
install.versions("ggplot2", "3.3.2", quiet=TRUE)
install.versions("fossil", "0.4.0", quiet=TRUE)
install.versions("sys", "3.4", quiet=TRUE)
install.versions("vegan", "2.5-7", quiet=TRUE)

#install Biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos=cran, ask=FALSE)
BiocManager::install(version = "3.11", ask=FALSE, quiet=TRUE)

#install biocmanager packages
BiocManager::install("dada2", update=FALSE, quiet=TRUE)
BiocManager::install("DESeq2", update=FALSE, quiet=TRUE)
BiocManager::install("biostrings", update=FALSE, quiet=TRUE)
BiocManager::install("phyloseq", update=FALSE, quiet=TRUE)
BiocManager::install("edgeR", update=FALSE, quiet=TRUE)
BiocManager::install("msa", update=FALSE, quiet=TRUE)
BiocManager::install("microbiome", update=FALSE, quiet=TRUE)
BiocManager::install("apeglm", update=FALSE, quiet=TRUE)


#load silva databases
download.file("https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz?download=1", destfile = "./silva_nr99_v138_train_set.fa.gz")
download.file("https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz?download=1", destfile = "./silva_species_assignment_v138.fa.gz")