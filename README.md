# Purpose

This repository analyzes the relationship between the microbiome and the transcriptome on IBD.

# Content

Most of the analysis are done on the intestinal_16S_RNAseq_metadb folder.
Some preliminar analysis are done on the main folder, including PCAs.
You can use the makefile for some of the steps.

# Dependencies

You need to have installed several packages from Bioconductor and CRAN mainly. 
If you install from github the package [llrs/integration-helper](https://github.com/llrs/integration-helper) written to help on this type of analysis you it should install most of the dependencies. 
With the notable exception of the package on github  [thomasp85/patchwork](https://github.com/thomasp85/patchwork).

# Data

The data is missing from this repository, and it is supposed to be located in three different folders (one for each kind of data).

 - [RNAseq](intestinal_RNAseq/)
 - [stool 16S](stool_16S/)
 - [intestinal 16S](intestinal_16S/)