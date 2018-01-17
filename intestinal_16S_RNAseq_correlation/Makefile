#!/usr/bin/make

ifeq ($(MAKELEVEL),1)
	R_OPTS=--no-site-file --no-environ --no-restore --no-save
endif

correlations.RDS: intestinal_16S_RNAseq_correlation.R ../helper_functions.R ../intestinal_16S/OTUs-Table-new-biopsies.csv intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt intestinal_RNAseq/table.counts.results intestinal_RNAseq/111217_metadata.csv
	@echo "Relating the RNAseq with the microbiota"
	R CMD BATCH $(R_OPTS) $(<F)
	
.PHONY: all 

all: correlations.RDS
