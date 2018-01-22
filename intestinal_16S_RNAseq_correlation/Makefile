#!/usr/bin/make

ifeq ($(MAKELEVEL),1)
	R_OPTS=--no-site-file --no-environ --no-restore --no-save
endif

pre_files=../intestinal_16S/OTUs-Table-new-biopsies.csv \
	../intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt \
	../intestinal_RNAseq/table.counts.results \
	../intestinal_RNAseq/111217_metadata.csv
	
correlations.RDS: intestinal_16S_RNAseq_correlation.R ../helper_functions.R $(pre_files)
	@echo "Relating the RNAseq with the microbiota"
	R CMD BATCH $(R_OPTS) $(<F)
	
.PHONY: all filter

post_files=../intestinal_16S_RNAseq_metadb/bootstrap_IBD.RData \
	../intestinal_16S_RNAseq_metadb/bootstrap.RData \
	../intestinal_16S_RNAseq_metadb/sgcca.RDS \
	../intestinal_16S_RNAseq_metadb/IBD.RDS \
	correlations.RDS

filter: filter.R ../helper_functions.R $(post_files)
	@echo "Filtering the correlations"
	R CMD BATCH $(R_OPTS) $(<F)

all: correlations.RDS filter
