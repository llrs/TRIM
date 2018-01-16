#!/usr/bin/make

# Overwrites parent make?

ifeq ($(MAKELEVEL),1)
	R_OPTS=--no-site-file --no-environ --no-restore --no-save
endif

TRIM.RDS sgcca.RData bootstrap.RData: intestinal_16S_RNAseq_metadb.R helper_functions.R $(pre_files)
	@echo "Relating the RNAseq with the microbiota"
	R CMD BATCH $(R_OPTS) $(<F)
	
# IBD only analysis
IBD.RData bootstrap_IBD.RData TRIM_IBD.RDS: IBD_patients.R helper_functions.R $(pre_files)
	@echo "Relating the RNAseq with the microbiota in IBD"
	R CMD BATCH $(R_OPTS) $(<F)

# Test which biological relation exists between the selected genes and microbiota
RNAseq_enrichment.csv gsea_RNAseq_pathways.csv Otus_genus_enrichment.csv gsea_otus_genus.csv: biological_relations.R helper_functions.R TRIM.RDS sgcca.RData bootstrap.RData
	@echo "Finding the meaning of these relationships"
	R CMD BATCH $(R_OPTS) biological_relations.R;\

# Test which biological relation exists between the selected genes and microbiota
RNAseq_enrichment_IBD.csv gsea_RNAseq_pathways_IBD.csv Otus_genus_enrichment_IBD.csv gsea_otus_genus_IBD.csv: biological_relations_IBD.R helper_functions.R IBD.RData bootstrap_IBD.RData TRIM_IBD.RDS
	R CMD BATCH $(R_OPTS) biological_relations_IBD.R

.PHONY: all biological_relations

biological_relations: RNAseq_enrichment_IBD.csv gsea_RNAseq_pathways.csv
	

all: RNAseq_enrichment_IBD.csv RNAseq_enrichment.csv TRIM_IBD.RDS TRIM.RDS
