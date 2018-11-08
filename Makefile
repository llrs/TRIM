#!/usr/bin/make

R_OPTS=--no-site-file --no-environ --no-restore --no-save

pre_files=stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt \
					stools_16S/OTUs-Table-refined-stools.tab \
					intestinal_16S/OTUs-Table-new-biopsies.csv \
					intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt \
					intestinal_RNAseq/metadata_25042018.csv \
					intestinal_RNAseq/table.counts.results

out_files=meta_coherent.csv \
					stools_16S/otus_coherent.csv \
					stools_16S/taxonomy.csv \
					intestinal_16S/otus_coherent.csv \
					intestinal_16S/taxonomy.csv \
					equivalent_otus.csv \
					equivalent_genus.csv

.PHONY: all STATegRa eqSpecies eqGenus ileum_integration colon_integration PCA \
stools_prevalence intestinal_prevalence prevalence upset intestinal_16S_RNAseq_correlation \
alpha cleaning variance GSV models_report 16S

all: eqGenus eqSpecies \
stool_intestinal_16S_integration/important_common_microrg.csv PCA Deconvolute \
STATegRa intestinal_16S_RNAseq_integration prevalence models_report

# Clean the input and prepare the output for integration
$(out_files): cleaning.R $(pre_files) 

cleaning: $(out_files)
	@echo "Preparing input data"
	R CMD BATCH $(R_OPTS) cleaning.R

	
# Handles the use of RGCCA in stools 16S and biopsies 16S
si16S=stool_intestinal_16S_integration
$(si16S)/%.csv: $(si16S)/stool_intestinal_integration.R $(out_files) 
	@echo "Integrating stools and biopsies 16S data" 
	@echo "\tUsing RGCCA"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

16S: $(si16S)/%.csv
	
# Handles the use of STATegRa in stools 16S and biopsies 16S
STATegRa: $(si16S)/STATegRa.R $(out_files) 
	@echo "Integrating stools and biopsies data" 
	@echo "\tUsing STATegRa"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
# Integrates via correlations the OTUs of the same species
eqSpecies: stool_intestinal_16S_correlation/stool_intestinal_16S_otus.R equivalent_otus.csv stool_intestinal_16S_metadb/important_common_microrg.csv 
	@echo "Analyse the microorganisms in common between stools and biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
# Integrates via correlations the OTUs of the same genus
eqGenus: stool_intestinal_16S_correlation/stool_intestinal_16S_genus.R equivalent_genus.csv stool_intestinal_16S_metadb/important_common_microrg.csv 
	@echo "Analyse the microorganisms in common between stools and biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Save
FOLDER=intestinal_RNAseq_metadb
$(FOLDER)/otus_table.RDS $(FOLDER)/otus_tax.RDS $(FOLDER)/expr.RDS $(FOLDER)/meta.RDS: $(FOLDER)/preprocessing.R
	@echo "Preprocessing"
	cd $(FOLDER); R CMD BATCH $(R_OPTS)/preprocessing.R

# Code to integrate clinical variables, biopsies 16S and stools 16S
stool_intestinal_metadb/%.csv: stool_intestinal_metadb/stool_intestinal_metadb.R $(out_files) 
	@echo "Integrating stools and biopsies 16S data \
	taking into account the metadata"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
 	
# Handles the integration of only the ileum between biopsies 16S and stools 16S
ileum_integration: ileum_integration/ileum_integration.R $(out_files) 
	@echo "Integrating stools and biopsies data from ileum"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the integration of only the ileum between biopsies 16S and stools 16S
colon_integration: colon_integration/colon_integration.R $(out_files) 
	@echo "Integrating stools and biopsies data from colon"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCAs	
PCA: PCAs.R $(pre_files) 
	@echo "PCAs of the data"
	R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for stools but controlling for clinical variables
stool_metadb: stool_metadb/stool_metadb.R $(pre_files) 
	@echo "PCAs controlling for clinical variables of 16S stools samples"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for biopsies 16S but controlling for clinical variables
intestinal_16S_metadb: intestinal_16S_metadb/intestinal_16S_metadb.R $(pre_files) 
	@echo "PCAs controlling for clinical variables of 16S biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for biopsies RNAseq but controlling for clinical variables
intestinal_RNAseq_metadb: $(FOLDER)/intestinal_RNAseq_metadb.R $(pre_files) 
	@echo "PCAs controlling for clinical variables of RNAseq biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

intestinal_16S_RNAseq_correlation: intestinal_16S_RNAseq_correlation/intestinal_16S_RNAseq_correlation.R intestinal_16S_RNAseq_correlation/filter.R intestinal_RNAseq_metadb
	cd $(<D); CMD BATCH $(R_OPTS) $(<F)
	
	
# Handles the integration between the biopsies
i16S_integration=intestinal_16S_RNAseq_integration

$(i16S_integration)/ctrls_%.RDS $(i16S_integration)/IBD_%.RDS: $(i16S_integration)/preprocessing.R $(pre_files)
	@echo "Preprocessing integration"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

$(i16S_integration)/otus_table.RDS $(i16S_integration)/otus_tax.RDS $(i16S_integration)/expr.RDS $(i16S_integration)/meta.RDS: $(i16S_integration)/preprocessing.R $(pre_files)
	@echo "Preprocessing integration"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
$(i16S_integration)/sgcca.RDS $(i16S_integration)/loo-model0.RDS $(i16S_integration)/bootstrap.RDS: $(i16S_integration)/intestinal_16S_RNAseq_integration.R $(i16S_integration)/otus_table.RDS $(i16S_integration)/otus_tax.RDS $(i16S_integration)/expr.RDS $(i16S_integration)/meta.RDS
	@echo "Integrating"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

$(i16S_integration)/IBD_patients.R $(i16S_integration)/Controls_patients.R: $(i16S_integration)/ctrls_%.RDS $(i16S_integration)/IBD_%.RDS

$(i16S_integration)/%.csv: $(i16S_integration)/biological_relations.R
	@echo "Biological meaning of the integration"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

intestinal_16S_RNAseq_integration: $(i16S_integration)/sgcca.RDS

# Test which part of the variance is explained by what
models.RData variance: intestinal_16S_RNAseq_metadb/variance.R $(pre_files)
	@echo "Testing the variance on clinical data"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

FOLDER2=intestinal_16S_RNAseq_metadb

$(FOLDER2)/otus_table.RDS: $(FOLDER2)/preprocessing.R $(pre_files)
$(FOLDER2)/otus_tax.RDS: $(FOLDER2)/preprocessing.R $(pre_files)
$(FOLDER2)/expr.RDS: $(FOLDER2)/preprocessing.R $(pre_files)
$(FOLDER2)/meta.RDS: $(FOLDER2)/preprocessing.R $(pre_files)
	@echo "Preprocessing"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the integration between the biopsies taking into account the metadata
$(FOLDER2)/sgcca.RDS: $(FOLDER2)/intestinal_16S_RNAseq_metadb.R $(FOLDER2)/IBD_patients.R $(FOLDER2)/Controls_patients.R $(FOLDER2)/otus_table.RDS $(FOLDER2)/otus_tax.RDS $(FOLDER2)/expr.RDS $(FOLDER2)/meta.RDS
	@echo "Relating the RNAseq with the microbiota"
	cd $(<D);\
	R CMD BATCH $(R_OPTS) $(<F);\
	R CMD BATCH $(R_OPTS) Controls_patients.R; \
	R CMD BATCH $(R_OPTS) IBD_patients.R

# Test which biological relation exists between the selected genes and microbiota
$(FOLDER2)/%.csv: $(FOLDER2)/biological_relations.R
	@echo "Finding the meaning of these relationships"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

biological_relations: $(FOLDER2)/%.csv $(i16S_integration)/%.csv
intestinal_16S_RNAseq_metadb: biological_relations







# Handles the calculation of the prevalence in intestinal 
intestinal_prevalence: intestinal_16S_conceptual/prevalence.R $(pre_files)
	@echo "Biopsies prevalence"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the calculation of the prevalence in stools	
stools_prevalence: stools_conceptual/prevalence.R $(pre_files)
	@echo "Stools prevalence"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Do both prevalences
prevalence: intestinal_prevalence stools_prevalence

# Report models Leave-one-out
$(FOLDER2)/models_report.R: $(FOLDER2)/sgcca.RDS $(i16S_integration)/sgcca.RDS $(i16S_integration)/loo-model0.RDS $(FOLDER2)/sgcca_model*.RDS $(FOLDER2)/model*_best.RDS
	
models_report: $(FOLDER2)/models_report.R
	@echo "Reporting model AVEs"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)




# GSV scores for pathways
GSV: intestinal_16S_pathways_metadb/GSVA_scores.R
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Analyse the data with PCA taking into account the categories
Deconvolute: stool_metadb intestinal_RNAseq_metadb intestinal_16S_metadb

upset: upsets.R $(pre-files)
	@echo "Making upset plots"
	R CMD BATCH $(R_OPTS) $(<F)

clean:
	find . -name "*.Rout" -type f -delete
	find . -name "*gsea*.csv" -type f -delete
	find . -name "*enrichment*.csv" -type f -delete
