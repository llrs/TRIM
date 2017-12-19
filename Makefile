R_OPTS=--vanilla
pre_files=stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt \
					stools_16S/OTUs-Table-refined-stools.tab \
					intestinal_16S/OTUs-Table-new-biopsies.csv \
					intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt \
					intestinal_RNAseq/111217_metadata.csv
out_files=meta_coherent.csv \
					stools_16S/otus_coherent.csv \
					stools_16S/taxonomy.csv \
					intestinal_16S/otus_coherent.csv \
					intestinal_16S/taxonomy.csv \
					equivalent_otus.csv \
					equivalent_genus.csv

all: eqGenus eqSpecies \
stool_intestinal_16S_integration/important_common_microrg.csv PCA Deconvolute \
STATegRa intestinal_16S_RNAseq_integration

# Clean the input and prepare the output for integration
$(out_files): cleaning.R $(pre_files) 
	@echo "Preparing input data"
	R CMD BATCH $(R_OPTS) $(<F)
	
#	Handles the use of RGCCA in stools 16S and biopsies 16S
stool_intestinal_16S_integration/*.csv: stool_intestinal_16S_integration/stool_intestinal_integration.R  $(out_files) helper_functions.R
	@echo "Integrating stools and biopsies 16S data" 
	@echo "\tUsing RGCCA"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
#	Handles the use of STATegRa in stools 16S and biopsies 16S
stool_intestinal_16S_integration/STATegRa.Rout: stool_intestinal_16S_integration/STATegRa.R $(out_files) helper_functions.R
	@echo "Integrating stools and biopsies data" 
	@echo "\tUsing STATegRa"
	cd $(<D); R CMD BATCH $(R_OPTS) STATegRa.R
	
STATegRa: stool_intestinal_16S_integration/STATegRa.Rout

# Integrates via correlations the OTUs of the same species
eqSpecies: stool_intestinal_correlation/stool_intestinal_16S_otus.R equivalent_otus.csv stool_intestinal_metadb/important_common_microrg.csv helper_functions.R
	@echo "Analyse the microorganisms in common between stools and biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
# Integrates via correlations the OTUs of the same genus
eqGenus: stool_intestinal_correlation/stool_intestinal_16S_genus.R equivalent_genus.csv stool_intestinal_metadb/important_common_microrg.csv helper_functions.R
	@echo "Analyse the microorganisms in common between stools and biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Code to integrate clinical variables, biopsies 16S and stools 16S
stool_intestinal_metadb/*.csv: stool_intestinal_metadb/stool_intestinal_metadb.R $(out_files)  helper_functions.R
	@echo "Integrating stools and biopsies 16S data \
	taking into account the metadata"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
 	
# Handles the integration of only the ileum between biopsies 16S and stools 16S
ileum_integration: ileum_integration/ileum_integration.R $(out_files) helper_functions.R
	@echo "Integrating stools and biopsies data from ileum"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the integration of only the ileum between biopsies 16S and stools 16S
colon_integration: colon_integration/colon_integration.R $(out_files) helper_functions.R
	@echo "Integrating stools and biopsies data from colon"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCAs	
PCA: PCAs.R $(pre_files) helper_functions.R
	@echo "PCAs of the data"
	R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for stools but controlling for clinical variables
stool_metadb: stool_metadb/stool_metadb.R $(pre_files) helper_functions.R
	@echo "PCAs controlling for clinical variables of 16S stools samples"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for biopsies 16S but controlling for clinical variables
intestinal_16S_metadb: intestinal_16S_metadb/intestinal_16S_metadb.R $(pre_files) helper_functions.R
	@echo "PCAs controlling for clinical variables of 16S biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	

# Handles the creation of PCA for biopsies RNAseq but controlling for clinical variables
intestinal_RNAseq_metadb: intestinal_RNAseq_metadb/intestinal_RNAseq_metadb.R $(pre_files) helper_functions.R
	@echo "PCAs controlling for clinical variables of RNAseq biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
# Handles the integration between the biopsies
intestinal_16S_RNAseq_integration: intestinal_16S_RNAseq_integration/intestinal_16S_RNAseq_integration.R $(pre_files) helper_functions.R
	@echo "Relating the RNAseq with the microbiota"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the integration between the biopsies taking into account the metadata
intestinal_16S_RNAseq_metadb: intestinal_16S_RNAseq_metadb/intestinal_16S_RNAseq_metadb.R $(pre_files) helper_functions.R
	@echo "Relating the RNAseq with the microbiota"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the calculation of the prevalence in intestinal 
intestinal_prevalence: intestinal_16S_conceptual/prevalence.R  helper_functions.R $(pre_files)
	@echo "Biopsies prevalence"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the calculation of the prevalence in stools	
stools_prevalence: stools_16S_conceptual/prevalence.R  helper_functions.R $(pre_files)
	@echo "Stools prevalence"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Do both prevalences
prevalence: intestinal_prevalence stools_prevalence

# Analyse the data with PCA taking into account the categories
Deconvolute: stool_metadb intestinal_RNAseq_metadb intestinal_16S_metadb

clean:
	rm *.Rout