R_OPTS=--vanilla
pre_files=stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt \
					stools_16S/OTUs-Table-refined-stools.tab \
					intestinal_16S/OTUs-Table-new-biopsies.csv \
					intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt
out_files=meta_coherent.csv \
					stools_16S/otus_coherent.csv \
					stools_16S/taxonomy.csv \
					intestinal_16S/otus_coherent.csv \
					intestinal_16S/taxonomy.csv \
					equivalent_otus.csv \
					equivalent_genus.csv

all: eqGenus eqSpecies stool_intestinal_integration/important_common_microrg.csv PCA

# Clean the input and prepare the output for integration
$(out_files): cleaning.R $(pre_files) 
	@echo "Preparing input data"
	R CMD BATCH $(R_OPTS) $(<F)
	
#	Handles the use of RGCCA in stools 16S and biopsies 16S
stool_intestinal_integration/*.csv: stool_intestinal_integration/stool_intestinal_integration.R  $(out_files) helper_functions.R
	@echo "Integrating stools and intestinal data" 
	@echo "\tUsing RGCCA"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
#	Handles the use of STATegRA in stools 16S and biopsies 16S
stool_intestinal_integration/STATegRa.Rout: stool_intestinal_integration/STATegRa.R $(out_files) helper_functions.R
	@echo "Integrating stools and intestinal data" 
	@echo "\tUsing STATegRa"
	cd $(<D); R CMD BATCH $(R_OPTS) STATegRa.R

# Integrates via correlations the OTUs of the same species
eqSpecies: stool_intestinal/stool_intestinal_otus.R equivalent_otus.csv stool_intestinal_metadb/important_common_microrg.csv helper_functions.R
	@echo "Analyse the microorganisms in common between stools and biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
# Integrates via correlations the OTUs of the same genus
eqGenus: stool_intestinal/stool_intestinal_genus.R equivalent_genus.csv stool_intestinal_metadb/important_common_microrg.csv helper_functions.R
	@echo "Analyse the microorganisms in common between stools and biopsies"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Code to integrate clinical variables, biopsies 16S and stools 16S
stool_intestinal_metadb/*.csv: stool_intestinal_metadb/stool_intestinal_metadb.R $(out_files)  helper_functions.R
	@echo "Integrating stools and intestinal data \
	taking into account the metadata"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
 	
# Handles the integration of only the ileum between biopsies 16S and stools 16S
ileum_integration: ileum_integration/ileum_integration.R $(out_files) helper_functions.R
	@echo "Integrating stools and intestinal data from ileum"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the integration of only the ileum between biopsies 16S and stools 16S
colon_integration: colon_integration/colon_integration.R $(out_files) helper_functions.R
	@echo "Integrating stools and intestinal data from colon"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCAs	
PCA: PCAs.R $(pre_files) helper_functions.R
	@echo "PCAs of the data"
	R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for biopsies but controlling for clinical variables
stool_metadb: stool_metadb/stool_metadb.R $(pre_files) helper_functions.R
	@echo "PCAs controlling for clinical variables"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)

# Handles the creation of PCA for stools but controlling for clinical variables
intestinal_metadb: intestinal_metadb/intestinal_metadb.R $(pre_files) helper_functions.R
	@echo "PCAs controlling for clinical variables"
	cd $(<D); R CMD BATCH $(R_OPTS) $(<F)
	
clean:
	rm *.Rout