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
					equivalent_otus.csv

all: stool_intestinal_integration stool_intestinal_metadb

# Clean the input and prepare it
$(out_files): $(pre_files) cleaning.R
	@echo "Preparing input data"
	R CMD BATCH $(R_OPTS) cleaning.R
	
# Code to integrate stools 16S and biopsies 16S
stool_intestinal_integration: stool_intestinal_integration.R $(out_files)
	@echo "Integrating stools and intestinal data" 
	R CMD BATCH $(R_OPTS) $(<F) stool_intestinal_integration/stool_intestinal_integration.Rout 

# Code to integrate biopsies, biopsies 16S and stools 16S
stool_intestinal_metadb: stool_intestinal_metadb.R $(out_files) 
	@echo "Integrating stools and intestinal data \
	taking into account the metadata"
	R CMD BATCH $(R_OPTS) $(<F) stool_intestinal_metadb/stool_intestinal_metadb.Rout 
 	
ileum_integration: ileum_integration.R $(out_files) 
	@echo "Integrating stools and intestinal data from ileum"
	R CMD BATCH $(R_OPTS) $(<F) ileum_integration/ileum_integration.Rout 

colon_integration: colon_integration.R $(out_files) 
	@echo "Integrating stools and intestinal data from colon"
	R CMD BATCH $(R_OPTS) $(<F) colon_integration/colon_integration.Rout 
	
PCA: PCAs.R $(pre_files)
	@echo "PCAs of the data"
	R CMD BATCH $(R_OPTS) $(<F)

clean:
	rm *.Rout