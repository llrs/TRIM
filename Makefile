R_OPTS=--vanilla

cleaning.R: stools_16S/OTUs_Table-refined-stools.tab stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt
	R CMD BATCH $(R_OPTS) cleaning.R
	
intestinal_16S/otus_coherent.csv: cleaning.R
	R CMD BATCH $(R_OPTS) cleaning.R
	
stools_16S/otus_coherent.csv: cleaning.R
	R CMD BATCH $(R_OPTS) cleaning.R
	
stool_intestinal_integration.R: intestinal_16S/otus_coherent.csv stools_16S/otus_coherent.csv helper_functions.R
	R CMD BATCH $(R_OPTS) stool_intestinal_integration.R
	