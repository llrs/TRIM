R_OPTS=--vanilla

cleaning.R: stools_16S/OTUs_Table-refined-stools.tab stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt intestinal_16S/OTUs-Table-new-biopsies.csv
	R CMD BATCH $(R_OPTS) cleaning.R
	
intestinal_16S/otus_coherent.csv stools_16S/otus_coherent.csv: cleaning.R
	R CMD BATCH $(R_OPTS) cleaning.R
	
stool_intestinal_integration.R: intestinal_16S/otus_coherent.csv stools_16S/otus_coherent.csv helper_functions.R meta_coherent.csv stools_16S/taxonomy.csv intestinal_16S/taxonomy.csv
	R CMD BATCH $(R_OPTS) stool_intestinal_integration.R

stool_intestinal_metadb.R: intestinal_16S/otus_coherent.csv stools_16S/otus_coherent.csv helper_functions.R meta_coherent.csv stools_16S/taxonomy.csv intestinal_16S/taxonomy.csv
	R CMD BATCH $(R_OPTS) stool_intestinal_metadb.R