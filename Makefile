R_OPTS=--vanilla

all: cleaning stool_intestinal_integration

# Clean the input and prepare it
cleaning:
	R CMD BATCH $(R_OPTS) cleaning.R

# Code to integrate stools 16S and biopsies 16S
stool_intestinal_integration: cleaning
	R CMD BATCH $(R_OPTS) stool_intestinal_integration.R

# Code to integrate biopsies, biopsies 16S and stools 16S
stool_intestinal_metadb: cleaning
	R CMD BATCH $(R_OPTS) stool_intestinal_metadb.R
 	
clean:
	rm *.Rout