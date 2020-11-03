library("MOSim")

# Load data ####
# Prepare ####
rnaseq_customdata <- omicData("RNA-seq", data = custom_rnaseq)
micro_customdata <- omicData("16S", data = custom_rnaseq)
rnaseq_customdata <- omicData("RNA-seq", data = custom_rnaseq)
# Simulate data according to our own data
data <- mosim(omics = list("RNA-seq" = rnaseq_customdata,
                           "micro" = micro_customdata),
              )
# Not able to simulate 16S
# Evaluate our method ####
# Compare with mcia ####
