# Prepare for the excel for NCBI
library("dplyr")
library("stringr")

meta <- readRDS("intestinal_16S_RNAseq_metadb/meta.RDS")

Samples_data <- meta %>% 
  mutate(
    "Sample name" = paste("Sample", seq_len(n())),
    title = `Sample Name_RNA`,
    "source name" = "biopsy",
    organism = "Homo sapiens",
    "characteristics: sex" = SEX,
    "characteristics: treatment" = Treatment,
    "characteristics: CDAI" = CDAI,
    "characteristics: PCR" = PCR,
    "characteristics: SESCD global" = SESCD_global,
    "characteristics: SESCD local" = SESCD_local,
    "characteristics: Exact location" = Exact_location,
    "characteristics: Time" = Time,
    "characteristics: Patient ID" = Patient_ID,
    "characteristics: Date sample" = DATE_SAMPLE,
    "characteristics: Afected area" = CD_Aftected_area,
    "characteristics: Involved" = Involved_Healthy,
    "characteristics: Active area" = Active_area,
    "characteristics: Transplant" = Transplant,
    "characteristics: Surgery" = Surgery,
    "characteristics: Age diagnostic" = AgeDiag,
    "characteristics: Age sample" = AGE_SAMPLE,
    molecule = "total RNA"
  )

tfastq <- read.table("GEO/transcriptome_fastq.txt", header = FALSE)
mfastq <- read.table("GEO/microbiome_fastq.txt", header = FALSE)

# transcriptome ####

files_transcriptome <- tfastq %>%
  mutate(name = str_replace(V1, "_[0-9]\\.fastq.+", "")) %>%
  select(name) %>%
  unique() %>%
  mutate(
    Left = paste0(name, "_1.fastq.gz"),
    Right = paste0(name, "_2.fastq.gz"),
    name = toupper(name)
  ) %>%
  filter(!grepl("-w", name)) %>%
  mutate(
    name = case_when(
      name == "16-TM29-TTR-CH" ~ "16-TM30-TTR-CH",
      name == "16-TM29-TTR-CIA" ~ "16-TM30-TTR-CIA",
      TRUE ~ name
    )
  ) %>% 
  filter(name  %in% meta$`Sample Name_RNA`) 

# Double check that the original names are  present and that all of them are in
stopifnot(nrow(files_transcriptome) == 158)
stopifnot(nrow(meta) == 158)
ntfastq <- c(files_transcriptome$Left, files_transcriptome$Right)
stopifnot(all(ntfastq %in% tfastq$V1))

# microbiome ####
files_microbiome <- mfastq %>% 
  mutate(name = str_replace(V1, pattern = "@[RF]\\.fastq", replacement = "")) %>% 
  select(name) %>% 
  unique() %>% 
  mutate(Left = paste0(name, "@F.fastq"),
         Right = paste0(name, "@R.fastq"),
         root = str_replace(name, "[0-9]*\\.", "")) %>% 
  filter(root %in% meta$Seq_code_uDNA)



# Double check that all are in.
stopifnot(nrow(files_microbiome) == 158)
stopifnot(all(c(files_microbiome$Left, files_microbiome$Right) %in% mfastq$V1))

# RAW FILES ####
# To plug in any data 
pivot_
mutate("instrument model" = "Illumina HiSeq 4000",
       "single or paire-end" = "paired",
       "filetype" = "fastq",
       "read length" = 255) # Check the size

# Pending the checksum


# Proceessed data files ####
# To plug on the data
mutate("file type" = "csv")