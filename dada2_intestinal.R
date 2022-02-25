library("dada2")
library("ggplot2")
library("dplyr")
library("tidyr")

# path <- "/home/lrevilla/Documents/projects/TRIM/intestinal_16S/HSCT/HSCT_biopsies"
# 
# r1 <- file.path(path, "HSCT_biopsies-R1.fastq.tar.gz")
# r2 <- file.path(path, "HSCT_biopsies-R2.fastq.tar.gz")
# i1 <- file.path(path, "HSCT_biopsies-I1.fastq.tar.gz")
# i2 <- file.path(path, "HSCT_biopsies-I2.fastq.tar.gz")

# Uncompress (only need to do it once)
# untar(r1, exdir = path)
# untar(r2, exdir = path)
# untar(i1, exdir = path)
# untar(i2, exdir = path)

# Uncompressed and one sample == one file are found here:
r1 <- list.files("intestinal_16S/Samples", pattern = "@R.fastq", full.names = TRUE)
r2 <- list.files("intestinal_16S/Samples", pattern = "@F.fastq", full.names = TRUE)
r1 <- r1[!grepl("Wasser", r1)]
r2 <- r2[!grepl("Wasser", r2)]

# Check quality if desired
# plotQualityProfile(r1) 
# plotQualityProfile(r2)

tempdir_r1 <- tempdir()
tempdir_r2 <- tempdir()
filtR1 <- paste0(tempdir_r1, "/filt_", basename(r1))
filtR2 <- paste0(tempdir_r2, "/filt_", basename(r2))

FaT <- filterAndTrim(fwd=r1, filt=filtR1, rev=r2, filt.rev=filtR2,
              trimLeft = 20, truncLen= 260, trimRight = 20,
              maxN=0, maxEE=2, matchIDs = TRUE,
              compress=TRUE, verbose=FALSE)

derepR1 <- derepFastq(filtR1, verbose=FALSE)
derepR2 <- derepFastq(filtR2, verbose=FALSE)
errR1 <- learnErrors(derepR1, multithread=FALSE)
errR2 <- learnErrors(derepR2, multithread=FALSE)
plotErrors(errR1, nominalQ = TRUE)
plotErrors(errR2, nominalQ = TRUE)
dadaR1 <- dada(derepR1, err=errR1, multithread=FALSE)
dadaR2 <- dada(derepR2, err=errR2, multithread=FALSE)
merger1 <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=FALSE)
merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=FALSE)
seqtab <- makeSequenceTable(merger1)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=FALSE)
# From https://web.stanford.edu/class/bios221/Pune/Lectures/Lecture_Day1_dada2_workflow.pdf
getN <-function(x){sum(getUniques(x))}
track <-cbind(FaT,
              sapply(dadaR1, getN),
              sapply(dadaR2, getN),
              sapply(merger1, getN),
              rowSums(seqtab.nochim))
colnames(track) <-c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

track %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample") %>% 
  group_by(sample) %>% 
  pivot_longer(cols = input:nonchim) %>% 
  mutate(perc = value/lag(value),
         perc = if_else(is.na(perc), 1, perc)) %>% 
  ungroup() %>% 
  mutate(name = forcats::fct_relevel(name, 
                                     c("input", "filtered", "denoisedR", "denoisedF", "merged", "nonchim"))) %>% 
  filter(!name %in% c("input", "denoisedR", "denoisedF")) %>% 
  ggplot() +
  geom_hline(yintercept = 0.5) +
  geom_violin(aes(name, perc)) +
  geom_point(aes(name, perc, group = sample, col = sample)) +
  geom_line(aes(name, perc, group = sample, col = sample)) +
  guides(col = FALSE) +
  scale_x_discrete(expand = expansion(add = 0.1)) +
  labs(x = element_blank()) +
  theme_minimal()

# saveRDS(seqtab.nochim, "intestinal_16S_RNAseq_metadb/dada2.RDS")
seqtab.nochim <- readRDS("intestinal_16S_RNAseq_metadb/dada2.RDS")
code_samples <- strcapture("filt_([0-9]+)\\.(.+)@", 
                rownames(seqtab.nochim), 
                data.frame(code = numeric(), sample = character()))


ASV <- colnames(seqtab.nochim)
counts_ASV <- seqtab.nochim
colnames(counts_ASV) <- NULL

ASV_counts <- t(counts_ASV)
cASV <- sort(colSums(ASV_counts))
barplot(log10(cASV))
abline(h = c(log10(500), log10(median(cASV))), col = c("red", "green"))

out <- assignTaxonomy(ASV, refFasta = "intestinal_16S/silva_nr99_v138_train_set.fa.gz",
                      outputBootstraps = TRUE,
               tryRC = TRUE, multithread = 6, verbose = TRUE)
summary(out$tax[, "Family"] %in% c(NA, "Mitochondria"))
saveRDS(out, "intestinal_16S/taxonomy_ASV.RDS")
out <- readRDS("intestinal_16S/taxonomy_ASV.RDS")
# taxonomy <- readRDS("data_out/taxonomy_ASV.RDS")
