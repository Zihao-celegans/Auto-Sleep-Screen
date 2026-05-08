library(plyr)
library(dplyr)
library(stringr)
library(SKAT)
library(fdrtool)

# Run this script from the SKAT_MMP/ root directory: Rscript scripts/run_pipeline.R
data_dir <- "inputs"
out_dir  <- "outputs"

# ── Step 1: Generate SSID ────────────────────────────────────────────────────
cat("Generating SSID...\n")
gene_variants <- read.table(file.path(data_dir, "gene_variants.txt"), stringsAsFactors = FALSE)
pattern <- "SN=[a-zA-Z0-9.]{1,}|CODING=[a-zA-Z0-9.]{1,}"
gene <- str_extract(gene_variants$V8, pattern)
gene <- sub("SN=|CODING=", "", gene)
SSID <- subset(data.frame(gene, gene_variants$V3), gene != "NA")
write.table(SSID, file.path(data_dir, "MMP.SSID"),
            row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)

# ── Step 2: Update MMP.fam phenotypes ───────────────────────────────────────
cat("Updating MMP.fam phenotypes...\n")
pheno_file <- read.csv(file.path(data_dir, "combined_phenotype.csv"),
                       header = TRUE, stringsAsFactors = FALSE)
fam_file <- read.table(file.path(data_dir, "MMP.fam"), stringsAsFactors = FALSE)
colnames(fam_file) <- c("Strain", "IID", "PID", "MID", "SEX", "PHENO")
fam_w_phenos <- left_join(fam_file, pheno_file, by = "Strain")
fam_w_phenos$PHENO <- fam_w_phenos$Phenotype
fam_w_phenos$Phenotype <- NULL
write.table(fam_w_phenos, file.path(data_dir, "MMP.fam"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)

# ── Step 3: Full SKAT run ────────────────────────────────────────────────────
cat("Generating SSD (full)...\n")
ssd_full  <- file.path(data_dir, "MMP.SSD")
info_full <- file.path(data_dir, "MMP.info")
if (file.exists(ssd_full))  file.remove(ssd_full)
if (file.exists(info_full)) file.remove(info_full)
Generate_SSD_SetID(
  file.path(data_dir, "MMP.bed"),
  file.path(data_dir, "MMP.bim"),
  file.path(data_dir, "MMP.fam"),
  file.path(data_dir, "MMP.SSID"),
  ssd_full, info_full
)

fam_final   <- read.table(file.path(data_dir, "MMP.fam"))
SSD.info    <- Open_SSD(ssd_full, info_full)
set.seed(100)
Null_Model  <- SKAT_Null_Model(fam_final$V6 ~ 1, out_type = "C")

cat("Running SKAT (full)...\n")
All_SKAT_Data <- SKAT.SSD.All(SSD.INFO = SSD.info, obj = Null_Model)

write.table(All_SKAT_Data$results,
            file.path(out_dir, "SKAT_all-pvals.results"),
            row.names = FALSE, col.names = TRUE, quote = FALSE)

qvals <- fdrtool(All_SKAT_Data$results$P.value,
                 statistic = "pvalue", cutoff.method = "fndr", plot = FALSE)
All_SKAT_Data$results$Q.value <- qvals$qval
All_SKAT_Data$results <- All_SKAT_Data$results[order(All_SKAT_Data$results$P.value), ]
write.table(All_SKAT_Data$results,
            file.path(out_dir, "SKAT_all-qvals.results"),
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# ── Step 4: Reduced SKAT run (genes with ≥4 variants) ───────────────────────
cat("Filtering to genes with >= 4 variants...\n")
SSID_full <- read.table(file.path(data_dir, "MMP.SSID"), stringsAsFactors = FALSE)
colnames(SSID_full) <- c("gene", "variant")
SSID_count <- plyr::count(SSID_full$gene)
genes_to_include <- data.frame(gene = SSID_count$x[SSID_count$freq >= 4])
reduced_SSID <- semi_join(SSID_full, genes_to_include, by = "gene")
cat("Gene sets in reduced SSID:", nrow(genes_to_include), "\n")
write.table(reduced_SSID, file.path(data_dir, "MMP-reduced.SSID"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Generating SSD (reduced)...\n")
ssd_red  <- file.path(data_dir, "MMP-reduced.SSD")
info_red <- file.path(data_dir, "MMP-reduced.info")
if (file.exists(ssd_red))  file.remove(ssd_red)
if (file.exists(info_red)) file.remove(info_red)
Generate_SSD_SetID(
  file.path(data_dir, "MMP.bed"),
  file.path(data_dir, "MMP.bim"),
  file.path(data_dir, "MMP.fam"),
  file.path(data_dir, "MMP-reduced.SSID"),
  ssd_red, info_red
)

SSD.info_reduced <- Open_SSD(ssd_red, info_red)
cat("Running SKAT (reduced)...\n")
All_SKAT_Data_reduced <- SKAT.SSD.All(SSD.INFO = SSD.info_reduced, obj = Null_Model)

qvals_red <- fdrtool(All_SKAT_Data_reduced$results$P.value,
                     statistic = "pvalue", cutoff.method = "fndr", plot = FALSE)
All_SKAT_Data_reduced$results$Q.value <- qvals_red$qval
All_SKAT_Data_reduced$results <- All_SKAT_Data_reduced$results[
  order(All_SKAT_Data_reduced$results$P.value), ]

write.table(All_SKAT_Data_reduced$results,
            file.path(out_dir, "SKAT_all_reduced_940.results"),
            row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Done. Results saved to", out_dir, "\n")
cat("Top 10 genes (reduced run):\n")
print(head(All_SKAT_Data_reduced$results[, c("SetID", "P.value", "Q.value")], 10))
