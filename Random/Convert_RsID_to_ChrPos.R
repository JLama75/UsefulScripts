# Install biomaRt if you haven't already
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)
library(readr)

# 1. Connect to the Ensembl SNP database for GRCh38 (the current default)
# Use "hsapiens_snp" for human variant data
snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# 2. Define your list of rsIDs
#rs12125734
#rs28613981
#table <- read_table("Khajawa_rsid.tsv", col_names = T)
table <- read_table("Khajawa_SIOP_rsID.tsv", col_names = T)
#rsids <- c("rs12757325", "rs12135187", "rs12125734", "rs28613981")
rsids <- table$rsID
# 3. Query the database
# 'refsnp_id' = rsID
# 'chr_name' = Chromosome
# 'chrom_start' = Genomic position (start)
results <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'),
                 filters = 'snp_filter',
                 values = rsids, 
                 mart = snp_mart)

# 4. View the results
results

table_pos <- merge(table, results, by.x="rsID", by.y="refsnp_id")
write.table(table_pos, "Khajawa_SIOP_rsID_converted.tsv", quote = F, row.names = F)

