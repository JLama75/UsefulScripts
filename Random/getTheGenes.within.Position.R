library(GenomicRanges)
library(dplyr)

download.file(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz",
  destfile = "gencode.v48.annotation.gtf.gz"
)
setwd("~/Documents/transferToNewLaptop/OGI/R_data/ChIP-seq_plot")
gtf <- rtracklayer::import("gencode.v48.annotation.gtf.gz")

# Filter for genes on chr1 within your region
genes_chr1 <- gtf %>%
  subset(type == "gene" &
           seqnames == "chr1" &
           start <= 103408522 &
           end >= 102576473)

unique(genes_chr1$gene_name)

protein_coding <- genes_chr1[genes_chr1$gene_type == "protein_coding", ]
unique(protein_coding$gene_name)

#########################################################################

genes_chr7 <- gtf %>%
  subset(type == "gene" &
           seqnames == "chr7" &
           start <= 139644985 &
           end >= 139423685)

unique(genes_chr7$gene_name)

protein_coding <- genes_chr7[genes_chr7$gene_type == "protein_coding", ]
unique(protein_coding$gene_name)
