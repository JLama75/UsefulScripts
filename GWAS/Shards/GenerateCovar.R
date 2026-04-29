library("readr")
library("readxl")
library("dplyr")
library("openxlsx")

data <- read_table("/PHShome/jl2251/Ines/Merge_Data/merged_AMD_GSA_CommVar.fam", col_names = F)
pheno <- read_excel("/PHShome/jl2251/Ines/pheno/AllBatch_AMD_pheno_777samples.xlsx")

colnames(data)[1] <- "FID"
colnames(data)[2] <- "IID"

data <- data %>% dplyr::select(FID, IID)

Batch1 <- read_table("/PHShome/jl2251/Ines/Husain-AMD/Husain-AMD_OmniExpress_Feb2017_Passing_BIRDSUITE_.FHG19.fam", col_names = F)
Batch2 <- read_table("/PHShome/jl2251/Ines/33112425786-6/33112425786-6.fam", col_names = F)

Batch1 <- Batch1 %>% dplyr::select(X1, X2)
colnames(Batch1) <- c("FID", "IID")
Batch2 <- Batch2 %>% dplyr::select(X1, X2)
colnames(Batch2) <- c("FID", "IID")

Batch2$StudyID <- gsub("-", "_", Batch2$IID)           # replace all - with _
#Batch2.subset <- merge(data, Batch2, by.x="FID", by.y = "StudyID")

data$Covariate <- ifelse(data$FID %in% Batch2$StudyID, "Batch2", "Batch1")
data %>% group_by(Covariate) %>% count()

Covar <- data %>% dplyr::select(IID, Covariate)
write.table(Covar, "/PHShome/jl2251/Ines/Shards/batches.tsv", row.names = F, quote = F)

pheno$seq_id <- paste0(pheno$FID, "_", pheno$IID)
pheno <- pheno %>% dplyr::select(-FID, -IID)

data <- data %>%  dplyr::select(FID, Covariate)
Batch2.subset <- merge(data, pheno, by.x="FID", by.y = "StudyID")
nrow(Batch2.subset[!duplicated(Batch2.subset$FID),]) #284
Batch2.subset <- Batch2.subset %>% dplyr::select(-seq_id)

Batch1.subset <- merge(data, pheno, by.x="FID", by.y = "seq_id")
nrow(Batch1.subset[!duplicated(Batch1.subset$FID),]) #492
Batch1.subset <- Batch1.subset %>% dplyr::select(-StudyID)

Batch_pheno <- rbind(Batch1.subset, Batch2.subset)
nrow(Batch_pheno[!duplicated(Batch_pheno$FID),])
# Find columns in Batch1 but not in Batch2
#setdiff(names(Batch1.subset), names(Batch2.subset))

write.xlsx(Batch_pheno, "AllBatch_AMD_pheno_777samples_FID.xlsx")

id_counts <- Batch_pheno %>%
  count(FID)
min(id_counts$n)
max(id_counts$n)

write.table(id_counts, "duplicatePheno_777samples.tsv", quote = F, row.names = F)


