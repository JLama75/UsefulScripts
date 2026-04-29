library(readr)
library(dplyr)
#Ploting PCA using QC'ed data before imputation. Then selecting for fame or fame cohort to plot PCA

setwd("/gpfs/fs1/data/Segre_Lab/users/jlama/GSA_new.All_040424/PCA/QC/FAME/FAME.EUR/")
ancestry = read_csv("inferred_ancestry.tsv") #1377

data = read_table("study_pca.eigenvec") #554
colnames(data)[1] <- "FID"

ancestry = ancestry %>% select(FID = IID, predicted_ancestry)
data = merge(ancestry, data, by = "FID")

fame = read_table("/gpfs/fs1/data/Segre_Lab/users/jlama/GSA_new.All_040424/pheno/pheno2_all/GWAS.phenotypeFile/FAME.final.052824.txt", col_names = T)
fame = merge(fame, data, by = "FID") #386

#fame.eur <- fame[fame$predicted_ancestry=="EUR",]#386
#fame.eur <- fame.eur[,c(1,7)]
#setwd("/gpfs/fs1/data/Segre_Lab/users/jlama/GSA_new.All_040424/PCA/QC/FAME/FAME.EUR/")
#write.table(fame.eur, "FAME.eur.sampleID.txt", quote = F, row.names = F, sep = "\t")

library(ggplot2)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (28), hjust = 0.5), 
                      legend.title = element_text( colour="black", size = (22)), 
                      legend.text = element_text( colour="black", size = (20)), 
                      axis.title = element_text( size = (22), colour = "black"),
                      axis.text = element_text( colour = "black", size = (18)),
                      legend.box.margin = margin(0, 0, 0, 10))

p1 = ggplot(fame, aes(x = PC1, y = PC2, color = responder, shape = responder)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(3, 1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(10, 10, 10, 10)) + mynamestheme + labs(color = "responder")
p2 = ggplot(fame, aes(x = PC1, y = PC3, color = responder, shape = responder)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(3, 1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(10, 10, 10, 10)
  ) + mynamestheme + labs(color = "responder")
p3 = ggplot(fame, aes(x = PC2, y = PC3, color = responder, shape = responder)) +
  geom_point(size = 3) + 
  scale_shape_manual(values = c(3, 1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(10, 10, 10, 10)
  ) + mynamestheme + labs(color = "responder")


library(gridExtra)
plot = grid.arrange(p1, p2, p3, ncol = 1)
ggsave("fame.QCedData.554SampleIDs.pdf", plot, width = 8, height = 12, units = "in", dpi = 300)
#ggsave("fame_responder.2.pdf", plot, width = 8, height = 12, units = "in", dpi = 300) #symbol- + and o

