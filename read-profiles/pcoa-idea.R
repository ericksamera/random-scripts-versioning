library(vegan)
library(ggplot2)
library(dplyr)

setwd("C:/Users/erick/Documents/GitHub/read-profile")

gene = "H19"
df <- read.csv(paste0("out/",gene,".csv"))# %>% filter(sample != "E59")

meta_df <- read.csv("simple_meta.csv")

species_data <- df[, -1]
rownames(species_data) <- df$sample

bc_dist <- vegdist(species_data, method = "bray")

pcoa_res <- cmdscale(bc_dist, eig = TRUE, k = 2)

# Creating a data frame for the PCoA results
pcoa_df <- as.data.frame(pcoa_res$points)
pcoa_df$Sample <- rownames(pcoa_df)

combined_df <- merge(pcoa_df, meta_df, by.x = "Sample", by.y = "sample")

total_variance <- sum(pcoa_res$eig)
variance_explained <- pcoa_res$eig[1:2] / total_variance * 100


ggplot(combined_df, aes(x = V1, y = V2, color = as.factor(coded_phenotype), label = Sample)) +
  geom_point() +
  geom_text(vjust = -1) +
  xlab(paste("PCoA 1 (", sprintf("%.2f", variance_explained[1]), "% variance)", sep = "")) +
  ylab(paste("PCoA 2 (", sprintf("%.2f", variance_explained[2]), "% variance)", sep = "")) +
  ggtitle(paste0("Principal Coordinates Analysis (PCoA) based on Bray-Curtis Dissimilarity (", gene, ")"))

#permanova_result <- adonis2(bc_dist ~ coded_phenotype, data = meta_df)
#print(permanova_result)

