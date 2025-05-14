# This script asess data quality by sample clustering and visualization
# Unsupervised Analyses
# Know how similar biological replicates are to each other and identifying outliers
# Updated 2023-12-10

# Get the necessary libraries:
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library('dplyr')
library(ggsci) # Color palettes inspired by scientific journals

# Set the run IO:
io <- list(
  'counts' = file.path('.', 'rdata', '2581-all-lines-counts.data'),
  'metadata' = file.path('.', 'rdata', '2581-all-lines-metadata.data'),
  'results.dir' = file.path('.', 'results')
)

# Read in the counts & metadata:
counts <- readRDS(io$counts)
metadata <- readRDS(io$metadata)

# Removing A2058 C18 & C21 clones (truncated proteins in WB):
valid.samples <- rownames(metadata[!(metadata$clone %in% c("C18", "C21")),])
counts <- counts[, valid.samples]
metadata <- metadata[valid.samples, ]

# Run DESeq2
# DeSeq will normalize the counts 
# DESeq assumes the input is raw counts
# It will internally correct for the library size when performing the analysis

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ line + condition 
)
# ~ line + condition:
# test for differential expression between conditions (Control and SIN3B KO) and cell lines

# Creating a DESeq object:
dds <- DESeq(dds)
# nrow(dds)
# nrow(counts)
# Number of rows = Number of genes. Should be equal in counts and dds object

# Hierarchical clustering with correlation heatmaps
# Assess the similarity in gene expression between the different samples in a dataset

# Log2 transform counts
vsd <- vst(dds, blind = TRUE) 

# PCA (Principal Component Analysis)
pcaData <- plotPCA(vsd, intgroup = c("line", "condition"), returnData = TRUE)
# returnData function return the data.frame of PC1 and PC2 with intgroup covariates
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot using `ggplot()` function and save to an object
annotated_pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = line)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()+
  scale_color_lancet()+
  ggtitle("PCA RNA Seq SIN3B KO vs Control")+
  theme(plot.title = element_text(hjust = 0.5))
annotated_pca_plot
