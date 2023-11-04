# load R
module purge

module load bioinformatics

module load r/gcc/4.3.1

R





# Load necessary libraries
library(DESeq2)
library(tximport)
library(readr)





quant_files <- c(
    "P1_H.m_15-23_young_fruits_quant.sf",
    "P1_H.m_21-29_mature_fruits_quant.sf",
    "P2-h.m_4-10_young_fruits_quant.sf",
    "P2-h.m_6-12_smooth_tissue_quant.sf",
    "P3_h.m_32-40_roots_quant.sf",
    "P1_H.m_16-24_warty_tissue_quant.sf",
    "P1_H.m_24-32_smooth_tissue_quant.sf",
    "P2-h.m_44-51_warty_tissue_quant.sf",
    "P3_h.m_28-36_tuber_tissue_quant.sf",
    "P3_h.m_33-41_young_fruits_quant.sf",
    "P1_H.m_18-26_flowers_quant.sf",
    "P1_H.m_26-34_tuber_tissue_quant.sf",
    "P2-h.m_45-54_tuber_tissue_quant.sf",
    "P3_h.m_29-37_warty_tissue_quant.sf",
    "P3_h.m_34-42_mature_fruits_quant.sf",
    "P1_H.m_19-27_roots_quant.sf",
    "P2-h.m_10-17_flowers_quant.sf",
    "P2-h.m_46-53_leaf_quant.sf",
    "P3_h.m_30-38_smooth_tissue_quant.sf",
    "P3_h.m_35-43_leaf_quant.sf",
    "P1_H.m_1_leaf_quant.sf",
    "P2-h.m_11-18_roots_quant.sf",
    "P2-h.m_5-11_mature_fruits_quant.sf",
    "P3_h.m_31-39_flowers_quant.sf"
)



sampleinfo <- read_csv("sampleinfo_hydno_phytum.csv")


tx2gene <- readRDS("tx2gene_hydno_mosel.RDS")

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)



dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design = ~ Tissue)

# Perform DESeq2 analysis
dds <- DESeq(dds)

vsd <- vst(dds,blind=TRUE)
plotPCA(vsd,intgroup="Tissue")

pca_plot <- plotPCA(vsd, intgroup = "Tissue")

# Add sample labels to the PCA plot
pca_plot <- pca_plot + geom_text(aes(label = quant_files), hjust = 0, vjust = 0, nudge_x = 0.1, nudge_y = 0.1)

# Print the PCA plot
print(pca_plot)
