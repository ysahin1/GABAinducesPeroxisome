library(ggplot2)
library(DESeq2)
GSE57950hz_raw_counts <- read.table("../HHZ_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
row.names(GSE57950hz_raw_counts) <- GSE57950hz_raw_counts$Geneid
GSE57950hz_raw_counts <- GSE57950hz_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE57950hz_raw_counts) <- c("GSM1398429",
                                  "GSM1398430",
                                  "GSM1398441",
                                  "GSM1398442")
GSE57950hz_grow_CG <- c("control","control","drought","drought")
GSE57950hz_sample_info <- data.frame(condition=GSE57950hz_grow_CG,
                                     row.names=names(GSE57950hz_raw_counts))
dds_GSE57950hz <- DESeqDataSetFromMatrix(countData = GSE57950hz_raw_counts,
                                         colData = GSE57950hz_sample_info,
                                         design = ~ condition)
dds_GSE57950hz <- dds_GSE57950hz[rowSums(counts(dds_GSE57950hz))> 100,]
dds <- DESeq(dds_GSE57950hz)
vsd <- varianceStabilizingTransformation(dds)
library(genefilter)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q95_wpn <- quantile( rowVars(wpn_vsd), .75)
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
library(tidyverse)
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)



tiff(filename = "Normalized_and_95_quantile_Expression.tiff", 
     units = "in", 
     width = 10, 
     height = 5, 
     res=600)
expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 75 quantile Expression",
    x = "Samples",
    y = "Normalized expression"
  )
dev.off()
input_mat = t(expr_normalized)
library(WGCNA)
allowWGCNAThreads()
powers = c(c(1:20), seq(from = 22, to = 40, by = 2))
sft =pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 5
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed hybrid",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 200,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

