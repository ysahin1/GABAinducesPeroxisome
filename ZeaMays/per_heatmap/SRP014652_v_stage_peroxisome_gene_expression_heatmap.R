library(DESeq2)
library(rafalib)
SRP014652v_raw_counts <- read.table("../counts/SRP014652_v_stage_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
row.names(SRP014652v_raw_counts) <- SRP014652v_raw_counts$Geneid
SRP014652v_raw_counts <- SRP014652v_raw_counts[,-c(1:18)] #c(1:7),10) exclude all information that do not include counts
names(SRP014652v_raw_counts) <- c(
                        "8th_V9_R1",
                        "8th_V9_R2",
                        "11th_V9_R1",
                        "11th_V9_R2",
                        "13th_V9_R1",
                        "13th_V9_R2",
                        "13th_V9_R3",
                        "im_V9_R1",
                        "im_V9_R2",
                        "im_V9_R3",
                        "13th_VT_R1",
                        "13th_VT_R2",
                        "13th_R2_R1",
                        "13th_R2_R2"
                        )
SRP014652v_grow_CG <- c("eighth_leaf_v9",
                        "eighth_leaf_v9",
                        "eleventh_leaf_v9",
                        "eleventh_leaf_v9",
                        "thirteenth_leaf_v9",
                        "thirteenth_leaf_v9",
                        "thirteenth_leaf_v9",
                        "immature_Leaves_v9",
                        "immature_Leaves_v9",
                        "immature_Leaves_v9",
                        "thirteenth_leaf_vT",
                        "thirteenth_leaf_vT",
                        "thirteenth_leaf_r2",
                        "thirteenth_leaf_r2")
SRP014652v_grow_CG <- factor(SRP014652v_grow_CG, levels = c("immature_Leaves_v9",
                                                            "eighth_leaf_v9",
                                                            "eleventh_leaf_v9",
                                                            "thirteenth_leaf_v9",
                                                            "thirteenth_leaf_vT",
                                                            "thirteenth_leaf_r2"))
SRP014652v_sample_info <- data.frame(condition=SRP014652v_grow_CG,
                          row.names=names(SRP014652v_raw_counts))
dds_SRP014652v <- DESeqDataSetFromMatrix(countData = SRP014652v_raw_counts,
                                   colData = SRP014652v_sample_info,
                                   design = ~ condition)

dim(assay(dds_SRP014652v))
dds_SRP014652v <- dds_SRP014652v[rowSums(counts(dds_SRP014652v))> 100,]#eleminate counts that are less than zero
vst_dds_SRP014652v <- vst(dds_SRP014652v, blind = TRUE)
vstc_dds_SRP014652v <- assay(vst_dds_SRP014652v)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv", stringsAsFactors = F)
zm_per$zm_locus
per_genes <-  zm_per$zm_locus
indx <- match(per_genes,rownames(vstc_dds_SRP014652v))
indx <- indx[!is.na(indx)]
length(indx)
perg_vst <- vstc_dds_SRP014652v[indx,]
dim(perg_vst)
indx_symbol <- match(rownames(perg_vst),zm_per$zm_locus)
per_symbol <- as.vector(zm_per$symbols[indx_symbol])
zscore <- (perg_vst-rowMeans(perg_vst))/(rowSds(as.matrix(perg_vst)))
dim(zscore)
library(RColorBrewer)

library(gplots) ##Available from CRAN
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.vector(SRP014652v_grow_CG))]
head(cbind(colnames(zscore),cols))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
jpeg("v_stage_pg_heatmap.jpeg",width = 2, height = 2, res = 600,units = "in")
heatmap.2(zscore, 
          labRow = FALSE,
          trace="none",
          cexRow=1,
          cexCol=1.5,
          lhei = c(2,7),
          lwid = c(1,3),
          keysize = 3,
          ColSideColors=cols,
          col=hmcol,
          margins = c(8,0))
dev.off()
######## hierarchical clustering##########

dis_vstc_dds_SRP014652v <- as.dist(1- cor(vstc_dds_SRP014652v, method ="pearson"))
hc <- hclust(dis_vstc_dds_SRP014652v)
############## DEG analysis ############################
dds_SRP014652v <- DESeq(dds_SRP014652v)
dds_SRP014652v_results <- results(dds_SRP014652v, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP014652v_results)
dn_SRP014652v <- rownames(subset(dds_SRP014652v_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP014652v)
