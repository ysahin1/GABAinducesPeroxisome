library(pheatmap)
library(DESeq2)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/os_zm_sym_per_homologs.txt")
perox_SRP014652v_de <- dn_SRP014652v[which(dn_SRP014652v %in% zm_per$zm_locus)]
nperDEG <- length(perox_SRP014652v_de)
# the number of peroxisome DEGs
nperDEG
# match the symbols
sym_idx <- match(perox_SRP014652v_de,zm_per$zm_locus)
sm_perox_SRP014652v <- as.vector(zm_per$symbols[sym_idx])

perox_vstc_SRP014652v <- vstc_dds_SRP014652v[perox_SRP014652v_de,]
dim(perox_vstc_SRP014652v)
DEperox_zscore_SRP014652v <- (perox_vstc_SRP014652v-rowMeans(perox_vstc_SRP014652v))/(rowSds(as.matrix(perox_vstc_SRP014652v)))
head(DEperox_zscore_SRP014652v)
dim(DEperox_zscore_SRP014652v)
vstc_dds_SRP014652v <- vstc_dds_SRP014652v[,c("SRR1619669","SRR1619670",
                                            "SRR1619633","SRR1619634",
                                            "SRR1619636","SRR1619637",
                                            "SRR1619667","SRR1619668")]
DEperox_zscore_SRP014652v <- DEperox_zscore_SRP014652v[,c("SRR531218",
                                                          "SRR531219",
                                                          "SRR531220",
                                                          "SRR531866",
                                                          "SRR531867",
                                                          "SRR531868",
                                                          "SRR531869",
                                                          "SRR531870",
                                                          "SRR531871",
                                                          "SRR531872",
                                                          "SRR531873",
                                                          "SRR531874",
                                                          "SRR940237",
                                                          "SRR940238",
                                                          "SRR940249",
                                                          "SRR940250",
                                                          "SRR940252",
                                                          "SRR940253",
                                                          "SRR940254",
                                                          "SRR940255",
                                                          "SRR940256",
                                                          "SRR940257",
                                                          "SRR940260",
                                                          "SRR940274",
                                                          "SRR940276",
                                                          "SRR940278"
)]
SRP014652v_coln <- data.frame(sample = c("topmost_leaf_v3",
                                         "topmost_leaf_v3",
                                         "topmost_leaf_v3",
                                         "base_stage2_leaf_v5",
                                         "base_stage2_leaf_v5",
                                         "base_stage2_leaf_v5",
                                         "tip_stage2_leaf_v7",
                                         "tip_stage2_leaf_v7",
                                         "tip_stage2_leaf_v7",
                                         "base_stage2_leaf_v7",
                                         "base_stage2_leaf_v7",
                                         "base_stage2_leaf_v7",
                                         "eighth_leaf_v9",
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
                                         "thirteenth_leaf_r2"))
SRP014652v_coln
row.names(SRP014652v_coln) <- colnames(vstc_dds_SRP014652v)
int_class <- match(rownames(DEperox_zscore_SRP014652v),
                   zm_per$zm_locus)
annotation_row = data.frame(df_Group = factor(as.vector(zm_per$class[int_class]), ordered = TRUE))
rownames(annotation_row) = rownames(DEperox_zscore_SRP014652v)
dim(DEperox_zscore_SRP014652v)
tiff(filename = "SRP014652v_peroxisome_class.tiff", units = "in", width = 15, height = 20, res=600)
pheatmap(DEperox_zscore_SRP014652v,
         show_rownames=T, 
         cluster_cols=F, 
         cluster_rows=F, 
         #scale="row",
         annotation_col = SRP014652v_coln,
         annotation_row = annotation_row,
         cex=1, 
         #clustering_distance_rows="euclidean", 
         labels_row =  paste(sm_perox_SRP014652v,perox_SRP014652v_de ,sep = "-"),
         clustering_distance_cols="euclidean", 
         clustering_method="complete", 
         border_color=F, 
         fontsize = 14)
dev.off()
