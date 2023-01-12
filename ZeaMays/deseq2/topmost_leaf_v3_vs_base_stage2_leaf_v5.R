library(pheatmap)
library(DESeq2)
library(dplyr)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv")
perox_SRP014652v_de <- dn_SRP014652v[which(dn_SRP014652v %in% zm_per$zm_locus)]
nperDEG <- length(perox_SRP014652v_de)
# the number of peroxisome DEGs
nperDEG

perox_vstc_SRP014652v <- vstc_dds_SRP014652v[perox_SRP014652v_de,]
DEperox_zscore_SRP014652v <- (perox_vstc_SRP014652v-rowMeans(perox_vstc_SRP014652v))/(rowSds(as.matrix(perox_vstc_SRP014652v)))
SRP014652v_coln <- data.frame(sample =  c("topmost_leaf_v3",
                                          "topmost_leaf_v3",
                                          "topmost_leaf_v3",
                                          "base_stage2_leaf_v5",
                                          "base_stage2_leaf_v5",
                                          "base_stage2_leaf_v5"))
row.names(SRP014652v_coln) <- colnames(vstc_dds_SRP014652v)
int_class <- match(rownames(DEperox_zscore_SRP014652v),
                   zm_per$zm_locus)
annotation_row = data.frame(df_Group = factor(as.vector(zm_per$functions[int_class]),
                                              ordered = F))
annotation_row$locus_id <- as.vector(zm_per$zm_locus[int_class])
annotation_row <-  arrange(annotation_row,df_Group)
annotation_row_2 <- data.frame(df_Group = annotation_row$df_Group)
rownames(annotation_row_2) <- annotation_row$locus_id
dim(DEperox_zscore_SRP014652v)

sym_idx <- match(rownames(annotation_row_2),zm_per$zm_locus)
sm_perox_SRP014652v <- as.vector(zm_per$symbols[sym_idx])
tiff(filename = "SRP014652v_peroxisome_function_topmost_leaf_v3_vs_base_stage2_leaf_v5.tiff", 
     units = "in", width = 22, height = 20, res=600)
pheatmap(DEperox_zscore_SRP014652v[rownames(annotation_row_2),],
         show_rownames=T, 
         cluster_cols=F, 
         cluster_rows=F, 
         scale="none",
         annotation_col = SRP014652v_coln,
         annotation_row = annotation_row_2,
         cex=1, 
         #clustering_distance_rows="euclidean", 
         labels_row =  paste(sm_perox_SRP014652v,rownames(annotation_row_2) ,sep = "-"),
         clustering_distance_cols="euclidean", 
         clustering_method="complete", 
         border_color=F, 
         fontsize = 14)
dev.off()
