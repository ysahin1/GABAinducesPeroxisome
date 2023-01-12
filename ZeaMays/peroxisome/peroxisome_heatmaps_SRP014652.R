library(pheatmap)
library(DESeq2)
library(dplyr)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv")
perox_SRP014652_de <- dn_SRP014652[which(dn_SRP014652 %in% as.vector(zm_per$zm_locus))]
nperDEG <- length(perox_SRP014652_de)
# the number of peroxisome DEGs
nperDEG# match the symbols
perox_vstc_SRP014652 <- vstc_dds_SRP014652[perox_SRP014652_de,]

DEperox_zscore_SRP014652 <- (perox_vstc_SRP014652-rowMeans(perox_vstc_SRP014652))/(rowSds(as.matrix(perox_vstc_SRP014652)))
SRP014652_coln <- data.frame(sample = SRP014652_grow_CG)
row.names(SRP014652_coln) <- colnames(vstc_dds_SRP014652)
int_class <- match(rownames(DEperox_zscore_SRP014652),
                   as.vector(zm_per$zm_locus))
annotation_row = data.frame(df_Group = factor(as.vector(zm_per$functions[int_class])
                                              ,ordered = F),
                            locus_id = as.vector(zm_per$zm_locus[int_class]))
rownames(annotation_row) = rownames(DEperox_zscore_SRP014652)
annotation_row <-  arrange(annotation_row,df_Group)
annotation_row_2 <- data.frame(df_Group = as.vector(annotation_row$df_Group))
rownames(annotation_row_2) <- annotation_row$locus_id
sym_idx <- match(rownames(annotation_row_2),as.vector(zm_per$zm_locus))
sm_perox_SRP014652 <- as.vector(zm_per$symbols[sym_idx])
jpeg(filename = "30_vs_12.jpeg", units = "in", width = 25, height = 18, res=600)
pheatmap(DEperox_zscore_SRP014652[rownames(annotation_row_2),],
         show_rownames=T, 
         cluster_cols=F, 
         cluster_rows=F, 
         scale="none",
         annotation_col = SRP014652_coln,
         annotation_row = annotation_row_2,
         cex=1, 
         #clustering_distance_rows="euclidean", 
         labels_row =  paste(sm_perox_SRP014652,rownames(annotation_row_2) ,sep = "-"),
         clustering_distance_cols="euclidean", 
         clustering_method="complete", 
         border_color=F, 
         fontsize = 20)
dev.off()
