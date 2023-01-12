library(pheatmap)
library(DESeq2)
library(dplyr)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv")
perox_GSE155932_de <- dn_GSE155932[which(dn_GSE155932 %in% as.vector(zm_per$os_locus))]
nperDEG <- length(perox_GSE155932_de)
# the number of peroxisome DEGs
nperDEG
# match the symbols
perox_vstc_GSE155932 <- vstc_dds_GSE155932[which( rownames(vstc_dds_GSE155932) %in% perox_GSE155932_de),]

DEperox_zscore_GSE155932 <- (perox_vstc_GSE155932-rowMeans(perox_vstc_GSE155932))/(rowSds(as.matrix(perox_vstc_GSE155932)))
GSE155932_coln <- data.frame(sample = GSE155932_grow_CG)
row.names(GSE155932_coln) <- colnames(vstc_dds_GSE155932)
int_class <- match(rownames(DEperox_zscore_GSE155932),
                   as.vector(zm_per$os_locus))
annotation_row = data.frame(df_Group = factor(as.vector(zm_per$functions[int_class])
                                              ,ordered = F),
                            locus_id = as.vector(zm_per$os_locus[int_class]))
rownames(annotation_row) = rownames(DEperox_zscore_GSE155932)
annotation_row <-  arrange(annotation_row,df_Group)
annotation_row_2 <- data.frame(df_Group = as.vector(annotation_row$df_Group))
rownames(annotation_row_2) <- annotation_row$locus_id
sym_idx <- match(rownames(annotation_row_2),as.vector(zm_per$os_locus))
sm_perox_GSE155932 <- as.vector(zm_per$symbols[sym_idx])
jpeg(filename = "s5_s4.jpeg", units = "in", width = 25, height = 20, res=600)
pheatmap(DEperox_zscore_GSE155932[rownames(annotation_row_2),],
         show_rownames=T, 
         cluster_cols=F, 
         cluster_rows=F, 
         scale="none",
         annotation_col = GSE155932_coln,
         annotation_row = annotation_row_2,
         cex=1, 
         #clustering_distance_rows="euclidean", 
         labels_row =  paste(sm_perox_GSE155932,rownames(annotation_row_2) ,sep = "-"),
         clustering_distance_cols="euclidean", 
         clustering_method="complete", 
         border_color=F, 
         fontsize = 20)
dev.off()

