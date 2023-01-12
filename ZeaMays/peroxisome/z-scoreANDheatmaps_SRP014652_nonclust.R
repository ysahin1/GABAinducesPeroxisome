library(pheatmap)
library(DESeq2)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/os_zm_sym_per_homologs.txt")
perox_SRP014652_de <- dn_SRP014652[which(dn_SRP014652 %in% zm_per$zm_locus)]
nperDEG <- length(perox_SRP014652_de)
# the number of peroxisome DEGs
nperDEG
# match the symbols
sym_idx <- match(perox_SRP014652_de,zm_per$zm_locus)
sm_perox_SRP014652 <- as.vector(zm_per$symbols[sym_idx])

perox_vstc_SRP014652 <- vstc_dds_SRP014652[perox_SRP014652_de,]
dim(perox_vstc_SRP014652)
DEperox_zscore_SRP014652 <- (perox_vstc_SRP014652-rowMeans(perox_vstc_SRP014652))/(rowSds(as.matrix(perox_vstc_SRP014652)))
head(DEperox_zscore_SRP014652)
dim(DEperox_zscore_SRP014652)
vstc_dds_SRP014652 <- vstc_dds_SRP014652[,c("SRR1619669","SRR1619670",
                                            "SRR1619633","SRR1619634",
                                            "SRR1619636","SRR1619637",
                                            "SRR1619667","SRR1619668")]
DEperox_zscore_SRP014652 <- DEperox_zscore_SRP014652[,c("SRR1619669","SRR1619670",
                                                        "SRR1619633","SRR1619634",
                                                        "SRR1619636","SRR1619637",
                                                        "SRR1619667","SRR1619668")]
SRP014652_coln <- data.frame(sample = c("0_DAP","0_DAP",
                                        "6_DAP","6_DAP",
                                        "12_DAP","12_DAP",
                                        #"18_DAP","18_DAP",
                                        #"24_DAP",
                                        "30_DAP","30_DAP"))
SRP014652_coln
row.names(SRP014652_coln) <- colnames(vstc_dds_SRP014652)
int_class <- match(rownames(DEperox_zscore_SRP014652),
                   zm_per$zm_locus)
annotation_row = data.frame(df_Group = factor(as.vector(zm_per$class[int_class]), ordered = TRUE))
rownames(annotation_row) = rownames(DEperox_zscore_SRP014652)
dim(DEperox_zscore_SRP014652)
tiff(filename = "SRP014652_peroxisome_class.tiff", units = "in", width = 13, height = 16, res=600)
pheatmap(DEperox_zscore_SRP014652,
         show_rownames=T, 
         cluster_cols=F, 
         cluster_rows=F, 
         #scale="row",
         annotation_col = SRP014652_coln,
         annotation_row = annotation_row,
         cex=1, 
         #clustering_distance_rows="euclidean", 
         labels_row =  paste(sm_perox_SRP014652,perox_SRP014652_de ,sep = "-"),
         clustering_distance_cols="euclidean", 
         clustering_method="complete", 
         border_color=F, 
         fontsize = 14)
dev.off()
