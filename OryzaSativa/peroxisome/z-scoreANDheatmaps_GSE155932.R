library(pheatmap)
library(DESeq2)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/os_zm_sym_per_homologs.txt")
perox_GSE155932_de <- dn_GSE155932[which(dn_GSE155932 %in% zm_per$os_locus)]
nperDEG <- length(perox_GSE155932_de)
# the number of peroxisome DEGs
nperDEG
# match the symbols
sym_idx <- match(perox_GSE155932_de,zm_per$os_locus)
sm_perox_GSE155932 <- as.vector(zm_per$symbols[sym_idx])

perox_vstc_GSE155932 <- vstc_dds_GSE155932[perox_GSE155932_de,]
dim(perox_vstc_GSE155932)
DEperox_zscore_GSE155932 <- (perox_vstc_GSE155932-rowMeans(perox_vstc_GSE155932))/(rowSds(as.matrix(perox_vstc_GSE155932)))
head(DEperox_zscore_GSE155932)


GSE155932_coln <- data.frame(sample = c("s2","s3","s3","s4","s4","s5","s5","s1","s1","s2"))
GSE155932_coln
row.names(GSE155932_coln) <- colnames(vstc_dds_GSE155932)

tiff(filename = "GSE155932_peroxisome.tiff", units = "in", width = 15, height = 25, res=600)
pheatmap(DEperox_zscore_GSE155932,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE155932_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  paste(sm_perox_GSE155932,perox_GSE155932_de ,sep = "-"),
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)
dev.off()
