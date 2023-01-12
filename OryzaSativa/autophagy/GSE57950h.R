library(genefilter)
library(pheatmap)
autophagy_gene_table <- read.delim("/media/yunus/TOSHIBA1TB/c3-c4/peroxisome-C3vsC4/autophagy/autophagy_homologs.csv", header=T, stringsAsFactors = F)
dim(autophagy_gene_table)
locus_IDs <- strsplit(as.vector(autophagy_gene_table$rice_autophagy_pythozome), "," )
nIDs <- sapply(1:length(locus_IDs), function(x) length(locus_IDs[[x]]))
os_autophagy_genes <- data.frame(gene_symbol = rep(autophagy_gene_table$gene_symbol, nIDs), 
                           locus_IDs = unlist(locus_IDs), 
                           stringsAsFactors=FALSE )
os_autophagy_genes
#####list comes from pythozome, no need to filter duplicates.
#autophagy_gene_table <- autophagy_gene_table[!duplicated(autophagy_gene_table$at_locus_id),]
#autophagy_gene_table$at_locus_id <- toupper(autophagy_gene_table$at_locus_id)

autophagy_GSE155932_de <- dn_GSE155932[which(dn_GSE155932 %in% os_autophagy_genes$locus_IDs)]
length(autophagy_GSE155932_de)

sm_atg_GSE155932 <- as.vector(os_autophagy_genes[which(as.vector(os_autophagy_genes$locus_IDs)%in% autophagy_GSE155932_de),])
dim(sm_atg_GSE155932)

autophagy_vstc_GSE155932 <- vstc_dds_GSE155932[autophagy_GSE155932_de,]
DEauto_zscore_GSE155932 <- (autophagy_vstc_GSE155932-rowMeans(autophagy_vstc_GSE155932))/(rowSds(as.matrix(autophagy_vstc_GSE155932)))
head(DEauto_zscore_GSE155932)

GSE155932_coln <- data.frame(sample = c("s2","s3","s3","s4","s4","s5","s5","s1","s1","s2"))
GSE155932_coln
row.names(GSE155932_coln) <- colnames(vstc_dds_GSE155932)

ind_idx <- match(rownames(DEauto_zscore_GSE155932),sm_atg_GSE155932$locus_IDs)
tiff("GSE155932_auto.tiff", units="in", width=6, height=10, res=600)
pheatmap(DEauto_zscore_GSE155932,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE155932_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_atg_GSE155932$gene_symbol[ind_idx],
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
  