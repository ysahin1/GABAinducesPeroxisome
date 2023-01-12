library(genefilter)
library(pheatmap)
autophagy_gene_table <- read.delim("/media/yunus/TOSHIBA1TB/c3-c4/peroxisome-C3vsC4/autophagy/autophagy_homologs.csv", header=T, stringsAsFactors = F)
dim(autophagy_gene_table)
locus_IDs <- strsplit(as.vector(autophagy_gene_table$zea_mays_), "," )
nIDs <- sapply(1:length(locus_IDs), function(x) length(locus_IDs[[x]]))
locus_IDs <- gsub(" ", "", unlist(locus_IDs), fixed = TRUE)
locus_IDs <- strsplit(locus_IDs, "_")
locus_IDs <- sapply(1:length(locus_IDs), function(x) locus_IDs[[x]][[1]] ) 
os_autophagy_genes <- data.frame(gene_symbol = rep(autophagy_gene_table$gene_symbol, nIDs), 
                                 locus_IDs = unlist(locus_IDs), 
                                 stringsAsFactors=FALSE )
os_autophagy_genes <- os_autophagy_genes[!duplicated(os_autophagy_genes$locus_IDs), ]
#####list comes from pythozome, no need to filter duplicates.
#autophagy_gene_table <- autophagy_gene_table[!duplicated(autophagy_gene_table$at_locus_id),]
#autophagy_gene_table$at_locus_id <- toupper(autophagy_gene_table$at_locus_id)

autophagy_SRP014652v_de <- dn_SRP014652v[which(dn_SRP014652v %in% os_autophagy_genes$locus_IDs)]
length(autophagy_SRP014652v_de)

sm_atg_SRP014652v <- as.vector(os_autophagy_genes[which(as.vector(os_autophagy_genes$locus_IDs)%in% autophagy_SRP014652v_de),])
dim(sm_atg_SRP014652v)

autophagy_vstc_SRP014652v <- vstc_dds_SRP014652v[autophagy_SRP014652v_de,]
DEauto_zscore_SRP014652v <- (autophagy_vstc_SRP014652v-rowMeans(autophagy_vstc_SRP014652v))/(rowSds(as.matrix(autophagy_vstc_SRP014652v)))
head(DEauto_zscore_SRP014652v)

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
                                         "thirteenth_leaf_r2") )
SRP014652v_coln
row.names(SRP014652v_coln) <- colnames(vstc_dds_SRP014652v)

ind_idx <- match(rownames(DEauto_zscore_SRP014652v),sm_atg_SRP014652v$locus_IDs)
tiff("SRP014652v_auto.tiff", units="in", width=6, height=10, res=600)
pheatmap(DEauto_zscore_SRP014652v,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = SRP014652v_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  paste(sm_atg_SRP014652v$gene_symbol[ind_idx],sm_atg_SRP014652v$locus_IDs[ind_idx] ,sep = "-"),
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
