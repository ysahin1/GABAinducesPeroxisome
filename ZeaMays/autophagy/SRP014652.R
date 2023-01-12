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

autophagy_SRP014652_de <- dn_SRP014652[which(dn_SRP014652 %in% os_autophagy_genes$locus_IDs)]
length(autophagy_SRP014652_de)

sm_atg_SRP014652 <- as.vector(os_autophagy_genes[which(as.vector(os_autophagy_genes$locus_IDs)%in% autophagy_SRP014652_de),])
dim(sm_atg_SRP014652)

autophagy_vstc_SRP014652 <- vstc_dds_SRP014652[autophagy_SRP014652_de,]
DEauto_zscore_SRP014652 <- (autophagy_vstc_SRP014652-rowMeans(autophagy_vstc_SRP014652))/(rowSds(as.matrix(autophagy_vstc_SRP014652)))
head(DEauto_zscore_SRP014652)

SRP014652_coln <- data.frame(sample = c("6_DAP","6_DAP",
                                        "12_DAP","12_DAP",
                                        #"18_DAP","18_DAP",
                                        #"24_DAP",
                                        "30_DAP","30_DAP",
                                        "0_DAP","0_DAP") )
SRP014652_coln
row.names(SRP014652_coln) <- colnames(vstc_dds_SRP014652)

ind_idx <- match(rownames(DEauto_zscore_SRP014652),sm_atg_SRP014652$locus_IDs)
tiff("SRP014652_auto.tiff", units="in", width=6, height=10, res=600)
pheatmap(DEauto_zscore_SRP014652,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = SRP014652_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  paste(sm_atg_SRP014652$gene_symbol[ind_idx],sm_atg_SRP014652$locus_IDs[ind_idx] ,sep = "-"),
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
