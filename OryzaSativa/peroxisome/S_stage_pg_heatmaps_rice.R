library(DESeq2)

GSE155932_raw_counts <- read.table("../counts/GSE155932_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
row.names(GSE155932_raw_counts) <- GSE155932_raw_counts$Geneid
GSE155932_raw_counts <- GSE155932_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE155932_raw_counts) <- c(
  "S2_R1",
  "S3_R1",
  "S3_R2",
  "S4_R1",
  "S4_R2",
  "S5_R1",
  "S5_R2",
  "S1_R1",
  "S1_R2",
  "S2_R2")
GSE155932_grow_CG <- c("s2",
                       "s3","s3",
                       "s4","s4",
                       "s5","s5",
                       "s1","s1",
                       "s2")
GSE155932_grow_CG <- factor(GSE155932_grow_CG, levels = c("s1",
                                                            "s2",
                                                            "s3",
                                                            "s4",
                                                            "s5"))
GSE155932_sample_info <- data.frame(condition=GSE155932_grow_CG,
                                     row.names=names(GSE155932_raw_counts))
dds_GSE155932 <- DESeqDataSetFromMatrix(countData = GSE155932_raw_counts,
                                         colData = GSE155932_sample_info,
                                         design = ~ condition)

dim(assay(dds_GSE155932))
dds_GSE155932 <- dds_GSE155932[rowSums(counts(dds_GSE155932))> 100,]#eleminate counts that are less than zero
vst_dds_GSE155932 <- vst(dds_GSE155932, blind = TRUE)
vstc_dds_GSE155932 <- assay(vst_dds_GSE155932)

rice_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv", stringsAsFactors = F)
per_genes <-  rice_per$os_locus
indx <- match(per_genes,rownames(vstc_dds_GSE155932))
indx <- indx[!is.na(indx)]
length(indx)
perg_vst <- vstc_dds_GSE155932[indx,]
dim(perg_vst)
indx_symbol <- match(rownames(perg_vst),rice_per$os_locus)
per_symbol <- as.vector(rice_per$symbols[indx_symbol])
zscore <- (perg_vst-rowMeans(perg_vst))/(rowSds(as.matrix(perg_vst)))



for (i in 1:ncol(zscore)) {
  n = 0
  
  for (ii in 1:length(zscore[,i])){
    if (ii > n) {
      
      rownames(zscore[ii,i])
      
    }
  }
}  










library(RColorBrewer)
library(rafalib)
library(gplots) ##Available from CRAN
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.vector(GSE155932_grow_CG))]
rows <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.vector(rice_per$functions[indx_symbol]))]
head(cbind(colnames(zscore),cols))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
jpeg("deneme.jpeg",width = 2, height = 2.5, res = 600,units = "in")
heatmap.2(zscore, 
          labRow = T,
          trace="none",
          cexCol=1.5,
          lhei = c(4,8),
          lwid = c(2,5),
          keysize = 3,
          ColSideColors=cols, 
          RowSideColors = rows,
          colRow=cols,
          col=hmcol,
          margins = c(5,5))
dev.off()
######## hierarchical clustering##########
annotation_row = data.frame(df_Group = factor(as.vector(rice_per$functions[indx_symbol])
                                              ,ordered = F),
                            locus_id = as.vector(rice_per$os_locus[indx_symbol]))
rownames(annotation_row) = rownames(zscore)
annotation_row <-  arrange(annotation_row,df_Group)
annotation_row_2 <- data.frame(df_Group = as.vector(annotation_row$df_Group))
rownames(annotation_row_2) <- annotation_row$locus_id

GSE155932_coln <- data.frame(sample = GSE155932_grow_CG)
row.names(GSE155932_coln) <- colnames(perg_vst)

pheatmap(zscore,
         show_rownames=T, 
         cluster_cols=T, 
         cluster_rows=T, 
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
dis_vstc_dds_GSE155932 <- as.dist(1- cor(vstc_dds_GSE155932, method ="pearson"))
hc <- hclust(dis_vstc_dds_GSE155932)
############## DEG analysis ############################
dds_GSE155932 <- DESeq(dds_GSE155932)
dds_GSE155932_results <- results(dds_GSE155932, 
                                  independentFiltering = TRUE, 
                                  alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE155932_results)
dn_GSE155932 <- rownames(subset(dds_GSE155932_results, 
                                 padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE155932)
