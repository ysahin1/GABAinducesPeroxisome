library(NMF)
library(pheatmap)
library(genefilter)
zm_per <- read.csv("../os_zm_sym_per_homologs.txt", header = T, sep = "\t")
perox_GSE48507_de <- dn_GSE48507[which(dn_GSE48507 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE48507_de)
nperDEG
head(perox_GSE48507_de)
article_genes <-  c(
  #"Zm00001d028834",
  #"Zm00001d002834",
  "Zm00001d036001",
  "Zm00001d027366",
  #"Zm00001d042338",
  #"Zm00001d023583",
  #"Zm00001d047461",
  #"Zm00001d009978",
  "Zm00001d054044",
  "Zm00001d027511",
  "Zm00001d014848"
)

article_perox_genes <- perox_GSE48507_de[match(article_genes,perox_GSE48507_de)]
article_perox_genes
sym_idx <- match(article_perox_genes,zm_per$zm_locus)

sm_perox_GSE48507 <- as.vector(zm_per$symbols[sym_idx])
perox_vstc_GSE48507 <- vstc_dds_GSE48507[article_perox_genes,]
dim(perox_vstc_GSE48507)
DEperox_zscore_GSE48507 <- (perox_vstc_GSE48507-rowMeans(perox_vstc_GSE48507))/(rowSds(as.matrix(perox_vstc_GSE48507)))
head(DEperox_zscore_GSE48507)
GSE48507_coln <- data.frame(sample = c("control", "control","drought_1","drought_1","drought_2","drought_2"))
colnames(DEperox_zscore_GSE48507)
GSE48507_coln
row.names(GSE48507_coln) <- colnames(vstc_dds_GSE48507)
sm_perox_GSE48507
pheatmap(DEperox_zscore_GSE48507,
         show_rownames=T, cluster_cols=T, cluster_rows=F, scale="row",annotation_col = GSE48507_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE48507,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

