library(DESeq2)
library(rafalib)
library(ggplot2)
library(plyr)
library(ggsci)
library(RColorBrewer)
SRP014652_raw_counts <- read.table("../counts/SRP014652_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
row.names(SRP014652_raw_counts) <- SRP014652_raw_counts$Geneid
SRP014652_raw_counts <- SRP014652_raw_counts[,-c(1:6, 11, 12, 13, 18)] #c(1:7),10) exclude all information that do not include counts
names(SRP014652_raw_counts) <- c(
                        "6_DAP_R1",
                        "6_DAP_R2",
                        "12_DAP_R1",
                        "12_DAP_R2",
                        "30_DAP_R1",
                        "30_DAP_R2",
                        "0_DAP_R1",
                        "0_DAP_R2")
SRP014652_grow_CG <- c("6_DAP","6_DAP",
                        "12_DAP","12_DAP",
                        "30_DAP","30_DAP",
                        "0_DAP","0_DAP")
SRP014652_grow_CG <- factor(SRP014652_grow_CG, levels = c("0_DAP",
                                                            "6_DAP",
                                                            "12_DAP",
                                                            "30_DAP"))
SRP014652_sample_info <- data.frame(condition=SRP014652_grow_CG,
                          row.names=names(SRP014652_raw_counts))
dds_SRP014652 <- DESeqDataSetFromMatrix(countData = SRP014652_raw_counts,
                                   colData = SRP014652_sample_info,
                                   design = ~ condition)
dds_SRP014652 <- dds_SRP014652[rowSums(counts(dds_SRP014652))> 100,]#eleminate counts that are less than zero
vst_dds_SRP014652 <- vst(dds_SRP014652, blind = TRUE)
vstc_dds_SRP014652 <- assay(vst_dds_SRP014652)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv", stringsAsFactors = F)
per_genes <-  zm_per$zm_locus
indx <- match(per_genes,rownames(vstc_dds_SRP014652))
indx <- indx[!is.na(indx)]
perg_vst <- vstc_dds_SRP014652[indx,]
indx_symbol <- match(rownames(perg_vst),zm_per$zm_locus)
per_symbol <- as.vector(zm_per$symbols[indx_symbol])
zscore <- (perg_vst-rowMeans(perg_vst))/(rowSds(as.matrix(perg_vst)))
fr_df <- data.frame()
for(i in 1:ncol(zscore)){
  ex = 0
  cn <- colnames(zscore)[i]
  for (ii in rownames(zscore)) {
    ge <- zscore[ii,i]
    print(ii)
    if (ge > ex) {
      gi <- match(ii,zm_per$zm_locus)
      cl <- zm_per$functions[gi]
      fr_df[ii,cn] <- cl
    }
    else (fr_df[ii,cn] <- "DOWN REGULATED")
  }
}
gdf <- data.frame()
for (i in 1:ncol(fr_df)){
  print(i)
  fr_tb <- as.data.frame(sort(table(fr_df[,i]),decreasing = TRUE))
  fr_tb <- fr_tb[-1,]
  fr_tb$sample <- colnames(fr_df)[i]
  gdf <- rbind(gdf,fr_tb)
}

df_sorted <- arrange(gdf, sample, Var1) 
df_cumsum <- ddply(df_sorted, "sample",
                   transform, label_ypos=cumsum(Freq))
mycolors <- as.vector(rbind(pal_nejm("default", alpha = 0.99)(8), 
                            pal_npg("nrc", alpha = 0.99)(8),
                            pal_jco("default", alpha = 0.99)(4)))
p1 <- ggplot(data=df_cumsum, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=0.4) +
  geom_text(aes(label=Freq), 
            #vjust=4,
            #hjust = "inward",
              position = position_stack(vjust = 0.001),
            color="black", 
            size=3.5)+
  scale_fill_manual(values = mycolors)+
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",
                                 size = 0.5, 
                                 linetype = "solid")
        )

p1
  library(ggpp)
leg <- get_legend(p1)
as_ggplot(leg)

scale_x_discrete(limits=c("0_DAP_R1",
                          "0_DAP_R2",
                          "6_DAP_R1",
                          "6_DAP_R2"))+

ggplot(data=gdf, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=0.4) +
  geom_text(aes(label=Freq), position = position_dodge(0.9), vjust=1.6, color="white", size=3.5)+
  theme_minimal()
library(gplots) ##Available from CRAN
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.vector(SRP014652_grow_CG))]
head(cbind(colnames(zscore),cols))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
jpeg("DPA_stage_pg_heatmap.jpeg",width = 2, height = 2, res = 600,units = "in")
heatmap.2(zscore, 
          labRow = FALSE,
          trace="none",
          cexRow=1,
          cexCol=1.5,
          lhei = c(2,7),
          lwid = c(1,3),
          keysize = 3,
          ColSideColors=cols,
          col=hmcol,
          margins = c(8,0))
dev.off()
######## hierarchical clustering##########

dis_vstc_dds_SRP014652 <- as.dist(1- cor(vstc_dds_SRP014652, method ="pearson"))
hc <- hclust(dis_vstc_dds_SRP014652)
############## DEG analysis ############################
dds_SRP014652 <- DESeq(dds_SRP014652)
dds_SRP014652_results <- results(dds_SRP014652, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP014652_results)
dn_SRP014652 <- rownames(subset(dds_SRP014652_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP014652)
