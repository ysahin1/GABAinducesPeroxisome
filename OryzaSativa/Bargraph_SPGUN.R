library(DESeq2)
library(rafalib)
library(ggplot2)
library(plyr)
library(ggsci)
library(RColorBrewer)
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
fr_df <- data.frame()
for(i in 1:ncol(zscore)){
  ex = 0
  cn <- colnames(zscore)[i]
  for (ii in rownames(zscore)) {
    ge <- zscore[ii,i]
    if (ge > ex) {
      gi <- match(ii,rice_per$os_locus)
      cl <- rice_per$functions[gi]
      fr_df[ii,cn] <- cl
    }
    else (fr_df[ii,cn] <- "DOWN REGULATED")
  }
}
gdf <- data.frame()
for (i in 1:ncol(fr_df)){
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
            position = position_stack(vjust = 0),
            color="black", 
            size=3.5)+
  scale_fill_manual(values = mycolors)+
  theme(legend.position = "panel",
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",
                                 size = 0.5, 
                                 linetype = "solid"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
jpeg(filename = "SPGUN.jpeg", res = "600", width = "5", height = "5", units = "in")
p1
dev.off()
