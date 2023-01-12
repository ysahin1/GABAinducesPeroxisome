library(DESeq2)
library(rafalib)
library(ggplot2)
library(plyr)
library(ggsci)
library(RColorBrewer)
SRP014652v_raw_counts <- read.table("../counts/SRP014652_v_stage_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
row.names(SRP014652v_raw_counts) <- SRP014652v_raw_counts$Geneid
SRP014652v_raw_counts <- SRP014652v_raw_counts[,-c(1:18)] #c(1:7),10) exclude all information that do not include counts
names(SRP014652v_raw_counts) <- c(
  "8th_V9_R1",
  "8th_V9_R2",
  "11th_V9_R1",
  "11th_V9_R2",
  "13th_V9_R1",
  "13th_V9_R2",
  "13th_V9_R3",
  "im_V9_R1",
  "im_V9_R2",
  "im_V9_R3",
  "13th_VT_R1",
  "13th_VT_R2",
  "13th_R2_R1",
  "13th_R2_R2"
)
SRP014652v_grow_CG <- c("eighth_leaf_v9",
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
                        "thirteenth_leaf_r2")
SRP014652v_grow_CG <- factor(SRP014652v_grow_CG, levels = c("immature_Leaves_v9",
                                                            "eighth_leaf_v9",
                                                            "eleventh_leaf_v9",
                                                            "thirteenth_leaf_v9",
                                                            "thirteenth_leaf_vT",
                                                            "thirteenth_leaf_r2"))

SRP014652v_sample_info <- data.frame(condition=SRP014652v_grow_CG,
                                     row.names=names(SRP014652v_raw_counts))
dds_SRP014652v <- DESeqDataSetFromMatrix(countData = SRP014652v_raw_counts,
                                         colData = SRP014652v_sample_info,
                                         design = ~ condition)


dds_SRP014652v <- dds_SRP014652v[rowSums(counts(dds_SRP014652v))> 100,]#eleminate counts that are less than zero
vst_dds_SRP014652v <- vst(dds_SRP014652v, blind = TRUE)
vstc_dds_SRP014652v <- assay(vst_dds_SRP014652v)
zm_per <- read.delim("/media/yunus/TOSHIBA1TB/meta-seq_leaf_development/DEG/maize/peroxisome_homologs/zm_os_peroxisome_homolog_table.csv", stringsAsFactors = F)
per_genes <-  zm_per$zm_locus
indx <- match(per_genes,rownames(vstc_dds_SRP014652v))
indx <- indx[!is.na(indx)]
perg_vst <- vstc_dds_SRP014652v[indx,]
indx_symbol <- match(rownames(perg_vst),zm_per$zm_locus)
per_symbol <- as.vector(zm_per$symbols[indx_symbol])
zscore <- (perg_vst-rowMeans(perg_vst))/(rowSds(as.matrix(perg_vst)))
fr_df <- data.frame()
for(i in 1:ncol(zscore)){
  ex = 0
  cn <- colnames(zscore)[i]
  for (ii in rownames(zscore)) {
    ge <- zscore[ii,i]
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
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",
                                 size = 0.5, 
                                 linetype = "solid"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

p1
library(ggpp)
leg <- get_legend(p1)
as_ggplot(leg)

scale_x_discrete(limits=c("0_DAP_R1",
                          "0_DAP_R2",
                          "6_DAP_R1",
                          "6_DAP_R2"))