library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE151455_raw_counts <- read.table("../counts/GSE151455_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE151455_raw_counts)
row.names(GSE151455_raw_counts) <- GSE151455_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE151455 <- GSE151455_raw_counts$Length
names(gene_length_GSE151455) <- GSE151455_raw_counts$Geneid
##################################################
GSE151455_raw_counts <- GSE151455_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE151455_raw_counts) <- c("GSM4578835",
                        "GSM4578836",
                        "GSM4578837",
                        "GSM4578844",
                        "GSM4578845",
                        "GSM4578846",
                        "GSM4578847",
                        "GSM4578848",
                        "GSM4578849")
GSE151455_grow_CG <- c("65_DAP","65_DAP","65_DAP",
                       "10_DAP","10_DAP","10_DAP","10_DAP","10_DAP","10_DAP")
#GSE151455_grow_CG <- factor(GSE151455_grow_CG, levels = c("control","drought"))
#GSE151455_grow_CG <- relevel(GSE151455_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE151455_raw_counts)
GSE151455_sample_info <- data.frame(condition=GSE151455_grow_CG,
                          row.names=names(GSE151455_raw_counts))

dds_GSE151455 <- DESeqDataSetFromMatrix(countData = GSE151455_raw_counts,
                                   colData = GSE151455_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE151455))

##################################
### filter out count datas and normalization
################################## 
dds_GSE151455 <- dds_GSE151455[rowSums(counts(dds_GSE151455))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE151455)) # number of genes counted
head(assay(dds_GSE151455))

######normalization
dds_GSE151455 <- estimateSizeFactors(dds_GSE151455)
sizeFactors(dds_GSE151455)
sfn_dds_GSE151455 <- counts(dds_GSE151455, normalized = TRUE)

########transformation
ln_dds_GSE151455 <- log2(sfn_dds_GSE151455+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE151455, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE151455 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE151455 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE151455 <- vst(dds_GSE151455, blind = TRUE)
vstc_dds_GSE151455 <- assay(vst_dds_GSE151455)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE151455,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE151455,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE151455 <- as.dist(1- cor(vstc_dds_GSE151455, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE151455),
      labels = colnames(dis_vstc_dds_GSE151455),
      xlab = "",
      main = "Hierarchical Clustering of GSE151455")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE151455))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE151455)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE151455)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE151455") + geom_text(aes(label=names(GSE151455_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE151455 <- DESeq(dds_GSE151455)
dds_GSE151455 <- estimateSizeFactors(dds_GSE151455)
dds_GSE151455 <- estimateDispersions(dds_GSE151455)
dds_GSE151455 <- nbinomWaldTest(dds_GSE151455)

dds_GSE151455_results <- results(dds_GSE151455, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE151455_results)

dn_GSE151455 <- rownames(subset(dds_GSE151455_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE151455)

resultsNames(dds_GSE151455)
