library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE155932_raw_counts <- read.table("../counts/GSE155932_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE155932_raw_counts)
row.names(GSE155932_raw_counts) <- GSE155932_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE155932 <- GSE155932_raw_counts$Length
names(gene_length_GSE155932) <- GSE155932_raw_counts$Geneid
##################################################
GSE155932_raw_counts <- GSE155932_raw_counts[,c(8,9,10,11)] #c(1:7),10) exclude all information that do not include counts
names(GSE155932_raw_counts) <- c(
                        #"GSM4716039",
                        #"GSM4716040",
                        #"GSM4716041",
                        #"GSM4716042",
                        #"GSM4716043",
                        #"GSM4716044",
                        #"GSM4716035",
                        #"GSM4716036",
                        #"GSM4716037",
                        #"GSM4716038",
                        "GSM4716039",
                        "GSM4716040",
                        "GSM4716041",
                        "GSM4716042")
GSE155932_grow_CG <- c(
                       #"s2","s2",
                       #"s4","s4",
                       #"s5","s5",
                       #"s1","s1",
                       "s3","s3",
                       "s4","s4")
#GSE155932_grow_CG <- factor(GSE155932_grow_CG, levels = c("control","drought"))
#GSE155932_grow_CG <- relevel(GSE155932_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE155932_raw_counts)
GSE155932_sample_info <- data.frame(condition=GSE155932_grow_CG,
                          row.names=names(GSE155932_raw_counts))

dds_GSE155932 <- DESeqDataSetFromMatrix(countData = GSE155932_raw_counts,
                                   colData = GSE155932_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE155932))

##################################
### filter out count datas and normalization
################################## 
dds_GSE155932 <- dds_GSE155932[rowSums(counts(dds_GSE155932))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE155932)) # number of genes counted
head(assay(dds_GSE155932))

######normalization
dds_GSE155932 <- estimateSizeFactors(dds_GSE155932)
sizeFactors(dds_GSE155932)
sfn_dds_GSE155932 <- counts(dds_GSE155932, normalized = TRUE)

########transformation
ln_dds_GSE155932 <- log2(sfn_dds_GSE155932+1)


### visulaise norm and log norm data#####
par(mfrow = c (1 ,2), cex.axis=0.5)
boxplot(sfn_dds_GSE155932, notch = TRUE,
        main = "URCOS4VS3", 
        ylab = " ",
        las = 2)

boxplot(ln_dds_GSE155932 , notch = TRUE ,
        main = "Log2-TSFNGAORS4VS3" ,
        ylab = " ",
        las=2)
################# check if your data homoskedastic
plot (ln_dds_GSE155932 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE155932 <- vst(dds_GSE155932, blind = TRUE)
vstc_dds_GSE155932 <- assay(vst_dds_GSE155932)



######visualize vst norm counts#####
msd_plot_1 <- meanSdPlot (ln_dds_GSE155932,
                          ranks = FALSE, # show the data on the original scale
                          plot = FALSE)
msd_plot_1$gg +
  ggtitle("SDNlog2-TRCIMVE") +
  ylab("SD")

msd_plot_2 <- meanSdPlot(vstc_dds_GSE155932,
                         ranks = FALSE, # show the data on the original scale
                         plot = FALSE)
msd_plot_2$gg +
  ggtitle ( "EODTIMVE") +
  ylab ( "SD")
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}
multiplot(msd_plot_1, msd_plot_2, cols=2)

######## hierarchical clustering##########

dis_vstc_dds_GSE155932 <- as.dist(1- cor(vstc_dds_GSE155932, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE155932),
      labels = colnames(dis_vstc_dds_GSE155932),
      xlab = "",
      main = "Hierarchical Clustering of GSE155932")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE155932))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE155932)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE155932)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE155932") + geom_text(aes(label=names(GSE155932_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE155932 <- DESeq(dds_GSE155932)
dds_GSE155932 <- estimateSizeFactors(dds_GSE155932)
dds_GSE155932 <- estimateDispersions(dds_GSE155932)
dds_GSE155932 <- nbinomWaldTest(dds_GSE155932)

dds_GSE155932_results <- results(dds_GSE155932, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE155932_results)

dn_GSE155932 <- rownames(subset(dds_GSE155932_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE155932)
write.csv(dds_GSE155932_results[dn_GSE155932,],"S4_vs_S3_degs.csv")
resultsNames(dds_GSE155932)
write.table(as.character(dds_GSE155932_results[dn_GSE155932,2]), 
            "os_fc.txt", 
            quote = F, 
            sep = "\t",
            row.names = as.vector( paste(dn_GSE155932,1, sep = ".")))
