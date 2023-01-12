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
GSE155932_raw_counts <- GSE155932_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE155932_raw_counts) <- c("GSM4716038",
                        "GSM4716039",
                        "GSM4716040",
                        "GSM4716041",
                        "GSM4716042",
                        "GSM4716043",
                        "GSM4716044",
                        "GSM4716035",
                        "GSM4716036",
                        "GSM4716037")
GSE155932_grow_CG <- c("s2",
                       "s3","s3",
                       "s4","s4",
                       "s5","s5",
                       "s1","s1",
                       "s2")
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
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE155932, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE155932 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE155932 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE155932 <- vst(dds_GSE155932, blind = TRUE)
vstc_dds_GSE155932 <- assay(vst_dds_GSE155932)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE155932,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE155932,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE155932 <- as.dist(1- cor(vstc_dds_GSE155932, method ="pearson"))

#######visiulaze HC ####
library(ggplot2)
library(ggdendro)
model <- hclust(dis_vstc_dds_GSE155932, "ave")
dhc <- as.dendrogram(model)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")
p_hc <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = ddata$labels, 
            aes(x = x, y = y, label = label), size = 3, hjust= 0,vjust = 0, check_overlap = TRUE) +
  coord_flip() + 
  scale_y_reverse(expand=c(0.3,0))
p_hc
P <- plotPCA(vst_dds_GSE155932)

P <- P + theme_bw()
print (P)
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
multiplot(p_hc,P, cols=2)
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

resultsNames(dds_GSE155932)
write.table(as.character(dds_GSE155932_results[dn_GSE155932,2]), 
            "os_fc.txt", 
            quote = F, 
            sep = "\t",
            row.names = as.vector( paste(dn_GSE155932,1, sep = ".")))
