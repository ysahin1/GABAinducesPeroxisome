library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
SRP014652_raw_counts <- read.table("../counts/SRP014652_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(SRP014652_raw_counts)
row.names(SRP014652_raw_counts) <- SRP014652_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_SRP014652 <- SRP014652_raw_counts$Length
names(gene_length_SRP014652) <- SRP014652_raw_counts$Geneid
##################################################
SRP014652_raw_counts <- SRP014652_raw_counts[,c(7,8,9,10)] #c(1:7),10) exclude all information that do not include counts
names(SRP014652_raw_counts) <- c(#"SRR1619633",
                                 #"SRR1619634",
                                 #"SRR1619636",
                                 #"SRR1619637",
                                 #"SRR1619660",
                                 #"SRR1619661",
                                 #"SRR1619664",
                                 #"SRR1619667",
                                 #"SRR1619668",
                                 #"SRR1619669",
                                 #"SRR1619670",
                                 #"SRR1626426"
                                 "SRR1619633",
                                 "SRR1619634",
                                 "SRR1619636",
                                 "SRR1619637")
SRP014652_grow_CG <- c(#"6_DAP","6_DAP",
                       #"12_DAP","12_DAP",
                       #"18_DAP","18_DAP",
                       #"24_DAP",
                       #"30_DAP","30_DAP",
                       #"0_DAP","0_DAP"
                       #"24_DAP"
                       "6_DAP","6_DAP",
                       "12_DAP","12_DAP")
#SRP014652_grow_CG <- factor(SRP014652_grow_CG, levels = c("control","drought"))
#SRP014652_grow_CG <- relevel(SRP014652_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(SRP014652_raw_counts)
SRP014652_sample_info <- data.frame(condition=SRP014652_grow_CG,
                                    row.names=names(SRP014652_raw_counts))

dds_SRP014652 <- DESeqDataSetFromMatrix(countData = SRP014652_raw_counts,
                                        colData = SRP014652_sample_info,
                                        design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_SRP014652))

##################################
### filter out count datas and normalization
################################## 
dds_SRP014652 <- dds_SRP014652[rowSums(counts(dds_SRP014652))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_SRP014652)) # number of genes counted
head(assay(dds_SRP014652))

######normalization
dds_SRP014652 <- estimateSizeFactors(dds_SRP014652)
sizeFactors(dds_SRP014652)
sfn_dds_SRP014652 <- counts(dds_SRP014652, normalized = TRUE)

########transformation
ln_dds_SRP014652 <- log2(sfn_dds_SRP014652+1)


### visulaise norm and log norm data#####
par(mfrow = c (1 ,2), cex.axis=0.5)
boxplot(sfn_dds_SRP014652, notch = TRUE,
        main = "URCO6V12", 
        ylab = " ",
        las = 2)

boxplot(ln_dds_SRP014652 , notch = TRUE ,
        main = "Log2-TSFNGAOR6V12" ,
        ylab = " ",
        las=2)
################# check if your data homoskedastic
plot (ln_dds_SRP014652 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_SRP014652 <- vst(dds_SRP014652, blind = TRUE)
vstc_dds_SRP014652 <- assay(vst_dds_SRP014652)



######visualize vst norm counts#####
msd_plot_1 <- meanSdPlot (ln_dds_SRP014652,
                          ranks = FALSE, # show the data on the original scale
                          plot = FALSE)
msd_plot_1$gg +
  ggtitle("SDNlog2-TRCIMVE") +
  ylab("SD")

msd_plot_2 <- meanSdPlot(vstc_dds_SRP014652,
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

dis_vstc_dds_SRP014652 <- as.dist(1- cor(vstc_dds_SRP014652, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_SRP014652),
     labels = colnames(dis_vstc_dds_SRP014652),
     xlab = "",
     main = "Hierarchical Clustering of SRP014652")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_SRP014652))
plot(pc$x[ ,1],pc$x[ ,2],
     col = colData(dds_SRP014652)[ ,1],
     main = "Principal Component Analysis")


P <- plotPCA(vst_dds_SRP014652)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of SRP014652") + geom_text(aes(label=names(SRP014652_raw_counts)))
print (P)

############## DEG analysis ############################
dds_SRP014652 <- DESeq(dds_SRP014652)
dds_SRP014652 <- estimateSizeFactors(dds_SRP014652)
dds_SRP014652 <- estimateDispersions(dds_SRP014652)
dds_SRP014652 <- nbinomWaldTest(dds_SRP014652)

dds_SRP014652_results <- results(dds_SRP014652, 
                                 independentFiltering = TRUE, 
                                 alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP014652_results)

dn_SRP014652 <- rownames(subset(dds_SRP014652_results, 
                                padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP014652)
write.csv(dds_SRP014652_results[dn_SRP014652,],"12DPA_vs_6DPA_degs.csv")
resultsNames(dds_SRP014652)
