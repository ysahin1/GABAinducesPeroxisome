library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
SRP014652v_raw_counts <- read.table("../counts/SRP014652_v_stage_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(SRP014652v_raw_counts)
row.names(SRP014652v_raw_counts) <- SRP014652v_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_SRP014652v <- SRP014652v_raw_counts$Length
names(gene_length_SRP014652v) <- SRP014652v_raw_counts$Geneid
##################################################
SRP014652v_raw_counts <- SRP014652v_raw_counts[,-c(1:18)] #c(1:7),10) exclude all information that do not include counts
names(SRP014652v_raw_counts) <- c(#"SRR531218",
                        #"SRR531219",
                        #"SRR531220",
                        #"SRR531866",
                        #"SRR531867",
                        #"SRR531868",
                        #"SRR531869",
                        #"SRR531870",
                        #"SRR531871",
                        #"SRR531872",
                        #"SRR531873",
                        #"SRR531874",
                        "SRR940237",
                        "SRR940238",
                        "SRR940249",
                        "SRR940250",
                        "SRR940252",
                        "SRR940253",
                        "SRR940254",
                        "SRR940255",
                        "SRR940256",
                        "SRR940257",
                        "SRR940260",
                        "SRR940274",
                        "SRR940276",
                        "SRR940278"
                        )
SRP014652v_grow_CG <- c(#"topmost_leaf_v3",
                        #"topmost_leaf_v3",
                        #"topmost_leaf_v3",
                        #"base_stage2_leaf_v5",
                        #"base_stage2_leaf_v5",
                        #"base_stage2_leaf_v5",
                        #"tip_stage2_leaf_v7",
                        #"tip_stage2_leaf_v7",
                        #"tip_stage2_leaf_v7",
                        #"base_stage2_leaf_v7",
                        #"base_stage2_leaf_v7",
                        #"base_stage2_leaf_v7",
                        "eighth_leaf_v9",
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
#SRP014652v_grow_CG <- factor(SRP014652v_grow_CG, levels = c("control","drought"))
#SRP014652v_grow_CG <- relevel(SRP014652v_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(SRP014652v_raw_counts)
SRP014652v_sample_info <- data.frame(condition=SRP014652v_grow_CG,
                          row.names=names(SRP014652v_raw_counts))

dds_SRP014652v <- DESeqDataSetFromMatrix(countData = SRP014652v_raw_counts,
                                   colData = SRP014652v_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_SRP014652v))

##################################
### filter out count datas and normalization
################################## 
dds_SRP014652v <- dds_SRP014652v[rowSums(counts(dds_SRP014652v))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_SRP014652v)) # number of genes counted
head(assay(dds_SRP014652v))

######normalization
dds_SRP014652v <- estimateSizeFactors(dds_SRP014652v)
sizeFactors(dds_SRP014652v)
sfn_dds_SRP014652v <- counts(dds_SRP014652v, normalized = TRUE)

########transformation
ln_dds_SRP014652v <- log2(sfn_dds_SRP014652v+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_SRP014652v, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_SRP014652v , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_SRP014652v [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_SRP014652v <- vst(dds_SRP014652v, blind = TRUE)
vstc_dds_SRP014652v <- assay(vst_dds_SRP014652v)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_SRP014652v,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_SRP014652v,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_SRP014652v <- as.dist(1- cor(vstc_dds_SRP014652v, method ="pearson"))

#######visiulaze HC ####
library(ggplot2)
library(ggdendro)
p_h <- plot(hclust(dis_vstc_dds_SRP014652v),
      labels = colnames(dis_vstc_dds_SRP014652v),
      xlab = "",
      main = "Hierarchical Clustering of SRP014652")

model <- hclust(dis_vstc_dds_SRP014652v, "ave")
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
P <- plotPCA(vst_dds_SRP014652v)

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
############## DEG analysis ############################
dds_SRP014652v <- DESeq(dds_SRP014652v)
dds_SRP014652v <- estimateSizeFactors(dds_SRP014652v)
dds_SRP014652v <- estimateDispersions(dds_SRP014652v)
dds_SRP014652v <- nbinomWaldTest(dds_SRP014652v)

dds_SRP014652v_results <- results(dds_SRP014652v, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP014652v_results)

dn_SRP014652v <- rownames(subset(dds_SRP014652v_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP014652v)

resultsNames(dds_SRP014652v)
