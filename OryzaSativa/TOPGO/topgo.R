library(topGO)
GO_slim_table <- read.delim2("/media/yunus/TOSHIBA1TB/reference_genomes/oryza_sativa/MSU7/GO_annotation.txt", header=FALSE, comment.char="#")
GOs_GSE155932<- data.frame(gene_id=GO_slim_table$V2,
                           go=GO_slim_table$V3)
geneID2GO <- by(GOs_GSE155932$go,
                GOs_GSE155932$gene_id,
                function(x) as.character(x))
# create named factor vector
all.genes <- sort(unique(as.character(GOs_GSE155932$gene_id)))
int.genes <-  rownames(subset(dds_GSE155932_results, padj < 0.05))# sig. genes
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes
#####3 create topGo object

############ BP ######################
go_BP_GSE155932<- new("topGOdata"
                      , ontology='BP'
                      , allGenes = int.genes
                      , nodeSize = 1
                      , annot = annFUN.gene2GO
                      , gene2GO = geneID2GO)

### make your analysis
resultFisher_BP_GSE155932<- runTest(go_BP_GSE155932, 
                                    algorithm = "classic", 
                                    statistic = "fisher")
number_bp_top_nodes <- length(score(resultFisher_BP_GSE155932))
result_BP_GSE155932<- GenTable(go_BP_GSE155932, 
                               classicFisher = resultFisher_BP_GSE155932, 
                               orderBy = "classicFisher", 
                               ranksOf = "classicFisher", 
                               topNodes = number_bp_top_nodes)
length(score(resultFisher_BP_GSE155932))
result_BP_GSE155932$adj.pvalue <- p.adjust(score(resultFisher_BP_GSE155932),
                                           method="BH")
result_BP_GSE155932$ontology <- rep("BP", 
                                    dim(result_BP_GSE155932)[1])
result_BP_GSE155932$ratio <- result_BP_GSE155932$Significant/result_BP_GSE155932$Annotated

# select those terms with a p-value < 0.05
result_BP_GSE155932$genes <- sapply(result_BP_GSE155932$GO.ID, function(x)
{
  genes<-genesInTerm(go_BP_GSE155932, x)
  genes[[1]][genes[[1]] %in% dn_GSE155932] # myGenes is the queried gene list
})
head(result_BP_GSE155932)
enriched_BP_GSE155932<- result_BP_GSE155932[result_BP_GSE155932$adj.pvalue < .05,]
#Number of Enriched BP_GO terms
dim(enriched_BP_GSE155932)[1]
write.csv(enriched_BP_GSE155932,"rice_S5_vs_S4.csv")
# Bu çalışmada GO enrichment açısından önemli olan gen sayısı
table(dn_GSE155932%in% unique(unlist(result_BP_GSE155932$genes)))[2]
################# MF #####################
enriched_TOPGOs_GSE155932<- rbind(enriched_BP_GSE155932)
dim(enriched_TOPGOs_GSE155932)
####################visualize results###################
########################################################
names(enriched_TOPGOs_GSE155932) <- c("category",
                                      "term",
                                      "numInCat",
                                      "Gene_Number",
                                      "Expected",
                                      "over_represented_pvalue",
                                      "Q_value",
                                      "ontology",
                                      "Rich_Factor",
                                      "Genes")
dottplot <- function(df, showCategory){
  df <- df[with(df, order(Rich_Factor, Q_value, decreasing = c(TRUE, FALSE))),]
  df <- head(df, n=showCategory)
  d_plot <- ggplot(df, aes_string(x="term", 
                                  y="Rich_Factor", 
                                  colour="Q_value",
                                  size="Gene_Number")) + 
    geom_point() +
    facet_grid(~ontology)+
    scale_color_gradient(low="#FF0000",
                         high="#000000") +
    coord_flip() +
    theme_bw(base_size=9)
  return(d_plot)
} 
library(ggplot2)
dottplot(enriched_TOPGOs_GSE155932, 
         showCategory = 107)



library(GOplot)
DEGs.table <- dds_GSE155932_results[dn_GSE155932,]
genes=DEGs.table
terms=enriched_TOPGOs_GSE155932
names(terms) <- c("Category",
                  "term",
                  "numInCat",
                  "numDEInCat",
                  "under_represented_pvalue",
                  "over_represented_pvalue",
                  "adj_pval",
                  "ontology",
                  "Rich_Factor",
                  "Genes"
                  
)
names(genes) <- c("baseMean",
                  "logFC",
                  "lfcSE",
                  "stat",
                  "P.Value",
                  "adj.P.Val")
genes$ID <- rownames(genes)
logFC <- sapply(as.vector(unlist(terms$Genes)), function(x) genes$logFC[match(x, genes$ID)])
count <- sapply(1:length(names(terms$Genes)), function(x) length(terms$Genes[[x]]))
#logFC[is.na(logFC)] <- 0 #Fold change değeri olmayanlara sıfır değeri verdim 

s <- 1; zsc <- c()
for (c in 1:length(count)){
  value <- 0
  e <- s + count[c] - 1
  value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
  zsc <- c(zsc, sum(value)/sqrt(count[c]))
  s <- e + 1
}

df_for_go_circle <- data.frame(category = rep(as.character(terms$Category), count), 
                               term = rep(as.character(terms$term), count),
                               genes = as.vector(unlist(terms$Genes)), 
                               logFC = as.vector(logFC), 
                               adj_pval = rep(terms$adj_pval, count),
                               zscore = rep(zsc, count), 
                               stringsAsFactors = FALSE)


enriched_BP_GSE155932[order(enriched_BP_GSE155932$adj.pvalue, decreasing = FALSE)[1:3],-c(10)]

unique(df_for_go_circle$term[order(df_for_go_circle$term)])
enriched_term <- "cellular process"
round(min(df_for_go_circle$logFC[df_for_go_circle$term == enriched_term]), digits = 2)
round(max(df_for_go_circle$logFC[df_for_go_circle$term == enriched_term]), digits = 2)


