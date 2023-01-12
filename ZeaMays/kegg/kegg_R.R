library(goseq)
ko_slim_table <- read.delim2("/media/yunus/TOSHIBA1TB/reference_genomes/maizeB73/annotations/maize/download.20210813.233442/phytozome/Zmays/RefGen_V4/annotation/Zmays_493_RefGen_V4.annotation_info.txt", header=FALSE, comment.char="#")
head(ko_slim_table)
zm_ko<- data.frame(gene_id=ko_slim_table$V2,
                   ko=ko_slim_table$V9)
zm_ko <- zm_ko[!(is.na(zm_ko$ko) | zm_ko$ko==""),]
zm_ko_1 <- strsplit(as.vector(zm_ko$ko), "," )
nIDs <- sapply(1:length(zm_ko_1), function(x) length(zm_ko_1[[x]]))
zm_ko <- data.frame(gene_symbol = rep(zm_ko$gene_id, nIDs), 
                    ko = unlist(zm_ko_1), 
                    stringsAsFactors=FALSE )
zm_ko <- zm_ko[!duplicated(zm_ko$gene_symbol), ]
gene_length_SRP014652
length(zm_ko$gene_symbol) # annotated gene number
factor_table_all_genes <- factor(as.integer(zm_ko$gene_symbol %in% dn_SRP014652))
length(factor_table_all_genes)
names(factor_table_all_genes) = zm_ko$gene_symbol # should be named vector
###########GOSEQ analizi####################################
length(gene_length_SRP014652)
length(gene_length_SRP014652[names(factor_table_all_genes)])

pwf <- nullp(factor_table_all_genes, 
             bias.data = gene_length_SRP014652[names(factor_table_all_genes)])
KO.wall <- goseq(pwf,gene2cat=zm_ko) #method = "Hypergeometric",use_genes_without_cat=TRUE

en_KO <- KO.wall$category[KO.wall$over_represented_pvalue<.05]
length(en_KO)
zm_ko[ which(zm_ko$ko %in% en_KO),]
write.table(zm_ko[ which(zm_ko$ko %in% en_KO),],file = "ko_zm_R.txt", quote = F,sep = "\t", row.names = F)
getGeneLists<- function(pwf, en_KO){
  out <- list()
  for(term in en_KO){
    tmp <- pwf[zm_ko[which(zm_ko$ko %in% term),1],]
    tmp <- rownames(tmp[tmp$DEgenes %in% 1,])
    out[[term]] <- as.vector(tmp)
  }
  out
}
goList <- getGeneLists(pwf, en_KO)
KO.wall$Genes <- sapply(KO.wall$category, function(x) goList[[x]])

entry.kegg.path <- list()
for(i in 1:length(kegg.ids$kegg_id))
{
  x <- kegg.ids$kegg_id[i]
  query <- tryCatch(KEGGREST::keggLink("pathway", kegg.ids$kegg_id[i]), error=function(e) NULL)
  entry.kegg.path[[x]] <- as.vector(query)
}





geneID2GO <- by(GOs_GSE57950h$go,
                GOs_GSE57950h$gene_id,
                function(x) as.character(x))