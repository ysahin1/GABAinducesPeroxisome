library(goseq)
ko_slim_table <- read.delim2("/media/yunus/TOSHIBA1TB/reference_genomes/oryza_sativa/MSU7/annotation/Osativa_323_v7.0.annotation_info.txt", header=T, comment.char="#")
head(ko_slim_table)
os_ko<- data.frame(gene_id=ko_slim_table$locusName,
                           ko=ko_slim_table$KO)
os_ko <- os_ko[!(is.na(os_ko$ko) | os_ko$ko==""),]
os_ko_1 <- strsplit(as.vector(os_ko$ko), "," )
nIDs <- sapply(1:length(os_ko_1), function(x) length(os_ko_1[[x]]))
os_ko <- data.frame(gene_symbol = rep(os_ko$gene_id, nIDs), 
                                 ko = unlist(os_ko_1), 
                                 stringsAsFactors=FALSE )
os_ko <- os_ko[!duplicated(os_ko$gene_symbol), ]5
gene_length_GSE155932
length(os_ko$gene_symbol) # annotated gene number
factor_table_all_genes <- factor(as.integer(os_ko$gene_symbol %in% dn_GSE155932))
length(factor_table_all_genes)
names(factor_table_all_genes) = os_ko$gene_symbol # should be named vector
###########GOSEQ analizi####################################
length(gene_length_GSE155932)
length(gene_length_GSE155932[names(factor_table_all_genes)])

pwf <- nullp(factor_table_all_genes, 
             bias.data = gene_length_GSE155932[names(factor_table_all_genes)])
KO.wall <- goseq(pwf,gene2cat=os_ko) #method = "Hypergeometric",use_genes_without_cat=TRUE

en_KO <- KO.wall$category[KO.wall$over_represented_pvalue<.05]
length(en_KO)
os_ko[ which(os_ko$ko %in% en_KO),]
write.table(os_ko[ which(os_ko$ko %in% en_KO),],file = "ko_os.txt", quote = F,sep = "\t", row.names = F)
  getGeneLists<- function(pwf, en_KO){
  out <- list()
  for(term in en_KO){
    tmp <- pwf[os_ko[which(os_ko$ko %in% term),1],]
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