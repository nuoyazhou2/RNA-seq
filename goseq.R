
setwd("/Users/Mingyang/Google Drive/Tong_RIP/goseq")
library(goseq)
library(KEGG.db)

cutoff_padj = 0.1
cutoff_logFC = log2(1.5)
supportedOrganisms()[supportedOrganisms()$Genome=="mm9",]

#file = "DESeq2_CTCF_3F_vs_input_3F.txt"

myfun = function(x) {
  symbolgenes_for_each_go = NULL
  for (i in seq(length(x))) {
    symbolgenes_for_each_go = c(symbolgenes_for_each_go, someenv[[x[i]]])
  }
  return (symbolgenes_for_each_go)
}

for (file in list.files("../DESeq2/", "^DESeq2.*txt$")) {
  table = read.table(paste0("../DESeq2/",file))
  factor = sub("DESeq2_(.*)_(3F|E14)_.*.txt", "\\1", file)
  cell = sub("DESeq2_(.*)_(3F|E14)_.*.txt", "\\2", file)
  
  # universe genes are all genes
  genes = as.integer(table$padj<=cutoff_padj & table$log2FoldChange>=cutoff_logFC)
  genes[is.na(genes)] = 0
  names(genes) = row.names(table)
  table(genes,useNA = "ifany")
  
  # universe genes are those with non-NAs padj
  # genes = as.integer(table[!is.na(table$padj),]$padj<=cutoff_padj)
  # names(genes) = row.names(table[!is.na(table$padj),])
  # table(genes,useNA = "ifany")
  
  pwf=nullp(genes,"mm9","geneSymbol")
  head(pwf)
  
  GO.wall=goseq(pwf,"mm9","geneSymbol")
  
  # GO.samp=goseq(pwf,"mm9","geneSymbol",method="Sampling",repcnt=1000)
  # head(GO.samp)
  # 
  # 
  # plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]), 
  #      xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)", 
  #      xlim=c(-3,0))
  # abline(0,1,col=3,lty=2)
  # 
  # 
  # GO.nobias=goseq(pwf,"mm9","geneSymbol",method="Hypergeometric")
  # head(GO.nobias)
  # 
  # plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.wall[,1]),2]), 
  #      xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)", 
  #      xlim=c(-3,0), ylim=c(-3,0))
  # abline(0,1,col=3,lty=2)
  
  # GO.MF=goseq(pwf,"mm9","geneSymbol",test.cats=c("GO:MF"))
  # head(GO.MF)
  
  GO.wall$padj = p.adjust(GO.wall$over_represented_pvalue,method="BH")
  # enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<1]
  # head(enriched.GO)
  
  # library(GO.db)
  # for(go in enriched.GO[1:10]){
  #   print(GOTERM[[go]])
  #   cat("--------------------------------------\n")
  # }
  genes2go=getgo(rownames(table[which(table$padj<=cutoff_padj & table$log2FoldChange>=cutoff_logFC),]),'mm9','geneSymbol')
  go2genes=goseq:::reversemapping(genes2go)
  GO.wall$DEgenes = sapply(go2genes[match(GO.wall$category,names(go2genes))], paste, collapse = ",")
  head(GO.wall)
  
  
  #pwf=nullp(genes,'mm9','geneSymbol')
  KEGG=goseq(pwf,'mm9','geneSymbol',test.cats="KEGG")
  KEGG$padj = p.adjust(KEGG$over_represented_pvalue,method="BH")
  
  pathway = stack(mget(KEGG$category, KEGGPATHID2NAME)) # Get pathway names patways
  KEGG$pathway = pathway$values
  
  KEGG2ENTREZ = as.list(org.Mm.egPATH2EG)
  de_genes = select(org.Mm.eg.db, keys=as.vector(rownames(table[which(table$padj<=cutoff_padj & table$log2FoldChange>=cutoff_logFC),])),
                          keytype="SYMBOL", columns=c("ENTREZID"))
  de_genes = de_genes[which(!is.na(de_genes$ENTREZID)),]
  
  someenv = new.env()
  for(i in seq(nrow(de_genes)))
  {
    someenv[[ de_genes[i,2] ]]<- de_genes[i,1]
  }

  temp = lapply(KEGG2ENTREZ[match(KEGG$category,names(KEGG2ENTREZ))],myfun)
  KEGG$DEgenes = sapply(temp, paste, collapse = ",")
  head(KEGG)
  
  write.table(as.data.frame(GO.wall),file=paste0("GO_",factor,"_",cell,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(as.data.frame(GO.wall[which(GO.wall$ontology == "CC"),]),file=paste0("GO_CC_",factor,"_",cell,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(as.data.frame(GO.wall[which(GO.wall$ontology == "MF"),]),file=paste0("GO_MF_",factor,"_",cell,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(as.data.frame(GO.wall[which(GO.wall$ontology == "BP"),]),file=paste0("GO_BP_",factor,"_",cell,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(as.data.frame(KEGG),file=paste0("KEGG_",factor,"_",cell,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
}

sessionInfo()


