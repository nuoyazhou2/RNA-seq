#!/home/rcf-proj/wl/caim/usr/R/bin/R
#for installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("ReportingTools")
library('DESeq2')
library("ggplot2")
library("RColorBrewer")
library("pheatmap")

directory <- "../htseq"
sampleFiles <- paste0(rep(paste0("htseq_sample",c(1:14),"_"), each=3),c(1:3),".txt")
sampleCondition <- c(rep("CTCF_E14",3), rep("SMC1_E14",3), rep("MED12_E14",3), rep("KLF4_E14",3), rep("FLAG_E14",3), rep("EZH2_E14",3), rep("input_E14",3),rep("CTCF_3F",3), rep("SMC1_3F",3), rep("MED12_3F",3), rep("KLF4_3F",3), rep("FLAG_3F",3), rep("EZH2_3F",3), rep("input_3F",3))
sampleName <- c(paste0(rep("CTCF_E14_R",3),c(1:3)), paste0(rep("SMC1_E14_R",3),c(1:3)), paste0(rep("MED12_E14_R",3),c(1:3)), paste0(rep("KLF4_E14_R",3),c(1:3)), paste0(rep("FLAG_E14_R",3),c(1:3)), paste0(rep("EZH2_E14_R",3),c(1:3)), paste0(rep("input_E14_R",3),c(1:3)), paste0(rep("CTCF_3F_R",3),c(1:3)), paste0(rep("SMC1_3F_R",3),c(1:3)), paste0(rep("MED12_3F_R",3),c(1:3)), paste0(rep("KLF4_3F_R",3),c(1:3)), paste0(rep("FLAG_3F_R",3),c(1:3)), paste0(rep("EZH2_3F_R",3),c(1:3)), paste0(rep("input_3F_R",3),c(1:3)))
sampleFactor <- sub("(.*)_E14|_3F","\\1",sampleCondition)
sampleCell <- c(rep("E14",21),rep("3F",21))
sampleTable <- data.frame(sampleName = sampleName, 
                          fileName = sampleFiles,
                          condition = sampleCondition,
			  factor = sampleFactor,
			  cell = sampleCell)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
#ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ] #pre-filtering to remove rows with only 0 or 1 read
# note that the "results" function below automatically performs independent filtering, so no need to do the filtering as above.
#ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref="NEGATIVE")
analysis = DESeq(ddsHTSeq)
resultsNames(analysis)
write.table(as.data.frame(counts(analysis,normalized=TRUE)),file="normalized_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")
write.table(as.data.frame(counts(analysis)),file="raw_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")

for (cell in c("E14","3F")) {
	for (protein in c("CTCF","SMC1","MED12","KLF4","FLAG","EZH2")) {
		write.table(as.data.frame(results(analysis, contrast=c("condition",paste0(protein,"_",cell),paste0("input","_",cell)))), file=paste0("DESeq2_",protein,"_",cell,"_vs_input_",cell,".txt"),quote=FALSE,row.names=TRUE,sep="\t")		
	}
}


#transformation
nt <- normTransform(ddsHTSeq) # defaults to log2(x+1)
rld <- rlog(ddsHTSeq, blind=FALSE)
#vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)

write.table(assay(nt),file="log2_transformed_normalized_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")
write.table(assay(rld),file="rld_transformed_normalized_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")


######################################################################
# data quality assessment
# Heatmap of the count matrix
select <- order(rowMeans(counts(analysis,normalized=TRUE)),decreasing=TRUE)[1:20] # 20 most highly expressed genes
#log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(ddsHTSeq)[,c("condition")])
rownames(df) <- colnames(rld)
colnames(df) <- "condition"

# 1. using default "euclidean" as clustering distance, and "complete" as clustering method.
pdf("regularized log transformation (euclidean complete).pdf")
pheatmap(assay(rld)[select,], clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
dev.off()

# 2. using "euclidean" as clustering distance, and "average" as clustering method. 
pdf("regularized log transformation (euclidean average).pdf")
pheatmap(assay(rld)[select,], clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "average", cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
dev.off()

# 3. using "correlation" as clustering distance, and "average" as clustering method. 
pdf("regularized log transformation (correlation average).pdf")
pheatmap(assay(rld)[select,], clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", clustering_method = "average", cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
dev.off()

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)), method="euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)
#rownames(sampleDistMatrix) <- paste(rld$condition, sep="") # note if there is "type" too
#colnames(sampleDistMatrix) <- paste(rld$condition, sep="") # note if there is "type" too

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("sample-to-sample distances (euclidean complete).pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, clustering_method = "complete")
dev.off()


# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)
#rownames(sampleDistMatrix) <- paste(rld$condition, sep="") # note if there is "type" too
#colnames(sampleDistMatrix) <- paste(rld$condition, sep="") # note if there is "type" too

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("sample-to-sample distances (euclidean average).pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, clustering_method = "average")
dev.off()


#sampleCor <- cor(assay(rld))
#sampleDist <- 1-sampleCor
#pdf("sample-to-sample distances (correlation average).pdf",height=8,width=8)
#pheatmap(sampleDist, clustering_distance_rows="correlation", clustering_distance_cols="correlation", clustering_method = "average")
#dev.off()



# Principal component plot of the samples
pdf("PCA.pdf",height=8,width=8)
data <- plotPCA(rld, intgroup=c("factor","cell"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=factor,shape=cell)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(axis.title=element_text(size=14))
dev.off()

pdf("PCA_default.pdf",height=8,width=8)
plotPCA(rld, intgroup=c("condition"))
dev.off()



######################################################################
# contrast: day7/day0
for (condition in c("ROSA")) {
        result = results(analysis, alpha=0.1, contrast=c("condition",condition,"NEGATIVE"))
        summary(result)

        resultOrdered <- result[order(result$padj),] #order our results table by the smallest adjusted p value

        #write.csv(as.data.frame(resultOrdered),file=paste("result_DESeq2_", condition, "_vs_day0.csv", sep=""))
	write.table(as.data.frame(result),file=paste("result_DESeq2_", condition, "_vs_NEGATIVE.txt", sep=""),quote=FALSE,row.names=TRUE,sep="\t")
        write.table(as.data.frame(resultOrdered),file=paste("result_DESeq2_", condition, "_vs_NEGATIVE_ordered.txt", sep=""),quote=FALSE,row.names=TRUE,sep="\t")

        #MA plot
        pdf(paste("MA_plot_", condition, "_vs_NEGATIVE.pdf", sep=""),height=8,width=8)
        plotMA(result, ylim=c(-2,2), alpha=0.1)
        dev.off()

        #count plot of the gene with min padj
        pdf(paste("count_plot_", condition, "_vs_NEGATIVE.pdf", sep=""),height=8,width=8)
        #plotCounts(ddsHTSeq, gene=which.min(result$padj), intgroup="condition", main=paste("count plot: ", condition, " vs. WT", sep=""))
        d <- plotCounts(ddsHTSeq, gene=which.min(result$padj), intgroup="condition", returnData=TRUE)
        print(ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)))
        dev.off()
}



library(corrplot)
M <- cor(assay(rld))

for (i in c(2:16)) {
	pdf(paste0("corr_hclust",i,".pdf"))
	corrplot(M, order="hclust", addrect=i, is.corr=FALSE, tl.col="black", tl.cex=0.7)
	dev.off()
}


library(ggplot2)
pdf("Klf4.pdf")
d <- plotCounts(ddsHTSeq,gene="Klf4", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition,y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + theme(axis.text.y=element_text(size=18), axis.text.x=element_text(size=18,angle=45,hjust=1), axis.title=element_text(size=18,face="bold")) + ggtitle("Klf4") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))
dev.off()

pdf("Sox2.pdf")
d <- plotCounts(ddsHTSeq,gene="Sox2", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition,y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + theme(axis.text.y=element_text(size=18), axis.text.x=element_text(size=18,angle=45,hjust=1), axis.title=element_text(size=18,face="bold")) + ggtitle("Sox2") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))
dev.off()

pdf("Pou5f1.pdf")
d <- plotCounts(ddsHTSeq,gene="Pou5f1", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition,y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + theme(axis.text.y=element_text(size=18), axis.text.x=element_text(size=18,angle=45,hjust=1), axis.title=element_text(size=18,face="bold")) + ggtitle("Pou5f1") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))
dev.off()

pdf("Myc.pdf")
d <- plotCounts(ddsHTSeq,gene="Myc", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition,y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + theme(axis.text.y=element_text(size=18), axis.text.x=element_text(size=18,angle=45,hjust=1), axis.title=element_text(size=18,face="bold")) + ggtitle("Myc") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))
dev.off()

pdf("Nanog.pdf")
d <- plotCounts(ddsHTSeq,gene="Nanog", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition,y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + theme(axis.text.y=element_text(size=18), axis.text.x=element_text(size=18,angle=45,hjust=1), axis.title=element_text(size=18,face="bold")) + ggtitle("Nanog") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))
dev.off()
