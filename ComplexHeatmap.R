setwd("Google Drive/Toshi_result/")
library("ComplexHeatmap")
library(dendextend)

#ROSA VS. NEGATIVE
data_DESeq <- read.table("alignment_sample/results_DESeq2.txt", header=T, comment.char="#")
rosa_genes <- row.names(data_DESeq[which(data_DESeq$log2FoldChange>0 & data_DESeq$padj<0.1),])
negative_genes <- row.names(data_DESeq[which(data_DESeq$log2FoldChange<0 & data_DESeq$padj<0.1),])

data_all <- read.table("alignment_sample/normalized_counts.txt", header=T)
data <- subset(data_all, row.names(data_all) %in% rosa_genes | row.names(data_all) %in% negative_genes)
data$label <- ""
for (i in 1:nrow(data)) {
  if (row.names(data[i,]) %in% rosa_genes) {
    data[i,]$label = 1
  }
  else if (row.names(data[i,]) %in% negative_genes) {
    data[i,]$label = 2
  }
  else {
    data[i,]$label = 3
  }
}

rnames <- row.names(data)
mat_data <- data.matrix(data[,1:6])
colnames(mat_data) <- colnames(data)[1:6]
rownames(mat_data) <- rnames        
base_mean = rowMeans(mat_data)
mat_scaled = t(apply(mat_data, 1, scale))
colnames(mat_scaled) <- colnames(data)[1:6]

#row_dend = hclust(dist(mat_scaled)) # row clustering
#col_dend = hclust(dist(t(mat_scaled))) # column clustering

type = gsub("_\\d+", "", colnames(mat_scaled))
ha = HeatmapAnnotation(df = data.frame(type = type))

# Define some graphics to display the distribution of columns
#.hist = anno_histogram(mat_scaled, gp = gpar(fill = "lightgreen"))
#.density = anno_density(mat_scaled, type = "line", gp = gpar(col = "blue"))
#ha_mix_top = HeatmapAnnotation(df = data.frame(type = type), hist = .hist, density = .density)

# Define some graphics to display the distribution of rows
#.hist = anno_histogram(mat_scaled, gp = gpar(fill = "lightyellow"), which="row")
#.density = anno_density(mat_scaled, type = "line", gp = gpar(col = "gold"), which="row")
#.violin = anno_density(mat_scaled, type = "violin", gp = gpar(fill = "lightgreen"), which = "row")
#.boxplot = anno_boxplot(mat_scaled, gp = gpar(fill = "deeppink"), which = "row")

ha_mix_right = HeatmapAnnotation(violin = .violin, bxplt = .boxplot, hist = .hist, density = .density, which = "row", width = unit(4, "cm"))

png("heatmaps1.png", width = 8*300, height = 8*300, res = 300)  
Heatmap(mat_scaled, 
        name = "mRNA", 
        split = data$label, 
        #km = 3,
        show_row_names = FALSE,
        show_row_dend = TRUE,
        column_names_gp = gpar(fontsize = 10),
        #cluster_rows = color_branches(row_dend, k = 3),
        #cluster_columns = color_branches(col_dend, k = 2),
        top_annotation = ha, 
        top_annotation_height = unit(4, "mm"),
        #col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
) +
  Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, column_names_gp = gpar(fontsize = 10), width = unit(5, "mm")) 
#+ ha_mix_right
dev.off()
densityHeatmap(scale(mat_scaled))



