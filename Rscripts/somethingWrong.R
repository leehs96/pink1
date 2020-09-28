#Mouse thyroid analysis

setwd("C:/Users/LEE HAN SAI/Desktop/work")
source("pink1.R")
library(DESeq2)
library(ggrepel)
library(tximport)
library(RColorBrewer)
library(pheatmap)

##file handling

#read files
files <- c("W9_rsem.genes.results",
           "W8_rsem.genes.results",
           "W7_rsem.genes.results",
           "W6_rsem.genes.results",
           "W5_rsem.genes.results",
           "P9_rsem.genes.results",
           "P8_rsem.genes.results",
           "P6_rsem.genes.results",
           "P4_rsem.genes.results",
           "P2_rsem.genes.results")
#add names to list
names(files) <- c("W9","W8","W7","W6","W5","P9","P8","P6","P4","P2")

#tximport datafiles in directory
txData <- tximport(files,'rsem')

#change seqLen 0s to 1s
txData$length[txData$length == 0] <-1

#read column data
columns <- read.table("column1.txt", header = T, sep = "\t" )

#data structure
dds <- DESeqDataSetFromTximport(txData, columns , ~Status)

#remove lowly exp. gene
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

#export normalized exp values
ddsN <- estimateSizeFactors(dds)
ddsN <- estimateDispersions(ddsN)

Ndds <- counts(ddsN, normalized = T)

write.csv(Ndds,"pink1_norm1.csv")

#DEG
ddsDE <- DESeq(dds)

#results change p-adj sig value 0.05
res <- results(ddsDE, contrast = c("Status","pink1 K/O", "WT") , alpha = 0.05)

#order by p-adj
resOrdered <- res[order(res$padj),]

write.csv(resOrdered,"ttest.csv")

##Plotting

#plot mean exp by lfc
plotMA(ddsDE)

#PCA
vsd <- vst(dds, blind = F)
plotPCA(vsd, intgroup = c("Species", "Status"))

#sample to sample dist.
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- vsd$Status
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
#heat-map sample sample distance
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         border_color = "white",
         display_numbers = T,
         cutree_rows = 2,
         cutree_cols = 2,
         color = colors
)

#ggplot

library(ggplot2)
library(reshape)

# exp file
normExp <- read.csv("pink1_norm_rmDUP.csv", row.names = 1)

#DESeq
statsDE <- read.csv("desep_results_4.csv", row.names = 1)
statsDE <- na.omit(statsDE)
statsDE <- statsDE[order(statsDE$padj),]
write.csv(statsDE,"desep_result_5.csv")

statsDE <- read.csv("desep_results_4.csv", row.names = 1)
statsDE$sig <- ifelse(statsDE$padj <= 0.05, "sig", "not")

#plotMA replicate
plotMA <- ggplot(statsDE, aes(x = log10(baseMean) , y = log2FoldChange, color = sig)) + 
  geom_point() +
  theme(legend.position = "none")
plotMA

#Volcano
statsDE$genelabels <- rownames(statsDE) %in% rownames(statsDE[1:10,])

volcano <- ggplot(statsDE, aes(x = log2FoldChange, y = -log10(padj), color = sig)) + 
  geom_point() +
  geom_text_repel(colour = 'black', size = 3, aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(statsDE),""))) +
  ggtitle("pink1 knock-out") +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  theme(legend.position = "none", plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

volcano

#heatmap stuff
top <- statsDE[1:500,]

write.csv(top,"top.csv")
top <- read.csv("top.csv")

top <- merge(top, normExp, by = 0)
top <- top[order(top$padj),]
keep <- c('Row.names','P2','P4','P6','P8','P9','W5','W6','W7','W8','W9')
top <- top[,names(top) %in% keep]
newTop <- top[,-1]
newTop <- newTop[order(newTop$padj),]

rownames(newTop) <- top[,1]

annotation <- data.frame(sampletype = rep(c("pink1 K/O","WT"),c(5,5)))
row.names(annotation) <- c('P2','P4','P6','P8','P9','W5','W6','W7','W8','W9')
pheatmap(log2(newTop+1), 
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         annotation = annotation,
         annotation_colors =,
         border_color = NA,
)

#top10

top10 <- newTop[1:10,]
top10m <-melt(as.matrix(top10))
names(top10m) <- c("gene","sample","exp")
top10m$sample <- ifelse(grepl("P", top10m$sample),"pink1 K/O","WT")



geneExp <- ggplot(top10m, aes(x = sample, y = log2(exp+1),color = sample)) +
  geom_point(size = 2, shape = "square") +
  facet_grid(~gene) +
  xlab(NULL) +
  ylab("Log2 fold change") +
  ggtitle("Top10") +
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1))
geneExp



pheatmap(log2(top10+1), 
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         annotation = annotation
)
