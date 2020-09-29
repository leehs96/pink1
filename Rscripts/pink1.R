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
files <- c("P2_rsem.genes.results",
           "P4_rsem.genes.results",
           "P6_rsem.genes.results",
           "P8_rsem.genes.results",
           "P9_rsem.genes.results",
           "W5_rsem.genes.results",
           "W6_rsem.genes.results",
           "W7_rsem.genes.results",
           "W8_rsem.genes.results",
           "W9_rsem.genes.results")
#add names to list
names(files) <- c("P2","P4","P6","P8","P9","W5","W6","W7","W8","W9")

#tximport datafiles in directory
txData <- tximport(files,'rsem')

#change seqLen 0s to 1s
txData$length[txData$length == 0] <-1

#read column data
columns <- read.table("column.txt", header = T, sep = "\t" )

#data structure
dds <- DESeqDataSetFromTximport(txData, columns , ~Status)

#remove lowly exp. gene
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

#export normalized exp values
ddsN <- estimateSizeFactors(dds)
ddsN <- estimateDispersions(ddsN)

Ndds <- counts(ddsN, normalized = T)

write.csv(Ndds,"pink1_norm.csv")

#DEG
ddsDE <- DESeq(dds)

#results change p-adj sig value 0.05
res <- results(ddsDE, contrast = c("Status","pink1 K/O", "WT") , alpha = 0.05)

#order by p-adj
resOrdered <- res[order(res$padj),]

write.csv(resOrdered,"desep_results.csv")

##Plotting

#plot mean exp by lfc
plotMA(ddsDE)

#PCA
vsd <- vst(dds, blind = F)
plotPCA(vsd, intgroup = c("Species", "Status"))

? plotP
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
         cluster_cols = F,
         cluster_rows = F,
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
statsDE <- read.csv("ttest.csv", row.names = 1)
statsDE <- na.omit(statsDE)
statsDE <- statsDE[order(statsDE$padj),]
write.csv(statsDE,"ttest1.csv")

statsDE <- read.csv("ttest1.csv", row.names = 1)
statsDE$sig <- ifelse(statsDE$padj <= 0.05, "sig", "not")

#plotMA replicate
plotMA <- ggplot(statsDE, aes(x = log10(baseMean) , y = log2FoldChange, color = sig)) + 
          geom_point() +
          theme(legend.position = "none")
plotMA

#Volcano
statsDE$genelabels <- rownames(statsDE) %in% rownames(statsDE[1:10,])



#set threshold
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

threshold <- statsDE$padj < padj.cutoff & abs(statsDE$log2FoldChange) > lfc.cutoff
length(which(threshold))
statsDE$threshold <- threshold
sigKO <- data.frame(subset(statsDE, threshold == TRUE & baseMean >= 100))

file <- read.csv("desep_results.csv",header = T, sep = ',')

threshold1 <- file$padj < padj.cutoff & abs(file$log2FoldChange) > lfc.cutoff & file$baseMean >= 100
length(which(threshold1))
file$threshold1 <- threshold1
sigKO1 <- data.frame(subset(file, threshold1 == TRUE))
write.csv(sigKO1,'DAVID.csv')



DownReg <- data.frame(subset(sigKO1, log2FoldChange < 0))
UpReg <- data.frame(subset(sigKO1, log2FoldChange >= 0))

write.csv(DownReg,'DAVID(Down).csv')
write.csv(UpReg,'DAVID(Up).csv')


volcano <- 
  ggplot(statsDE, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  geom_point() +
  geom_text_repel(
    colour = 'black',
    size = 3,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      label = ifelse(genelabels == T, rownames(statsDE),""))
    ) +
  ggtitle("Pink1 knock-out") +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
    )

volcano

#heatmap stuff

filtered <- read.csv("DAVID.csv", header = T, sep = ",")
top <-  filtered      #statsDE[1:500,]

test <- read.csv("pink1_norm.csv", header = T, sep = ",")
top <- merge(top, test, by = 1)
top <- top[order(top$padj),]
keep <- c('X','P2','P4','P6','P8','P9','W5','W6','W7','W8','W9')
top <- top[,names(top) %in% keep]

newTop <- top[,-1]

rownames(newTop) <- top[,1]

annotation <- data.frame(sampletype = rep(c("pink1 K/O","WT"),c(5,5)))
row.names(annotation) <- c('P2','P4','P6','P8','P9','W5','W6','W7','W8','W9')
pheatmap(log(newTop+1), 
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         annotation = annotation,
         border_color = NA,
         cluster_cols= FALSE
         )

#top10

top10 <- newTop[1:10,]
top10m <-melt(as.matrix(top10))
names(top10m) <- c("gene","sample","exp")
top10m$sample <- ifelse(grepl("P", top10m$sample),"pink1 K/O","WT")

top20 <- newTop[1:20,]
top20m <-melt(as.matrix(top20))
names(top20m) <- c("gene","sample","exp")
top20m$sample <- ifelse(grepl("P", top20m$sample),"pink1 K/O","WT")

order_data <- data.frame(first_column = c("Rab6b","Bcat1","Gpnmb","Stra6","Lama3","Pink1","Gm49980","Rbp1","Sox17","Hs3st4"), 
                         second_column = c(1,2,3,4,5,6,7,8,9,10))

top10m$gene <- factor(top10m$gene, levels = top10m$gene[order(order_data$second_column)])

geneExp <- ggplot(top10m, aes(x = sample, y = exp ,color = sample)) +
  scale_y_log10() +
  geom_point(size = 2) +
  facet_grid(~gene) +
  xlab(NULL) +
  ylab("Normalized counts") +
  ggtitle("Top10") +
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1))
geneExp


geneExp1 <- ggplot(top10m) +
  scale_y_log10() +
  geom_point(aes(x = gene, y = exp ,color = sample), size = 2) +
  xlab("Genes") +
  ylab("Normalized counts") +
  ggtitle("Top10 DE genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
geneExp1

geneExp20 <- ggplot(top20m) +
  scale_y_log10() +
  geom_point(aes(x = gene, y = exp ,color = sample), size = 2) +
  xlab("Genes") +
  ylab("Normalized counts") +
  ggtitle("Top20 DE genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
geneExp20


geneExp20_1 <- 
  ggplot(
    top20m, aes(x = sample, y = log2(exp+1),color = sample)) +
    geom_point(size = 2, shape = "square") +
    facet_grid(~gene) +
    xlab(NULL) +
    ylab("Log2 fold change") +
    ggtitle("Top20 DE genes") +
    theme(
      axis.text.x = element_text(
        face = "bold", 
        size = 10,
        angle = 45,
        hjust = 1
        )
      )

geneExp20_1

pheatmap(
  log2(top10+1), 
  show_rownames = T,
  show_colnames = F,
  border_color = NA,
  color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
  annotation = annotation
)

#mouse to human
musGenes <- read.csv('mGene.csv', header = T , sep = ',')
musGenes <- as.list(musGenes)

musGenes <-c("Hmmr", "Tlx3", "Cpeb4")
musGenes <- c("Rab6b","Bcat1","Gpnmb","Stra6","Lama3","Pink1","Gm49980","Rbp1","Sox17","Hs3st4")
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = x ,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human, uniqueRows=T
    )
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(humanx)
  return(humanx)
}


convertMouseGeneList(musGenes)

