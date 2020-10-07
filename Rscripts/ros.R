setwd("C:/Users/LEE HAN SAI/Desktop/Work/Pink1")
ros_exp <- read.csv("ttest1_upper.csv", row.names = 1)
ros_list <- read.csv("ros_signal.csv", header = T)


ros_exp$ros <- rownames(ros_exp) %in% ros_list$ros

ros_exp <- data.frame(subset(ros_exp, ros_exp$ros == TRUE))


#set threshold
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

threshold <- ros_exp$padj < padj.cutoff & abs(ros_exp$log2FoldChange) > lfc.cutoff

ros_exp$threshold <- threshold
ros_exp <- data.frame(subset(ros_exp, baseMean >= 100))
#volcano
library(ggplot2)
library(ggrepel)

ros_exp$genelabels <- TRUE

volcano <- 
  ggplot(ros_exp, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  geom_point(size = 4) +
  geom_text_repel(
    colour = 'black',
    size = 5,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      label = ifelse(threshold == T, rownames(ros_exp),""))
  ) +
  ggtitle("ROS related gene expression") +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10,2))

volcano


#boxplot, heatmap
library(pheatmap)

boxplot <- data.frame(subset(proinf_exp,threshold == T))

exp <- read.csv("pink1_norm_rmDUP_upper.csv", header = T, sep = ",", row.names = 1)


boxplot <- merge(boxplot, exp, by = 0)
boxplot <- boxplot[order(boxplot$padj),]
keep <- c('Row.names','P2','P4','P6','P8','P9','W5','W6','W7','W8','W9')
boxplot <- boxplot[,names(boxplot) %in% keep]

newBox <- boxplot[,-1]

rownames(newBox) <- boxplot[,1]

annotation <- data.frame(sampletype = rep(c("pink1 K/O","WT"),c(5,5)))
row.names(annotation) <- c('P2','P4','P6','P8','P9','W5','W6','W7','W8','W9')
pheatmap(log(newBox+1), 
         
         show_colnames = F,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         annotation = annotation,
         border_color = NA,
         cluster_cols= FALSE
)


newBoxm <-melt(as.matrix(newBox))
names(newBoxm) <- c("gene","sample","exp")
newBoxm$sample <- ifelse(grepl("P", newBoxm$sample),"pink1 K/O","WT")

ggplot(newBoxm, aes(x=gene, y=log2(exp)+1, fill=sample)) + 
  geom_boxplot(coef = 30) +
  ggtitle("Pro-inflammatory gene") +
  xlab("pro-inflammatory related gene") +
  ylab("normalized counts")
