library(edgeR)
library(GenomicFeatures)
library(RMariaDB)
library(tximport)

#setwd('C:/Users/LEE HAN SAI/Desktop/Work/Pink1/Rscripts')
#mus.ens <-  makeTxDbFromGFF('Mus_musculus.GRCm38.101.gtf')

#all.genes <- genes(mus.ens)
#genes.length <- width(all.genes)
#names(genes.length) <- all.genes$gene_id

#write.csv(genes.length, "Mus_musculus.GRCm38.genelength_ensembl.txt")

#read.table('Mus_musculus.GRCm38.genelength_ensembl.txt', header = TRUE, sep = "\t" )

setwd('C:/Users/LEE HAN SAI/Desktop/Work/Pink1/Rscripts')

#import htseq files
P2 <- read.table("P2_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
P4 <- read.table("P4_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
P6 <- read.table("P6_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
P8 <- read.table("P8_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
P9 <- read.table("P9_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
W5 <- read.table("W5_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
W6 <- read.table("W6_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
W7 <- read.table("W7_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
W8 <- read.table("W8_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')
W9 <- read.table("W9_HTseq_count.txt", sep = "\t", header = T, row.names = 'gene')

htseq <- cbind(P2,P4,P6,P8,P9,W5,W6,W7,W8,W9)


Total_counts <- htseq

gene_length <- read.csv(
  "Mus_musculus.GRCm38.genelength_ensembl.csv",
  sep = ",",
  header = T,
  row.names = 'gene'
)

y_all <- DGEList(Total_counts) 
y_all <- calcNormFactors(y_all)
rpkm_y_all <- rpkm(
  y_all,
  gene.length = gene_length[rownames(y_all),],
  normalized.lib.sizes = TRUE
)
rpkm_y_all <- round(rpkm_y_all, digits = 2)

#Output file
write.csv(
  rpkm_y_all,
  file = "pink1_HTSeq_edgeR_FPKM.csv", 
  sep = ",", 
  row.names = rownames(rpkm_y_all),
  col.names = colnames(rpkm_y_all), 
  quote = FALSE
)
