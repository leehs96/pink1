library(cmapR)
library(CePa)
setwd('//10.80.51.8/hslee/gtex')

#gtex_parse <- read.gct('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct')


gtex_tpm <- read.table(
  'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.txt',
  header = T,
  sep = '\t'
  )

gtex_parse <- system.file("extdata","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", package ="cmapR")
gtex <- parse_gctx(gtex_tpm)

gtex_anno <- read.csv('gtex_anno.csv', header = T, sep = ',')
