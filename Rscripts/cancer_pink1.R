setwd("C:/Users/LEE HAN SAI/Desktop/Work/Pink1/cancer")
library(ggplot2)
library(reshape)
library(ggpubr)
library(rstatix)
library(dplyr)

# file input
fpkm <- read.csv("cancer_fpkm_pink.csv", sep = ',', header = T, row.names = 1)
list <- read.csv("cancer_list.csv", sep = ',', header = T, row.names = 1)
LT <- read.csv("LT.csv", sep = ',', header = T)
#fpkm <- rbind(fpkm,list)
#fpkm <- na.omit(fpkm)

melt_fpkm <- melt(as.matrix(fpkm))
melt_fpkm <- na.omit(melt_fpkm)
melt_list <- melt(as.matrix(list))
names(melt_fpkm) <- c("PINK1","status","exp_fpkm")

#splitting normal, tumor, subset 
normal_all <-  melt_fpkm[grep("N", melt_fpkm$status),]
tumor_all <- melt_fpkm[-grep("N", melt_fpkm$status),]
FA_miFTC <- melt_fpkm[grep("FA|FTC", melt_fpkm$status),]
FVPTC <- melt_fpkm[grep("FV", melt_fpkm$status),]
PTC <- melt_fpkm[grep("PTC", melt_fpkm$status),]

#normal LT vs no-LT

normal_all$status <- gsub(".N","", normal_all$status)
normal1 <- select(normal_all, status, exp_fpkm )
normal_all_LT <- merge(normal1, LT, by = 1)

normal_all_LT$LT <- ifelse(grepl("1", normal_all_LT$LT),"LT","no-LT")

##boxplot

#normal_box <-ggplot(normal_all_LT, aes(x=LT, y=log2(exp_fpkm)+1, fill=LT)) + 
#              geom_boxplot() +
#              ggtitle("PINK1 expression in normal thyroid tissue") +
#              xlab('') +
#              ylab('normalized counts') +
#              theme(legend.position = "none")

normal_box <- ggboxplot(normal_all_LT, x = "LT", y = "exp_fpkm", fill = "LT") + theme(legend.position = "none")
normal_box
ttt<- ToothGrowth

##statistical test

stat.test_norm <- normal_all_LT %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm

              
stat.test_norm <- stat.test_norm %>% add_xy_position(x = LT)
#star

normal_box +
  stat_pvalue_manual(stat.test_norm, label = "p.adj.signif") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

normal_box +
  stat_pvalue_manual(stat.test_norm, label = "T-test, p = {p}",
                     vjust = -1, bracket.nudge.y = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))



