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
LT_cancer <- read.csv("LT_cancer표기.csv", sep = ',', header = T)
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

#boxplot - normal

#normal_box <-ggplot(normal_all_LT, aes(x=LT, y=log2(exp_fpkm)+1, fill=LT)) + 
#              geom_boxplot() +
#              ggtitle("PINK1 expression in normal thyroid tissue") +
#              xlab('') +
#              ylab('normalized counts') +
#              theme(legend.position = "none")

normal_box <- ggboxplot(
  normal_all_LT,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
    ) +
  ylab("normalized counts") +
  xlab("normal tissue") +
  ggtitle("PINK1 expression")
normal_box


##statistical test

stat.test_norm1 <- normal_all_LT %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm1

              
stat.test_norm1 <- stat.test_norm1 %>% add_xy_position(x = LT)
#star

normal_box +
  stat_pvalue_manual(
    #stat.test_all[1,],
    stat.test_norm1,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

normal_box +
  stat_pvalue_manual(stat.test_norm1, label = "T-test, p = {p.adj}",
                     vjust = -1, bracket.nudge.y = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------------------------------------------------------#

#boxplot - tumor all

tumor1 <- select(tumor_all, status, exp_fpkm )
tumor_all_LT <- merge(tumor1, LT, by = 1)

tumor_all_LT$LT <- ifelse(grepl("1", tumor_all_LT$LT),"LT","no-LT")


tumor_box <- ggboxplot(
  tumor_all_LT,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("Thyroid cancer") +
  ggtitle("PINK1 expression")
tumor_box

##statistical test

stat.test_norm2 <- tumor_all_LT %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm2


stat.test_norm2 <- stat.test_norm2 %>% add_xy_position(x = LT)
#star

tumor_box +
  stat_pvalue_manual(
    stat.test_norm2,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
    ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

tumor_box +
  stat_pvalue_manual(
    stat.test_norm2,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
    ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


#-----------------------------------------------------------------------------#

#boxplot - tumor FA/miFTC

tumor2 <- select(FA_miFTC, status, exp_fpkm )
tumor_FA_LT <- merge(tumor2, LT, by = 1)

tumor_FA_LT$LT <- ifelse(grepl("1", tumor_FA_LT$LT),"LT","no-LT")



tumor2_box <- ggboxplot(
  tumor_FA_LT,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("FA/miFTC") +
  ggtitle("PINK1 expression")
tumor2_box

##statistical test

stat.test_norm3 <- tumor_FA_LT %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm3


stat.test_norm3 <- stat.test_norm3 %>% add_xy_position(x = LT)
#star

tumor2_box +
  stat_pvalue_manual(
    stat.test_norm3,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

tumor2_box +
  stat_pvalue_manual(
    stat.test_norm3,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------------------------------------------------------#

#boxplot - tumor FVPTC

tumor3 <- select(FVPTC, status, exp_fpkm )
tumor_FVPTC_LT <- merge(tumor3, LT, by = 1)

tumor_FVPTC_LT$LT <- ifelse(grepl("1", tumor_FVPTC_LT$LT),"LT","no-LT")


tumor3_box <- ggboxplot(
  tumor_FVPTC_LT,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("FVPTC") +
  ggtitle("PINK1 expression")
tumor3_box

##statistical test

stat.test_norm4 <- tumor_FVPTC_LT %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm4


stat.test_norm4 <- stat.test_norm4 %>% add_xy_position(x = LT)
#star

tumor3_box +
  stat_pvalue_manual(
    stat.test_norm4,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

tumor3_box +
  stat_pvalue_manual(
    stat.test_norm4,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------------------------------------------------------#

#boxplot - tumor PTC

tumor4 <- select(PTC, status, exp_fpkm )
tumor_PTC_LT <- merge(tumor4, LT, by = 1)

tumor_PTC_LT$LT <- ifelse(grepl("1", tumor_PTC_LT$LT),"LT","no-LT")


tumor4_box <- ggboxplot(
  tumor_PTC_LT,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("PTC") +
  ggtitle("PINK1 expression")
tumor4_box

##statistical test

stat.test_norm5 <- tumor_PTC_LT %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm5


stat.test_norm5 <- stat.test_norm5 %>% add_xy_position(x = LT)
#star

tumor4_box +
  stat_pvalue_manual(
    stat.test_norm5,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

tumor4_box +
  stat_pvalue_manual(
    stat.test_norm5,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------------------------------------------------------#

#boxplot - tumor FA only

tumor_FA_only_LT <- merge(tumor2, LT_cancer, by = 1)

tumor_FA_only_LT$LT <- ifelse(grepl("1", tumor_FA_only_LT$LT),"LT","no-LT")

FA_only <- data.frame(subset(tumor_FA_only_LT, cancer == "FA"))

FA_only_box <- ggboxplot(
  FA_only,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("FA only") +
  ggtitle("PINK1 expression")
FA_only_box

##statistical test

stat.test_norm_FA <- FA_only %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm_FA


stat.test_norm_FA <- stat.test_norm_FA %>% add_xy_position(x = LT)
#star

FA_only_box  +
  stat_pvalue_manual(
    stat.test_norm_FA,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

FA_only_box +
  stat_pvalue_manual(
    stat.test_norm_FA,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------------------------------------------------------#

#boxplot - tumor miFTC only


miFTC_only <- data.frame(subset(tumor_FA_only_LT, cancer == "miFTC"))

miFTC_only_box <- ggboxplot(
  miFTC_only,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("normalized counts") +
  xlab("miFTC only") +
  ggtitle("PINK1 expression")
miFTC_only_box

##statistical test

stat.test_norm_miFTC <- miFTC_only %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_norm_miFTC


stat.test_norm_miFTC <- stat.test_norm_miFTC %>% add_xy_position(x = LT)
#star

miFTC_only_box  +
  stat_pvalue_manual(
    stat.test_norm_miFTC,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

miFTC_only_box +
  stat_pvalue_manual(
    stat.test_norm_miFTC,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------------------------------------------------------#

##summary
#tumor_all_LT 

#normal_all_LT$cancer <- "normal" 
#tumor_FVPTC_LT$cancer <- "FVPTC" 
#tumor_PTC_LT$cancer <- "PTC"

normal_all_LT
tumor_FVPTC_LT
tumor_PTC_LT
FA_only
miFTC_only


all_box <- rbind(
  normal_all_LT,
  tumor_FVPTC_LT,
  tumor_PTC_LT,
  FA_only,
  miFTC_only
)

#ordering for facet
library(Hmisc)
library(gdata)

all_box$cancer <- factor(all_box$cancer, levels = c("normal","FA","miFTC","FVPTC","PTC"))

all_box_plot <- ggboxplot(
  all_box,
  x = "LT",
  y = "exp_fpkm",
  fill = "LT",
  facet.by = "cancer"
  ) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  facet_wrap(~ cancer, ncol = 5)+
  ylab("normalized counts") +
  xlab("") +
  ggtitle("PINK1 expression")
all_box_plot

##statistical test

stat.test_all <- all_box %>%
  group_by(cancer) %>%
  t_test(exp_fpkm ~ LT) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test_all

stat.test_all <- stat.test_all %>% add_xy_position(x = "LT")
#star

all_box_plot  +
  stat_pvalue_manual(
    stat.test_all,
    label = "p.adj.signif",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#p-value

miFTC_only_box +
  stat_pvalue_manual(
    stat.test_norm_miFTC,
    label = "T-test, p = {p}",
    vjust = -0.5,
    bracket.nudge.y = 2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))



#-----------------------------------------------------------------------------#
#comparisons ref

all_box_plot_1 <- ggboxplot(
  all_box,
  x = "cancer",
  y = "exp_fpkm",
  fill = "cancer",
  palette = c("#00798c","#d1495b","#edae49","#66a182","#2e4057")
) +
  ylab("normalized counts") +
  xlab("comparisons aginst reference(normal)") +
  ggtitle("PINK1 expression") +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )
  
all_box_plot_1
stat.test_all_1 <- all_box %>% t_test(exp_fpkm ~ cancer, ref.group = "normal") 
stat.test_all_1 <- stat.test_all_1 %>% add_xy_position(x = "cancer")
stat.test_all_1
all_box_plot_1 + stat_pvalue_manual(stat.test_all_1, label = "p.adj.signif", tip.length = 0.01)


#comparisons all(basemean)

all_box_plot_2 <- ggboxplot(
  all_box,
  x = "cancer",
  y = "exp_fpkm",
  fill = "cancer",
  palette = c("#00798c","#d1495b","#edae49","#66a182","#2e4057")
) +
  ylab("normalized counts") +
  xlab("comparisons aginst all(basemean)") +
  ggtitle("PINK1 expression") +
  theme(
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )


stat.test_all_2 <- all_box %>% t_test(exp_fpkm ~ cancer, ref.group = "all") 
stat.test_all_2 <- stat.test_all_2 %>% add_xy_position(x = "cancer")
stat.test_all_2
all_box_plot_2 + stat_pvalue_manual(stat.test_all_2, label = "p.adj.signif", tip.length = 0.01, y.position = 40 )
