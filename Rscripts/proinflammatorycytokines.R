statsDE <- read.csv("ttest1.csv", row.names = 1)
proinf <- read.csv("proinflammatorycytokines.csv", header = T)


statsDE$proInf <- rownames(statsDE) %in% proinf$cytokine

