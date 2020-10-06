statsDE <- read.csv("ttest1.csv", row.names = 1)
proinf <- read.csv("proinflammatorycytokines.csv", header = T)

statsDE$proInf <- ifelse(row.names(statsDE) == , "pro", "not")
