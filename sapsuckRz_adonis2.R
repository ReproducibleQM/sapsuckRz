#Prepare data for NMMDS
library(cdlTools)
sites.crops0 <- sites.crops[which(rowSums(sites.crops[5:14])>0),] #removes rows with no aphid obs
env.nmds <- sites.crops0[,c(3:14)]

env.nmds <- env.nmds[,3:ncol(env.nmds)]
colnames(sites.crops0) <- updateNamesCDL(colnames(sites.crops0))


adonis. <- adonis(vegdist(env.nmds) ~ sites.crops0$Corn, permutations = 100)



aadonis.
