setwd("/share/data1/zhangy2/project/63.elderlycohort.2025/aa.analyses/14.votu.alpha")
rm(list=ls())

myread <- function(inf){
  library(data.table)
  dt = fread(inf, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

dt = myread("../00.data/prok.profile")

data = data.frame(shann = vegan::diversity(t(dt), index = "shann"))
samp = read.table("../00.data/sample.info", sep="\t", header=T)
sampf <- samp %>% select(Sample_ID, NCBI_BioProject_ID)
dm <- merge(data, sampf, by.x='row.names', by.y='Sample_ID')
write.table(dm, "prok.diversity.tsv", sep="\t")
