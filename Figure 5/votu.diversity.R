setwd("/share/data1/zhangy2/project/63.elderlycohort.2025/aa.analyses/14.votu.alpha")
rm(list=ls())

library(dplyr)
library(ggplot2)

myread <- function(inf){
  library(data.table)
  dt = fread(inf, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

# dt = myread("../00.data/votu.profile")
# alpha = data.frame(shann = vegan::diversity(t(dt)))
# obs = data.frame(obs = colSums(dt>0))
# diver <- data.frame(obs=obs, alpha=alpha)

# write.table(diver, "diversity_data.tsv", sep="\t", quote=F)

dt <- read.table("./diversity_data.tsv", sep="\t", header=T)
samp_map = read.table("../00.data/sample.info", sep="\t", header=T, check.names=F)

dtm <- merge(samp_map, dt, by.x='Sample_ID', by.y='row.names')

dtf <- dtm %>%
  filter(Sample_type == "airway" | Sample_site  == "Saliva") %>%
  mutate(x = paste(NCBI_BioProject_ID, "(",Sample_site,")", sep=""))

site.count <- unique(dtf[,c("NCBI_BioProject_ID", "Sample_site")]) %>% group_by(Sample_site) %>% count() %>% rename(site.count = n)
dtf <- merge(dtf, site.count, by='Sample_site')

## 挑选大于10个样本的项目
proj.sample.count = dtf %>% group_by(x) %>% count() %>% filter(n >= 10) %>% rename(proj.sample.count = n)

dtf <- merge(dtf, proj.sample.count, by='x')

dtf <- dtf %>%
  arrange(desc(site.count), desc(proj.sample.count))

x = read.table("./sample_type.color.tsv", sep="\t", comment.char = "")

type.colors = x$V2
names(type.colors) = x$V1

pdt <- rbind(data.frame(x=dtf$x, y=dtf$obs, group="OBS", type=dtf$Sample_site),
data.frame(x=dtf$x, y=dtf$shann, group="Shannon", type=dtf$Sample_site))

pdt$x = factor(pdt$x, levels=unique(dtf$x))



p1 <- ggplot(subset(pdt, group=="OBS"), aes(x=x, y=y, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=type.colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))+
  scale_y_continuous(trans="sqrt")

p2 <- ggplot(subset(pdt, group=="Shannon"), aes(x=x, y=y, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=type.colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))

p <- ggpubr::ggarrange(plotlist=list(p1,p2), ncol=1, common.legend = T)

ggsave("votu.diversity.pdf", p, width=12, height=12)

