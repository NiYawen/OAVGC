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

## 只挑选重要的几个部位
samp_map <- samp_map %>% filter(Sample_site %in% c("Saliva","Dental plaque","Tongue","BALF","Sputum","Oropharynx"))

dtf <- merge(samp_map, dt, by.x='Sample_ID', by.y='row.names') %>%
  filter(!Sample_site %in% c("Oral mix", "Airway mix")) %>%
  mutate(x = paste(BioProject_ID, "(",Sample_site,")", sep=""))


## 每个部位不同的项目数量
site.count <- unique(dtf[,c("BioProject_ID", "Sample_site")]) %>% group_by(Sample_site) %>% count() %>% rename(site.count = n)
dtf <- merge(dtf, site.count, by='Sample_site')

## 挑选大于20个样本的项目部位
proj.sample.count = dtf %>% group_by(x) %>% count() %>% filter(n >= 20) %>% rename(proj.sample.count = n)

dtf <- merge(dtf, proj.sample.count, by='x')

dtf <- dtf %>%
  arrange(desc(site.count), Sample_site, desc(proj.sample.count))

x = read.table("./sample_type.color.tsv", sep="\t", comment.char = "")


type.colors = x$V2
names(type.colors) = x$V1

pdt <- rbind(data.frame(x=dtf$x, y=dtf$obs, group="OBS", type=dtf$Sample_site, site.type=dtf$Sample_type),
             data.frame(x=dtf$x, y=dtf$shann, group="Shannon", type=dtf$Sample_site, site.type=dtf$Sample_type))

pdt$x = factor(pdt$x, levels=unique(dtf$x))
pdt$site.type = factor(pdt$site.type, levels=c("oral cavity", "airway"))

numFunc = function(x){
  if(x < 0.001){formatC(x, digits = 1, width = 1, format = "e", flag = "0")}
  else if(x<0.05){formatC(x, digits = 3, width = 1, format = "f", flag = "0")}
  else{NA}
}

xs = c("Saliva", "Dental plaque", "Tongue", "BALF","Sputum","Oropharynx")
comp = combn(as.character(unique(pdt$type)), 2, list)

pdt$type = factor(pdt$type, levels=xs)
pdt = pdt %>% arrange(type)
pdt$x = factor(pdt$x, levels=unique(pdt$x))

p1 <- ggplot(subset(pdt, group=="OBS"), aes(x=x, y=y, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=type.colors)+
  facet_grid(. ~ site.type, scales="free_x", space ="free_x" )+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))
  #scale_y_continuous(trans="log10")

p12 <- subset(pdt, group=="OBS") %>%
  group_by(type, x) %>%
  summarise(value=mean(y)) %>%
  ggplot(., aes(x=type, y=value, fill=type))+
  geom_boxplot()+
  geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=numFunc)+
  scale_fill_manual(values=type.colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))



p2 <- ggplot(subset(pdt, group=="Shannon"), aes(x=x, y=y, fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=type.colors)+
  facet_grid(. ~ site.type, scales="free_x", space ="free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))

p22 <- subset(pdt, group=="Shannon") %>%
  group_by(type, x) %>%
  summarise(value=mean(y)) %>%
  ggplot(., aes(x=type, y=value, fill=type))+
  geom_boxplot()+
  geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=numFunc)+
  scale_fill_manual(values=type.colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))


p <- ggpubr::ggarrange(plotlist=list(p1,p12,p2,p22), nrow=2, ncol=2, widths=c(3,1), common.legend = T)

p

# ggsave("votu.diversity.v2.pdf", p, width=15, height=12)

