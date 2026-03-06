setwd("/share/data1/zhangy2/project/63.elderlycohort.2025/aa.analyses/14.votu.alpha")
rm(list=ls())


alpha.vir = read.table("./votu.lifestyle.tsv", sep="\t")[,c('Row.names','x', 'obs','shann','lifestyle')]

alpha.prok = read.table("./prok.diversity.tsv", sep="\t")[,c('Row.names', 'shann')]
colnames(alpha.prok)[2] = "prok.shann"

samp <- read.table("../00.data/sample.info", sep="\t", header=T, check.names = F)[,c("Sample_ID",'Sample_site', "Sample_type")]

head(alpha.prok)
head(alpha.vir)
xs = c("Saliva", "Dental plaque", "Tongue", "BALF","Sputum","Oropharynx")

dm <- merge(alpha.vir, alpha.prok[,c(1,2)], by='Row.names') %>%
  merge(., samp, by.x='Row.names', by.y='Sample_ID', all.x=T)

dmf <- dm %>%
  filter(Sample_site %in% xs)

res.cor = rbind()
for (life in unique(dmf$lifestyle)){
  
  for (site in unique(dmf$Sample_site)){
    
    tmp.dt = subset(dmf, lifestyle==life & Sample_site ==site)
    x = tmp.dt$prok.shann
    y = tmp.dt$shann
    z = cor.test(x,y,method ="spearman")
    tmp = data.frame(site=site,lifestyle=life, pval=z$p.value, cor=z$estimate[[1]])
    res.cor = rbind(res.cor, tmp)
  }
}

p1 <- ggplot(data=dmf, aes(y=prok.shann, x=shann))+
  facet_wrap(lifestyle ~ Sample_site, scales='free', ncol=6)+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="Shannon vOTU", y="Shannon (Bacteria)")


all.vir = read.table("./diversity_data.tsv", sep="\t")
dm1 <- merge(alpha.prok, all.vir, by.x='Row.names', by.y='row.names') %>%
  merge(samp, by.x='Row.names', by.y='Sample_ID') %>%
  merge(unique(alpha.vir[,c("Row.names","x")]), by='Row.names') %>%
  filter(Sample_site %in% xs)

p2 <- ggplot(dm1, aes(y=prok.shann, x=shann))+
  facet_wrap(. ~ Sample_site, scales='free', ncol=6)+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x="Shannon vOTU", y="Shannon (Bacteria)")
p2

for (site in unique(dm1$Sample_site)){
  
  tmp.dt = subset(dm1, Sample_site ==site)
  x = tmp.dt$prok.shann
  y = tmp.dt$shann
  z = cor.test(x,y,method ="spearman")
  tmp = data.frame(site=site,lifestyle='all', pval=z$p.value, cor=z$estimate[[1]])
  res.cor = rbind(res.cor, tmp)
}

p <- ggpubr::ggarrange(plotlist=list(p1,p2), ncol=1, nrow=2, heights=c(2,1))
# write.table(res.cor,"alpha.virus.prok.pdf.corr.tsv", sep="\t")
# ggsave("alpha.virus.prok.pdf", p, width=12, height=7)
