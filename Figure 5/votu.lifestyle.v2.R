setwd("/share/data1/zhangy2/project/63.elderlycohort.2025/aa.analyses/14.votu.alpha")
rm(list=ls())

source("/share/data1/zhangy2/scripts/R_my_functions/zy_alpha_diversity.R")
map_samp = read.table("../00.data/sample.info", sep="\t", header=T)
study.ord = read.table("study.order.tsv", sep="\t")

No.site <- map_samp %>%
  group_by(Sample_site) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(No.site = n)

No.pro <- map_samp %>%
  group_by(BioProject_ID) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(No.pro = n)


sampf <- map_samp %>%
  mutate(x = paste(BioProject_ID, Sample_site, sep=".")) %>%
  select(Sample_ID, BioProject_ID, Sample_site, DataType, x) %>%
  merge(., No.pro, by='BioProject_ID', all.x=T) %>%
  merge(., No.site, by='Sample_site', all.x=T) %>%
  arrange(Sample_site, desc(No.pro), desc(No.site)) %>%
  filter(No.pro >= 20)


get_xxx <- function(dt.virus, map_samp){
  
  ids = intersect(colnames(dt.virus), map_samp$Sample_ID)
  map_samp = subset(map_samp, Sample_ID %in% ids)
  obs = data.frame(obs = colSums((dt.virus>0)+0))
  
  shann = apply(dt.virus, 2, function(x){-sum(x*log(x), na.rm = T)})
  shann = data.frame(shann = shann)
  
  mydiver = merge(obs, shann, by='row.names')
  plot_dt = merge(mydiver, map_samp, by='Row.names', by.y='Sample_ID')
  
  plot_dt
}

myread <- function(inf){
  library(data.table)
  dt = fread(inf, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}


load("../00.data/color.RData")

# dt1 = myread("../00.data/votu.profile.Temperate")
# dt2 =  myread("../00.data/votu.profile.Virulent")
# 
# plot_dt1 = get_xxx(dt1/100, sampf)
# plot_dt2 = get_xxx(dt2/100, sampf)
# 
# plot_dt1$lifestyle = "Temperate"
# plot_dt2$lifestyle = "Virulent"
# 
# plot_dt = rbind(plot_dt1, plot_dt2)
# 
# pdt <- plot_dt
# # write.table(pdt, "votu.lifestyle.tsv", sep="\t")


pdt = read.table("votu.lifestyle.tsv", sep="\t")
pdt$x = paste(pdt$BioProject_ID,"(", pdt$Sample_site,")", sep="")

pdt$x = factor(pdt$x, levels=study.ord$V1)
pdt = subset(pdt, !is.na(x))

numFunc = function(x){
  if(x < 0.001){formatC(x, digits = 1, width = 1, format = "e", flag = "0")}
  else if(x<0.05){formatC(x, digits = 3, width = 1, format = "f", flag = "0")}
  else{NA}
}

## shannon 指数

xs = c("Saliva", "Dental plaque", "Tongue", "BALF","Sputum","Oropharynx")
comp = combn(xs, 2, list)
pdt <- pdt %>%
  filter(Sample_site %in% xs) %>%
  mutate(Sample_site = factor(Sample_site, levels=xs))


p11 <- ggplot(pdt, aes(x=x, y=shann, fill=lifestyle) )+
  geom_boxplot()+
  theme_bw()+
  facet_grid(. ~ lifestyle, scales="free_x")+
  scale_fill_manual(values=colors.lifestyle)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
p11

p12 <- pdt %>%
  group_by(x, Sample_site, lifestyle) %>%
  summarise(sum_shann = mean(shann)) %>%
  ggplot(., aes(x=Sample_site, y=sum_shann, fill=lifestyle))+
  geom_boxplot()+
  geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=numFunc,tip_length = 0)+
  
  facet_grid(. ~ lifestyle, scales="free_x")+
  scale_fill_manual(values=colors.lifestyle)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(x=NULL, y="Mean shannon")

## obs 指数
p21 <- ggplot(pdt, aes(x=x, y=obs, fill=lifestyle) )+
  geom_boxplot()+
  theme_bw()+
  facet_grid(. ~ lifestyle, scales="free_x")+
  scale_fill_manual(values=colors.lifestyle)+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle=90, hjust=1))


p22 <- pdt %>%
  group_by(x, Sample_site, lifestyle) %>%
  summarise(value = mean(obs)) %>%
  ggplot(., aes(x=Sample_site, y=value, fill=lifestyle))+
  geom_boxplot()+
  geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=numFunc,tip_length = 0)+
  facet_grid(. ~ lifestyle, scales="free_x")+
  scale_fill_manual(values=colors.lifestyle)+
  scale_y_log10(trans='log10', breaks=c(1,10,100,1000))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(x=NULL, y="Mean OBS")

p <- ggpubr::ggarrange(plotlist=list(p11,
                                     p12,
                                     p21,
                                     p22),
                       ncol=2, nrow=2,
                       common.legend = T, widths=c(4,1), heights = c(1,1))
p

# ggsave(filename="votu.lifestyle.pdf", p, width=15, height=8)