setwd("/share/data1/zhangy2/project/63.elderlycohort.2025/aa.analyses/14.votu.alpha")
rm(list=ls())

source("/share/data1/zhangy2/scripts/R_my_functions/zy_alpha_diversity.R")
map_samp = read.table("../00.data/sample.info", sep="\t", header=T)

No.site <- map_samp %>%
  group_by(Sample_site) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(No.site = n)

No.pro <- map_samp %>%
  group_by(NCBI_BioProject_ID) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(No.pro = n)


sampf <- map_samp %>%
  mutate(x = paste(NCBI_BioProject_ID, Sample_site, sep=".")) %>%
  select(Sample_ID, NCBI_BioProject_ID, Sample_site, DataType, x) %>%
  merge(., No.pro, by='NCBI_BioProject_ID') %>%
  merge(., No.site, by='Sample_site') %>%
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

dt1 = myread("../00.data/votu.profile.Temperate")
dt2 =  myread("../00.data/votu.profile.Virulent")

plot_dt1 = get_xxx(dt1/100, sampf)
plot_dt2 = get_xxx(dt2/100, sampf)

plot_dt1$lifestyle = "Temperate"
plot_dt2$lifestyle = "Virulent"

plot_dt = rbind(plot_dt1, plot_dt2)

pdt <- plot_dt
pdt$x = factor(pdt$x, levels=unique(sampf$x))

# write.table(pdt, "votu.lifestyle.tsv", sep="\t")

## shannon 指数
p12 <- ggplot(pdt, aes(x=x, y=shann, fill=lifestyle) )+
  geom_boxplot()+
  theme_bw()+
  #facet_grid(. ~ DataType, scales="free_x")+
  scale_fill_manual(values=colors.lifestyle)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
p12

## obs 指数
p13 <- ggplot(pdt, aes(x=x, y=obs, fill=lifestyle) )+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=colors.lifestyle)+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle=90, hjust=1))
p13


p <- ggpubr::ggarrange(plotlist=list(p12+theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()),
                                     p13),
                       ncol=1, nrow=2,
                       common.legend = T, heights = c(2,3))
p

# ggsave(filename="votu.lifestyle.pdf", p, width=12, height=12)