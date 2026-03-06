

### Figure 1b
dt = read.table("../00.data/virus.info", sep="\t", header=T, check.names=F)

dt$quality <- ifelse( dt$completeness == "100", 
                      ifelse(  dt$quality_level == "ITR (high-confidence)" | dt$quality_level == "DTR (high-confidence)",
                               "High-confidence complete","Low-confidence complete"),
                      "Near-complete")


colors.complete = structure(c("#3bc9db", "#f7b731","#d9d9d9"),
                            names=c("High-confidence complete", "Low-confidence complete","Near-complete"))


source("/share/data1/zhangy2/scripts/R_my_functions/zy_pie.R")

pdt <- data.frame(table(dt$quality))
p <- zy_pie(pdt, value="Freq", fill="Var1", col=colors.complete)




### Figure 1c
dt = read.table("../00.data/virus.ckv.gz", sep="\t", header=F, check.names=F)
dt$quality <- ifelse( dt$V10 == "100", 
                      ifelse(  dt$V11 == "ITR (high-confidence)" | dt$V11 == "DTR (high-confidence)",
                               "High-confidence complete","Low-confidence complete"),
                      "Near-complete")


colors.complete = structure(c("#3bc9db", "#f7b731","#d9d9d9"),
                            names=c("High-confidence complete", "Low-confidence complete","Near-complete"))

p <- ggplot(dt, aes(x=V2/1000, fill=quality,color=quality))+
  geom_histogram(bins=100, alpha=0.6)+
  scale_x_log10(breaks=c(5, 6, 10, 20, 30, 60, 100, 200, 300, 600))+
  scale_y_continuous(breaks=c(2500, 5000, 7500))+
  scale_fill_manual(values=colors.complete)+
  scale_color_manual(values=colors.complete)+
  labs(x="Length (Kbp)", y = "Number of viral genoems")+
  theme_bw()+
  theme(panel.grid.minor.x=element_blank())

p

### Figure 1d
dt = read.table("./roc.tsv", header=T, check.names=F)

pdt <- dt %>%
  group_by(group, nspecies) %>%
  summarise(value = mean(nclust))
head(pdt)

p <- ggplot(pdt, aes(x=nspecies/1000, y=value/1000, color=group))+
  geom_line()+
  theme_bw()+
  labs(x="Number of viral genomes(x1000)", y="Number of vOTUs(x1000)")


## Figure 1e
source("/share/data1/zhangy2/scripts/R_my_functions/zy_pie.R")

dt = read.table("./oavgc.shared.pub.all.tsv", sep="\t", header=T)
p3 <- zy_pie(dt, value="value", fill="group")
p3


dt = read.table("./oavgc.shared.pub.all.tsv.detail", sep="\t", header=T)
dtf <- dt %>% arrange(shared)
dtf$name = factor(dtf$name, levels=unique(dtf$name))


p4 <- ggplot(dtf, aes(x=shared, y=name, fill=group))+
  geom_bar(stat='identity')+
  geom_text(aes(label=shared))+
  theme_bw()
p4


p <- ggpubr::ggarrange(plotlist=list(p1,p2,p3,p4,p5), ncol=2, nrow=3)
