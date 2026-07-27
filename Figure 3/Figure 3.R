
## Figure 3a
pdt <- data.frame(anno=c(202429,36929,69120), lab = c("high","low","unk"))
source("/share/data1/zhangy2/scripts/R_my_functions/zy_pie.R")

p1 <- zy_pie(pdt, value="anno", fill="lab")


## Figure 3b
dt = read.table("annoted.genes.all", sep="\t", header=T, check.names=F)
dt$rate = round(dt$anno / 68708, digits=2)
dt$label = paste(dt$rate, "%", sep="")
dt$db = factor(dt$db, levels=rev(dt$db))

p <- ggplot(dt, aes(x=rate, y=db))+
  geom_bar(stat="identity")+
  geom_text(aes(label=anno))+
  theme_bw()

## Figure 3c

## Figure 3d
library(circlize)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(tidyr) # extract



align_heat <- function(heat_matrix, select_names){
  no = setdiff(select_names, colnames(heat_matrix))
  heat_matrix[,no] = 0
  heat_matrix = heat_matrix[,select_names]
  heat_matrix[is.na(heat_matrix)] = 0
  
  a = heat_matrix[rowSums(heat_matrix)!=0,]
  
  return(a)
}


## lifestyle
virus.map = read.table("../00.data/votu.group", sep="\t", header = T)

family.style = virus.map %>%
  filter(family != "") %>%
  group_by(family, lifestyle) %>%
  count() %>%
  dcast(formula = family ~ lifestyle, value.var="n")


family.style[is.na(family.style)] = 0

family.style <- family.style %>%
  mutate(tot = rowSums(select(., Mix, Temperate, Virulent)),
         mix = Mix/tot*100,
         temp = Temperate / tot*100,
         vir = Virulent / tot * 100)
rownames(family.style) = family.style$family


#### 宿主
dt = read.table("./family.host.phylum-genus", sep="\t")

dtf = dt %>%
  mutate(rate = V3/V4) %>%
  filter(rate >= 0.5, V4>=80, V2 != "multiple") %>%
  extract(V2,c("phylum","genus"),"(.*)xxxx(.*)") %>%
  mutate(assign = V3, other = V4-V3)

dtf <- dtf %>%
  merge( dtf %>% group_by(phylum) %>% summarise(pc = sum(V3)), by='phylum') %>%
  merge( dtf %>% group_by(genus) %>% summarise(gc = sum(V3)), by='genus') %>%
  arrange(desc(pc), desc(gc), desc(V3))

host.color = read.table("./family.host.color", sep="\t", comment.char = "")

colors.taxo = host.color$V2
names(colors.taxo) = host.color$V1
  

colors.lifestyle = c("#d9d9d9", "#df5173", "#27537c" )
names(colors.lifestyle) = c("Mix", "Temperate", "Virulent")
colors.vir.rate = c("#666666", "#cccccc")
names(colors.vir.rate) = c("current taxo", "Others")

top_anno = HeatmapAnnotation(
  # lifestyle = anno_customize(rep("pie",86), graphics = graphics),
  lifestyle = anno_barplot(family.style[dtf$V1,c("mix","temp",'vir')], gp = gpar(fill=colors.lifestyle), which = 'row', bar_width=1),
  
  Novotu = anno_barplot(dtf[,c("V3","other")], gp = gpar(fill = colors.vir.rate ), which='row', bar_width = 1),
  
  genus = dtf$genus,
  genus_name = anno_text(dtf$genus),
  phylum = dtf$phylum,
  phylum_name = anno_text(dtf$phylum),
  
  annotation_height = unit(c(4, 4, 0.5, 4, 0.5, 4), "cm"),
  
  col = list(
    genus = colors.taxo,
    phylum = colors.taxo
  )
)

############ 主图
kegg = read.table("./family.kegg.pathway.levelC.metabolism", sep="\t", header=T, check.names=F, row.names=1, comment.char = "", quote = "")
vc = read.table("family.votu.count", sep="\t", row.names=1)
vc = vc[colnames(kegg), ]

## 计算比例
kegg = t(t(kegg)/vc)
heat_matrix = align_heat(kegg, dtf$V1)

##### 过滤，进行fisher检验
tmp1 = melt(heat_matrix)
tmp1$group = tmp1$value > 0
tmp2 <- merge(tmp1, dtf, by.x="Var2",by.y="V1") %>%
  group_by(Var1, group,phylum ) %>%
  count()

tmp.res = matrix(NA, nrow=nrow(heat_matrix), ncol=2, dimnames = list(rownames(heat_matrix), NULL))
for(mod in rownames(heat_matrix)){
  
  tmp.tmp <- subset(tmp2, Var1==mod)
  
  tmp.tmp1 <- dcast(tmp.tmp, phylum ~ group, value.var="n")
  tmp.tmp1[is.na(tmp.tmp1)] = 0
  
  tmp.tmp1 = tmp.tmp1[rowSums(tmp.tmp1[,c(2,3)])> 5,] # 筛选，至少有5个family的phylum才计算pvalue
  max_rate = max(tmp.tmp1$`TRUE` / tmp.tmp1$`FALSE`)
  
  p <- fisher.test(tmp.tmp1[,c(2,3)])$p.value
  tmp.res[mod,] = c(p, max_rate) 
}

show.modules <- as.data.frame(tmp.res) %>% filter(V1<0.05) # 筛选pvalue
show.modules <- as.data.frame(tmp.res) %>% filter(V2>=0.5)

heat_matrix = heat_matrix[rownames(show.modules),] # 至少在10个subfamily中有这个通路的基因才行

heat_matrix = apply(heat_matrix, 2, function(x){log10(x+1e-5)} )
hist(heat_matrix)

mycol = colorRamp2(breaks =c(0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5, 6 ), 
                   # colors=c("white", "#f0f9e8", "#a8ddb5", "#7bccc4", "#4eb3d3", "#08589e" )
                   colors = c("white", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081")
)

mycol = colorRamp2(breaks =c(-5, -4, -3, -2, -1, 1 ), colors=c("white","#ffffd9", "#c7e9b4",  "#41b6c4", "#225ea8", "#081d58" ))
nh = nrow(heat_matrix)
nw = ncol(heat_matrix)
myunit = 5
max_nchar = max(nchar(rownames(heat_matrix)))


levelc_desc = read.table("/share/data1/Database/KEGG/20230401/tmp/pathway.levelC.f", sep="\t", quote = "", comment.char = "", row.names=1)

right_anno = rowAnnotation(
  levelc = anno_text(levelc_desc[rownames(heat_matrix),'V2'])
)

ha <- Heatmap(heat_matrix,
              #row_split = as.factor(keggf[rownames(heat_matrix),"B"]),
              row_title_rot = 0,
              show_row_dend = F,
              top_annotation = top_anno,
              row_names_gp = gpar(fontsize=8),
              row_split = as.factor(levelc_desc[rownames(heat_matrix),'V3']),
              
              #clustering_distance_row = "binary",
              
              right_annotation = right_anno,
              column_split = as.factor(dtf$phylum),
              column_title_rot = 45,
              
              cluster_columns = F,
              col = mycol,
              row_gap=unit(0,'mm'), column_gap = unit(0,'mm'), border = T,
              rect_gp = gpar(col = "grey", lwd = 0.1), # 内部线条颜色
              height = nh * unit(myunit,"mm"), width = nw * unit(myunit,"mm"), # 保持单元格是方的
              row_names_max_width = unit(max_nchar, "char"), # 有时候legend和标签文字会重叠
              
              
)


pdf("heatmap.kegg_levelC.pdf", width=30 ,height=50)
draw(ha)
dev.off()

