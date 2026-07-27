library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(vegan)
library(Maaslin2)
library(pROC)
library(survival)
library(survminer)

#Figure 7a-c---------------------------------------------------------------------
vartype <- read.csv("vartype.csv")
dtvir = read.csv("dtvir.example.csv",row.names=1,check.names=F)
dist.dtvir <- read.csv("dist.dtvir.example.csv",row.names=1,check.names=F)

var_factor <- filter(vartype,discrete=="1") %>% .$variable
var_number <- filter(vartype,discrete=="0") %>% .$variable
for (col in var_factor) {  
  if (col %in% names(meta)) {  
    meta[[col]] <- as.factor(meta[[col]])  
  }  
} 

for (col in var_number) {  
  if (col %in% names(meta)) {  
    meta[[col]] <- as.numeric(meta[[col]])  
  }  
} 
str(meta)

align_data=function(dt,dist.dt,metadata){
  tmp=intersect(colnames(dt),rownames(metadata))
  metadata = metadata[tmp,]
  dist.dt = dist.dt[tmp,tmp]
  dt =dt[,tmp] %>% .[rowSums(.)>0,]
  return(list(dt = dt, dist.dt = dist.dt, metadata = metadata))
}

op.meta <- meta %>% select(-c(ID,Sputumsample)) %>% filter(!is.na(OPsample)) %>% column_to_rownames("OPsample")
op.list <- align_data(dt=dtvir,dist.dt=dist.dtvir,metadata=op.meta)
op.meta <- op.list$metadata; op.dist.vir=op.list$dist.dt;op.vir=op.list$dt

spu.meta <- meta %>% select(-c(ID,OPsample)) %>% filter(!is.na(Sputumsample)) %>% column_to_rownames("Sputumsample")
spu.list <- align_data(dt=dtvir,dist.dt=dist.dtvir,metadata=spu.meta)
spu.meta <- spu.list$metadata; spu.dist.vir=spu.list$dist.dt;spu.vir=spu.list$dt

##adonis after adjustment
desired_combinations <- list(  
  c("op.meta", "op.dist.vir"),  
  c("spu.meta", "spu.dist.vir")  
)  

confounders <- c("Age","Age_class","Gender","BMI","Weight_class","Season_of_sampling",
                 "Year_of_sampling","Month_of_sampling","Antibiotics","Antiviral_drugs")

get_adjusted_r2 <- function(adonis_object) {
  n_observations = tail(adonis_object$Df,1)+1 
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
  adjusted_r2
}

for (combination in desired_combinations){
  #combination = desired_combinations[[1]]
  k <- combination[1]   
  j <- combination[2]  
  print(paste("Processing:", k, j))  
  metaf <- eval(parse(text = k))  
  dist.dt <- eval(parse(text=j)) %>% as.dist()
  
  res = data.frame(matrix(NA, ncol=4, nrow=ncol(metaf), dimnames = list(colnames(metaf), c('name','r2','pvalue','adj2'))))
  ci = 1
  for(i in colnames(metaf)){
    #i= colnames(metaf)[2]
    if (i %in% confounders) next 
    print(paste("Running:", i))
    #print(ci)
    x = metaf[,i]
    if(length(unique(x, na.rm=T)) == 1){ci = ci+1; next}
    ados = adonis2(dist.dt ~ metaf[,i]+metaf[,"Age"]+metaf[,"Gender"]+
                     metaf[,"BMI"]+metaf[,"Season_of_sampling"]+
                     metaf[,"Antibiotics"]+metaf[,"Antiviral_drugs"],by="terms",
                   parallel=40, na.action = na.omit) 
    r2 = ados$R2[1]
    p = ados$`Pr(>F)`[1]
    adjr = get_adjusted_r2(ados)
    res[i,1] = i  
    res[i,2:4] = c(r2,p,adjr)
    ci = ci + 1
  }
  res_new = res
  write.csv(res_new, file=paste("table",j,"adonis_adjust.csv", sep="_"),quote=F, row.names=F)
}

op.ado <- read.csv("table_op.dist.vir_adonis_adjust.csv")
spu.ado <- read.csv("table_spu.dist.vir_adonis_adjust.csv")

###subset R2
all_combinations <- list(  
  c("op.meta", "op.ado","op.dist.vir"),  
  c("spu.meta", "spu.ado","spu.dist.vir")  
) 

for (com in all_combinations){
  #com = all_combinations[[1]]
  i = com[1]
  j = com[2]
  k = com[3]
  print(paste("Processing:", i, j, k))  
  metaf <- eval(parse(text=i))
  res_new <- eval(parse(text=j))
  dist.dt <- eval(parse(text=k))
  
  res_plot <- res_new %>% subset(.$name %in% colnames(meta))
  res_plot <-  res_plot %>% filter(pvalue < 0.05) 
  res_p <- merge(res_plot,vartype,by.x="name",by.y="variable") %>% arrange(subset,adj2) 
  
  vsub <- unique(res_p$subset)
  res = data.frame(matrix(NA, ncol=4, nrow=length(vsub), dimnames = list(vsub, c('name','r2','pvalue','adj2'))))
  ci = 1
  for (m in vsub){
    #m= vsub[1]
    var <- filter(res_p,subset == m) %>% .[,"name"]
    formula= paste0("dist.dt~",paste(var,collapse="+"))
    print(formula)
    ados = adonis2(as.formula(formula), data=metaf, by=NULL, parallel=40, na.action = na.omit) 
    r2 = ados$R2[1]
    p = ados$`Pr(>F)`[1]
    adjr = get_adjusted_r2(ados)
    res[m,1] = m  
    res[m,2:4] = c(r2,p,adjr)
    ci = ci + 1
  }
  res_new = res
  write.csv(res_new, file=paste("subset_r2",j,"adonis_adjust.csv", sep="_"),quote=F, row.names=F)
}

##plot
var.uri <- filter(vartype,subset=="Respiratory_Tract_Infection") %>% pull(variable)
ado.sig.sub <- ado.sig[!ado.sig %in% var.uri]
op.res.sub <- filter(op.res,name %in% ado.sig.sub)
spu.res.sub <- filter(spu.res,name %in% ado.sig.sub)

ado.sig.uri <- ado.sig[ado.sig %in% var.uri]
op.res.uri <- filter(op.res,name %in% ado.sig.uri)
spu.res.uri <- filter(spu.res,name %in% ado.sig.uri)

var.sub.order <-c("Clinical Indices","Diseases","Elderly Health","Environment","Oral Hygiene","Lifestyle",
                  "Medicine","Respiratory Tract Infection","Demography","Symptom","Vaccines")
color.subset <- brewer.pal(length(var.sub.order),"Paired")
color.subset <- setNames(color.subset,var.sub.order)

plot_adonis_bar_adjusted <- function(plot_data,
                                     r2_show = "adj2",
                                     type = "type",
                                     pval = "pvalue",
                                     name = "name",
                                     color.map = NULL,
                                     color_col = "subset",
                                     title = NULL,
                                     name_order=NULL,
                                     star_x=NULL,
                                     type_x=NULL,
                                     limit_x=NULL
) {
  plot_data[[r2_show]] <- ifelse(plot_data[[r2_show]] < 0, 0, 100 * plot_data[[r2_show]])
  
  if(is.null(name_order)){
    plot_data <- plot_data %>%
      arrange(desc(is.na(.data[[r2_show]])), .data[[r2_show]])
    name_order <- plot_data[[name]]
  }else{
    name_order <- name_order[name_order %in% plot_data[[name]]]
  }
  plot_data[[name]] <- factor(plot_data[[name]], levels = name_order)
  
  plot_data <- plot_data %>%
    mutate(p_star = case_when(
      .data[[pval]] < 0.01 ~ "**",
      .data[[pval]] < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  max_val <- max(plot_data[[r2_show]], na.rm = TRUE)
  if(is.null(star_x)){
    star_x <- max_val * 1.05  
  }
  if(is.null(limit_x)){
    limit_x <- max_val * 1.25 
  }
  
  p <- ggplot(plot_data, aes(y = .data[[name]])) + 
    geom_segment(aes(x = 0, xend = .data[[r2_show]], 
                     yend = .data[[name]], color = .data[[color_col]]), size = 1.2) +
    geom_point(aes(x = .data[[r2_show]], color = .data[[color_col]]), size = 3) +
    geom_text(aes(x = star_x, label = p_star), 
              hjust = 0, vjust = 0.7, color = "black", size = 5) +
    expand_limits(x = limit_x) +
    coord_cartesian(clip = "off") + 
    scale_color_manual(values = color.map) +
    scale_fill_manual(values = color.map) +
    labs(x = expression("Effect size (adjusted R"^2*",%)"), y = NULL, title = title) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      legend.position = "right",
      plot.margin = margin(5, 50, 5, 5) 
    )
  
  return(list(plot = p, name_order = name_order,
              star_x =star_x,type_x= type_x,limit_x=limit_x))
  
}

p.op <- plot_adonis_bar_adjusted (plot_data=op.res.sub,
                                  r2_show = "adj2_adjusted",
                                  pval = "pvalue_adjusted",
                                  name = "name",
                                  color.map = color.subset,
                                  color_col = "subset",
                                  title = "OP adonis") 

p.spu <- plot_adonis_bar_adjusted (plot_data=spu.res.sub,
                                   r2_show = "adj2_adjusted",
                                   pval = "pvalue_adjusted",
                                   name = "name",
                                   color.map = color.subset,
                                   color_col = "subset",
                                   title = "Sputum adonis",
                                   name_order=p.op$name_order,
                                   star_x=p.op$star_x,
                                   limit_x=p.op$limit_x) 
p <- p.op$plot + p.spu$plot+plot_layout(nrow=1)

###subset
p.op.subset <- plot_adonis_bar_adjusted (plot_data=op.res.subset.sub,
                                         r2_show = "adj2_adjusted",
                                         pval = "pvalue_adjusted",
                                         name = "name",
                                         color.map = color.subset,
                                         color_col = "name",
                                         title = "OP adonis subset") 

##spu
p.spu.subset <- plot_adonis_bar_adjusted (plot_data=spu.res.subset.sub,
                                          r2_show = "adj2_adjusted",
                                          pval = "pvalue_adjusted",
                                          name = "name",
                                          color.map = color.subset,
                                          color_col = "name",
                                          title = "Sputum adonis subset",
                                          name_order=p.op.subset$name_order,
                                          star_x=p.op.subset$star_x,
                                          limit_x=p.op.subset$limit_x) 
p.subset <- p.op.subset$plot + p.spu.subset$plot+plot_layout(nrow=1)

#Figure 7d----------------------------------------------------------------------
plotdata <- read.csv("plot_r2_data.csv")

color.r2 <- c("#A6CEE3","#1F78B4","#FB9A99","#FDBF6F","#CAB2D6", "#808001","burlywood3")
pr <- ggplot(plotdata,aes(x=source,y=mean_r2,fill=type))+
  geom_bar(stat="identity",position="dodge",width=0.6)+
  theme_bw()+
  scale_fill_manual(values=color.r2)+
  labs(fill="type")+
  xlab("")+
  ylab(expression("RFVC R"^2))

#Figure 7e----------------------------------------------------------------------
host.stat <- result.anno %>% group_by(host.genus,source) %>% summarise(count=n( ))
host_plot <- host.stat %>%
  mutate(host.genus = ifelse(count == 1, "Other genera", host.genus)) %>%
  group_by(source, host.genus) %>%
  summarise(count = sum(count), .groups = "drop") 

nodes_unique <- nodes %>% 
  distinct(host.genus.label, nodes.color) 
color.host <- setNames(nodes_unique$nodes.color, nodes_unique$host.genus.label)

plist <- list()
plotlist <- c("Lung function","Elderly health")
for (i in plotlist){
  #i="Lung function"
  plotsub <- host_plot %>% filter(source==i) %>% 
    arrange(count) %>%
    mutate(host.genus = factor(host.genus, levels = c("Other genera",setdiff(host.genus, "Other genera"))))
  
  color.map <- color.host[as.character(unique(plotsub$host.genus))]
  
  p <- ggplot(plotsub, aes(x = host.genus, y = count,fill=host.genus)) +
    geom_col() +
    scale_fill_manual(values=color.map)+
    coord_flip() +
    labs(x = "Host genus", y = "Number of vOTUs",title=i) +
    theme_bw(base_size = 14)+
    theme(legend.position = "none") 
  plist[[i]] <- p
}
p <- wrap_plots(plist)

#Figure 7f----------------------------------------------------------------------
maaslin_result = function(mapdt=NA,taxdt=NA,ID="Sample",group="URTI_in_the_next_12_months",
                          covar=c("Age","BMI","Gender"),reference = NULL,
                          random_effects = NULL,factor_var=c("URTI_in_the_next_12_months","Gender"),
                          numeric_var = c("Age","BMI"),compare_groups = list(c("A", "B"), c("A", "C")),
                          min.abun=0.01, min.count=NULL,min.prev=0.1,group_factor=T,ncore=8,
                          dir="maaslin/",
                          name="op.urti"){
  
  if (!is.null(covar)) {
    group_all <- c(group,covar)
  }else{
    group_all <- group
  }
  
  if (!is.null(random_effects)) {
    group_all_r <- c(group_all,random_effects)
  }else{
    group_all_r <- group_all
  }
  
  mapdt <- mapdt[,c(ID,group_all_r)] %>% na.omit()
  
  taxfile <- filter_species(tax=taxdt,sample_map=mapdt,group=group,ID=ID,
                            min.abun=min.abun,min.count=min.count,min.prev=min.prev,group_factor = group_factor)
  
  meta <- taxfile$sample_map
  dtfil <- taxfile$tax_new
  print(colSums(dtfil))
  print(mean(colSums(dtfil)))
  
  dt.input <- t(dtfil) %>% as.data.frame()
  meta.input <- meta %>% column_to_rownames(ID)
  
  if (!is.null(factor_var)) {
    meta.input[factor_var] <- lapply(meta.input[factor_var], as.factor)
  }
  
  if (!is.null(numeric_var)) {
    meta.input[numeric_var] <- lapply(meta.input[numeric_var], as.numeric)
  }
  print(str(meta.input))
  
  #MASLIN2
  fit_data <- Maaslin2(
    input_data = dt.input, 
    input_metadata = meta.input, 
    output = dir, 
    min_abundance = 0,
    min_prevalence = 0, 
    random_effects = random_effects,
    fixed_effects = group_all, 
    normalization = "NONE", 
    reference = reference, 
    plot_heatmap = F, plot_scatter = F,  
    analysis_method = "LM", 
    correction = "BH",
    standardize = T, 
    cores = ncore, 
    transform = "LOG"  
  )
  
  res = fit_data$results
  resf = subset(res, metadata==group)
  tmpname <- data.frame(feature.name=colnames(dt.input),feature=make.names(colnames(dt.input)))
  resf1 <- resf %>% left_join(tmpname,by="feature")
  
  if(!is.null(compare_groups)){
    abundance_means <- dt.input %>%
      rownames_to_column(var = ID) %>%
      left_join(meta[, c(ID, group)], by = ID) %>% 
      dplyr::group_by(!!sym(group)) %>%
      dplyr::summarise(
        across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),.groups = "drop"
      ) %>%
      pivot_longer(-!!sym(group), names_to = "feature.name", values_to = "mean_abundance") %>%
      pivot_wider(names_from = !!sym(group), values_from = mean_abundance)
    
    # fold change
    for (comp in compare_groups) {
      if (length(comp) != 2) next
      
      group1 <- comp[1]
      group2 <- comp[2]
      
      fd_col <- paste0("FC_", group1, "_vs_", group2)
      fc_col  <- paste0("log2_FC_", group1, "_vs_", group2)
      enr_col <- paste0("Enriched_", group1, "_vs_", group2)
      
      abundance_means <- abundance_means %>%
        dplyr::mutate(
          !!fd_col :=(!!sym(group1) + 1e-6) / (!!sym(group2) + 1e-6),
          !!fc_col := log2((!!sym(group1) + 1e-6) / (!!sym(group2) + 1e-6)),
          !!enr_col :=
            ifelse(
              !!sym(fc_col) > 0,
              paste0(group1, "-enriched"),
              paste0(group2, "-enriched")
            )
        )
    }
    resf2 <- resf1 %>% left_join(abundance_means,by="feature.name")
  }else{
    resf2 <- resf1
  }
  
  write.csv(resf2, file = file.path(dir, paste0("maaslin_", name, ".csv")))
  return(list(result=resf2,metadata=meta,tax=dtfil))
}

op.ma <- maaslin_result(mapdt=meta.old.dt,taxdt=dtvir,ID="OPsample",group=group.var,
                        covar=c("Age","BMI","Gender","Center","Season_of_sampling","Antibiotics","Antiviral_drugs"),reference = c("Season_of_sampling,0"),
                        random_effects = NULL,factor_var=c(group.var,"Gender","Center","Season_of_sampling","Antibiotics","Antiviral_drugs"),
                        numeric_var = c("Age","BMI"),compare_groups = list(c("1", "0")),
                        min.abun=0.01, min.count=NULL,min.prev=0.1,group_factor=T,ncore=8,
                        dir=dir.path.op,
                        name=name.op)
op.result <- op.ma$result

spu.ma <- maaslin_result(mapdt=meta.old.dt,taxdt=dtvir,ID="Sputumsample",group=group.var,
                         covar=c("Age","BMI","Center","Gender","Season_of_sampling","Antibiotics","Antiviral_drugs"),reference = c("Season_of_sampling,0"),
                         random_effects = NULL,factor_var=c(group.var,"Gender","Center","Season_of_sampling","Antibiotics","Antiviral_drugs"),
                         numeric_var = c("Age","BMI"),compare_groups = list(c("1", "0")),
                         min.abun=0.01, min.count=NULL,min.prev=0.1,group_factor=T,ncore=8,
                         dir=dir.path.spu,
                         name=name.spu)
spu.result <- spu.ma$result

#plot
host.summ.op <- filter(op.result,pval < 0.05) %>% left_join(vir.anno,by=c("feature"="votu"))
host.stat.op <- host.summ.op %>% group_by(Enriched_1_vs_0,host.genus) %>% summarise(count=n( ))
host_plot.op <- host.stat.op %>%
  mutate(host.genus = ifelse(count == 1, "other genera", host.genus)) 

plotdata <- host_plot.op %>%
  group_by(host.genus) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  arrange(total_count) %>%
  mutate(host.genus = factor(host.genus, levels = c("other genera",setdiff(host.genus,"other genera"))))

total_species = sum(plotdata$count)
title=paste0("OP (n=", total_species, ")")

p1 <- ggplot(plotdata,
             aes(x = host.genus,
                 y = count,
                 fill = Enriched_1_vs_0)) +
  geom_col(position = "stack") +
  coord_flip() +
  labs(x = "Host genus",
       y = "Number of vOTUs",
       title = title,
       fill = "Enrichment") +
  scale_fill_manual(values=c("#a6cee3","#fb9a99"))+
  theme_bw(base_size = 14)
p1

#Figure 7g----------------------------------------------------------------------
map_name <- function(roc.list){
  oc = c() # old name
  nc = c() # new name
  for(rb in names(roc.list)){
    # b = signif(ci(roc.list[[rb]], of="auc")*100, digits=3)
    b = sprintf("%0.1f",ci(roc.list[[rb]], of="auc")*100) 
    c = paste(rb, " (",b[2], "%)\t95% CI: " , b[1],"%-",b[3],"%", sep="") 
    oc = c(oc,rb)
    nc = c(nc,c)
  }
  names(nc) = oc 
  nc
}
plot_roc <- function(
    dt, pred=NA, true=NA, group=NULL, levels=NA,
    fill=FALSE, title=NA,
    cols = NA, conf_level=0.95, boot_n=2000
){
  old_scipen <- getOption("scipen")
  old_digits <- getOption("digits")
  options(scipen=0)
  options(digits=7)
  
  roc.list <- list()
  
  if(is.null(group)){  
    if(sum(is.na(cols))==1) cols <- "darkblue"
    if(sum(is.na(levels))==1) levels <- unique(dt[[true]])
    
    y_true <- factor(dt[[true]], levels = levels)
    y_pred <- as.numeric(dt[[pred]])
    
    roc.list[['AUC']] <- roc(response = y_true, predictor = y_pred, levels = levels, quiet=TRUE)
    
  } else { 
    if(length(group) > 1){  
      old_group <- group
      dt <- my_join(dt, group, " - ")
      group <- "zy_tmp_group"
    }
    grps <- unique(dt[[group]])
    if(sum(is.na(cols))==1) cols <- 1:length(grps)
    
    for(g in grps){
      temp_dt <- dt[dt[[group]]==g, , drop=FALSE] %>% droplevels()
      
      y_true <- factor(temp_dt[[true]], levels = levels)
      y_pred <- as.numeric(temp_dt[[pred]])
      
      roc.list[[as.character(g)]] <- roc(response = y_true, predictor = y_pred, levels = levels, quiet=TRUE)
    }
  }
  
  new_name_map <- map_name(roc.list)
  
  p <- ggroc(roc.list) +
    theme_bw() +
    geom_segment(data = data.frame(x=0, y=1), 
                 aes(x=x, y=y, xend=1, yend=0),
                 color="#d9d9d9", lwd=.4, inherit.aes=FALSE) +
    scale_color_manual(values=cols, labels=new_name_map) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid = element_line(linetype="dashed", color="black", linewidth=0.2),
      panel.border = element_rect(color="black", linewidth=0.5),
      axis.ticks = element_line(color="black", linewidth=0.5),
      axis.ticks.length = unit(2,"mm"),
      axis.text = element_text(color="black"),
      legend.title = element_blank()
    ) +
    labs(title = title) +
    theme(axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.title=element_text(size=20),
          plot.title=element_text(size=24),
          legend.text = element_text(size=20))
  
  if(isTRUE(fill)){
    if(!"ci.se" %in% ls("package:pROC")){
      stop("pROC::ci.se() not found. Please install/upgrade pROC >= 1.18")
    }
    ci.list <- lapply(roc.list, function(rocobj){
      pROC::ci.se(rocobj, specificities=seq(0,1,0.01), conf.level=conf_level, boot.n=boot_n) %>%
        data.frame(check.names=F) %>%
        data.table::setDT(keep.rownames=TRUE)
    })
    data_ci <- dplyr::bind_rows(ci.list, .id="plot_group")
    data_ci$rn <- as.numeric(data_ci$rn)
    
    
    p <- p +
      ggplot2::geom_ribbon(data=data_ci,
                           aes(x = rn, ymin=`2.5%`, ymax=`97.5%`, fill=plot_group),
                           inherit.aes=FALSE, alpha=0.3, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values=cols)
  }
  options(scipen=old_scipen)
  options(digits=old_digits)
  
  return(list(plot=p, ROC=roc.list, labels=new_name_map))
}

p.vir <- plot_roc(
  dt=virpred, pred="mean_prob", true="label", group="Type", levels=c("1","0"),
  fill=FALSE, title=NULL,
  cols = NA, conf_level=0.95, boot_n=2000
)

#Figuer 7h-j--------------------------------------------------------------------
survival_event_plot <- function(
    event.dt,
    event_sample_col = "Sample",
    prob.dt,
    prob_col = "prob",
    prob_sample_col = "Sample",
    title = NULL,
    color.map = NULL
){
  medscore <- median(prob.dt[[prob_col]], na.rm = TRUE)
  prob.dt$score <- ifelse(
    prob.dt[[prob_col]] > medscore,
    "high_risk",
    "low_risk"
  )
  
  colnames(prob.dt)[match(prob_sample_col, colnames(prob.dt))] <- "Sample"
  colnames(event.dt)[match(event_sample_col, colnames(event.dt))] <- "Sample"
  
  surv_data <- left_join(prob.dt, event.dt, by = "Sample")
  
  surv_data$Months <- as.numeric(surv_data$Months)
  surv_data$Event  <- as.numeric(surv_data$Event)
  surv_data$score <- factor(surv_data$score,levels = c("low_risk", "high_risk"))

  km_fit <- survfit(Surv(Months, Event) ~ score, data = surv_data, model = TRUE)
  pval_res <- surv_pvalue(km_fit, data = surv_data)
  
  cox_fit <- coxph(Surv(Months, Event) ~ score, data = surv_data)
  cox_sum <- summary(cox_fit)
  
  p.forest <- ggforest(cox_fit, data = as.data.frame(surv_data), main = title)
  p.forest
  
  HR_res <- tibble(
    HR = round(cox_sum$coef[1, "exp(coef)"], 3),
    CI_lower = round(cox_sum$conf.int[1, "lower .95"], 3),
    CI_upper = round(cox_sum$conf.int[1, "upper .95"], 3),
    pvalue = signif(cox_sum$coef[1, "Pr(>|z|)"], 3)
  )
  
  if (is.null(color.map)) {
    color.map <- c("#45bcec", "#f591b6")
  } else {
    color.map <- color.map[1:2]
  }
  
  pl <- ggsurvplot(
    km_fit,
    data = surv_data,
    fun = "event",
    conf.int = FALSE,
    xlab = "Time (months)",
    ggtheme = theme_pubr(),
    pval = round(pval_res$pval, 3),
    title = title,
    palette = color.map
  )

  return(list(
    HR = HR_res,
    km_fit=km_fit,
    forest=p.forest,
    plot = pl
  ))
}

op.vir.test <- survival_event_plot(
  event.dt=event.new.op.dt,
  event_sample_col = "Sample",
  prob.dt=prob.op.vir,
  prob_col = "mean_prob",
  prob_sample_col = "Sample",
  title = paste0("OP virome ",group.var),
  color.map = NULL
)

spu.vir.test <- survival_event_plot(
  event.dt=event.new.spu.dt,
  event_sample_col = "Sample",
  prob.dt=prob.spu.vir,
  prob_col = "mean_prob",
  prob_sample_col = "Sample",
  title = paste0("Spu virome ",group.var),
  color.map = NULL
)

