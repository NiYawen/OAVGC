library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(vegan)
library(metafor)
library(pROC)
library(ggpubr)
library(ggrepel)
library(randomForest)
library(ComplexHeatmap)
library(circlize)

meta.disease <- read.csv("meta.disease.csv",check.names=F)
vir.disease.example <- read.csv("vir.disease.example.csv",check.names=F,row.names=1)

color.disease <- c("#ee9d9e","#9ee4d8","#f7e29b","#c2b1d9","#ae2233","#aac7f2","#f1cda4",
                   "#f0d8e1","#4e728f","#e28c68","#a5cbb3","#cb889a")
disease_order <- c("Dental caries","Periodontitis","Upper respiratory tract infection", 
                   "Pneumonia","COVID19","COPD",
                   "Cystic fibrosis","Pancreatic ductal carcinoma","Colorectal cancer",
                   "Rheumatoid arthritis","Type 1 diabetes mellitu","Autism spectrum disorder")
color.disease <- setNames(color.disease,disease_order)

type_order=c("Saliva","Dental plaque","Subgingival plaque","Tongue","Nasopharynx" ,"Tonsils",
             "Cough swabs","Oropharynx","Sputum")

meta.disease$Disease_name <- factor(meta.disease$Disease_name,levels= disease_order)
meta.disease$Sample_site <- factor(meta.disease$Sample_site,levels=type_order)
meta.disease <- meta.disease[order(meta.disease$Disease_name,meta.disease$Sample_site),]
project_order <- unique(meta.disease$Project_samplesite)

#Figure 6a-------------------------------------------------------------------------
disease.diversity.change=function(dt=NA,sample_map=NA, ID="Sample",index="shannon",
                                  Group="Group",Project="Project",sample.color=NA,title=NA,
                                  SampleType="SampleType",DiseaseType="disease_type",
                                  disease_order=disease_order,type_order=NA,dir=NA){
  
  sample_map <- filter(sample_map, !is.na(!!sym(Group)))
  inter <- intersect(colnames(dt),sample_map[,ID])
  if (length(inter) == 0) {  
    stop("Error: No overlapping columns found between dt and sample_map.")  
  } 
  dt <- dt[,colnames(dt) %in% inter]
  dt = dt[rowSums(dt)!=0,]
  sample_map <- filter(sample_map, sample_map[,ID] %in% inter)
  
  
  #alpha
  alpha=data.frame()
  
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }
  
  if(tolower(index) == "shannon"){
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  if(tolower(index) == "total_abund"){
    alpha = data.frame(alpha=colSums(dt))
  }  
  
  alpha <- alpha %>% rownames_to_column(ID)
  dm = merge(alpha,sample_map, by=ID)
  

  vprojtype <- unique(sample_map[,Project])
  res_list = list()
  for (proj in vprojtype){
    #proj=vprojtype[1]
    print(proj)
    tmpf <- dm %>% filter(!!sym(Project) == proj) %>% unique()
    if (length(unique(tmpf[[Group]])) > 1) {
      d1 = filter(tmpf,!!sym(Group) == "Control") %>% pull(alpha)
      d2 = filter(tmpf,!!sym(Group) != "Control") %>% pull(alpha)
      p = wilcox.test(d1,d2)$p.value
      res <- data.frame(Project=proj, pval=p)
      res_list = append(res_list, list(res))
    }
  }
  res_pval <- do.call("rbind", res_list)
  
  #foldchange
  dm[[Group]][dm[[Group]] != "Control"] <- "Disease"
  
  data.fc <- dm %>%
    dplyr::group_by(!!sym(Project),!!sym(Group)) %>%
    dplyr::summarise(alpha = mean(alpha)) %>%
    pivot_wider(names_from = !!sym(Group),values_from = alpha) %>% 
    dplyr::mutate(fold_change=(Disease-Control)/Control,
                  text_pos = ifelse(fold_change < 0, fold_change-0.1, fold_change+0.1)) %>%
    merge(res_pval, by.x=Project,by.y="Project") %>% 
    left_join(distinct(dm[,c(Project,DiseaseType,SampleType)]),by=Project)
  
  if(!is.na(dir)){
    write.csv(data.fc,paste(dir,"table.csv",sep="."))
  }
  
  #data.fc[,DiseaseType] <- gsub("_"," ",data.fc[,DiseaseType])
  data.fc[,DiseaseType] <- factor(data.fc[,DiseaseType],levels= disease_order)
  data.fc[,SampleType] <- factor(data.fc[,SampleType],levels=type_order)
  data.fc <- data.fc[order(data.fc[,DiseaseType],data.fc[,SampleType]),]
  project_order <- unique(data.fc[,Project])
  data.fc[,Project] <- factor(data.fc[,Project],levels=rev(project_order))
  
  data.fc <- data.fc %>%
    mutate(shape = ifelse(pval < 0.01, "**",ifelse(pval < 0.05, "*", NA)))
  data.fc$fold_change <- data.fc$fold_change*100
  data.fc$text_pos <- data.fc$text_pos*100
  
  
  alpha_plot <- ggplot(data=data.fc, aes(x=fold_change, y=!!sym(Project), color=!!sym(DiseaseType), fill=!!sym(DiseaseType)))+
    # geom_segment(aes(x=0,xend=fold_change, y=!!sym(Project),yend=!!sym(Project)))+
    # geom_point(size=3)+
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black", fill=NA),
          text=element_text(size=20))+
    geom_text(aes(label=shape, x=text_pos), nudge_y=-0.2,size=10,color="black")+ 
    scale_color_manual(values=color.disease)+
    scale_fill_manual(values=color.disease)+
    ggtitle(paste(title,"(%increase/decrease)",sep="\n"))+
    ylab("")+xlab("")
  
  return(alpha_plot)
}


#vir.shannon
p1 <- disease.diversity.change(dt=vir.disease,sample_map=meta.disease,ID="Sample_ID",index="shannon",Group="Disease",
                               Project="Project_samplesite",sample.color=color.disease,title="Shannon of vir",
                               SampleType="Sample_site",DiseaseType="Disease_name",disease_order=disease_order,
                               type_order=type_order,dir="shannon.vir")

p2 <- disease.diversity.change(dt=vir.disease,sample_map=meta.disease,ID="Sample_ID",index="obs",Group="Disease",
                               Project="Project_samplesite",sample.color=color.disease,title="Obs of vir",
                               SampleType="Sample_site",DiseaseType="Disease_name",disease_order=disease_order,
                               type_order=type_order,dir="shannon.vir")

alpha_plot_vir <- wrap_plots(p1,p2, nrow=1)


#Figure 6b----------------------------------------------------------------------
meta_metafor <- function(dt, sample="Sample", group = "group", group_pair = c("Disease", "Control"),
                         proj = "proj", measure = "SMD", method = "REML") {

  dt <- dt %>% 
    dplyr::rename(proj = all_of(proj), group = all_of(group), sample=all_of(sample))
  
  feature <- setdiff(colnames(dt), c("sample", "group", "proj"))
  meta_outp <- rbind()
  x = 1
  nfeature = length(feature)
  for (i in feature) {
    #i=feature[1]
    cat("\r", x,"/",nfeature)
    x = x+1
    tib <- dt %>% 
      subset(group%in%group_pair[1]) %>% 
      dplyr::select(all_of(i), proj) %>% 
      dplyr::rename(index = all_of(i)) %>% 
      dplyr::group_by(proj) %>% 
      dplyr::summarise(d_Mean = mean(index), d_Sd = sd(index), d_N = n())
    tib2 <- dt %>% 
      subset(group%in%group_pair[2]) %>% 
      dplyr::select(all_of(i), proj) %>% 
      dplyr::rename(index = all_of(i)) %>% 
      dplyr::group_by(proj) %>% 
      dplyr::summarise(c_Mean = mean(index), c_Sd = sd(index), c_N = n())

    meta_in <- merge(tib, tib2, by= "proj")
    smd_meta <- escalc(measure = measure, data = meta_in, append = T,
                       m1i = d_Mean, m2i = c_Mean, 
                       sd1i = d_Sd, sd2i = c_Sd, 
                       n1i = d_N, n2i = c_N)
    non_na <- smd_meta %>% dplyr::filter(!is.na(yi)) %>% nrow()
    if(non_na !=0 ){
      smd_rma = tryCatch({
        rma(yi, vi, method = method, data = smd_meta)
      }, error = function(e){
        message("reset")
        rma(yi, vi, method = method, data = smd_meta, control=list(stepadj=0.5, maxiter=1000))
      })
      
      
      # merge each data
      smd_meta <- smd_rma$data %>% 
        add_column(measure = measure,
                   model = "RM",  
                   method_tau2 = method,
                   val_tau2 = as.numeric(smd_rma$tau2),
                   I2 = paste0(round(smd_rma$I2, digits = 2), "%"), 
                   Q = smd_rma$QE, 
                   Q_pval = smd_rma$QEp, 
                   feature = i,
                   estimate = as.numeric(smd_rma$beta),
                   ci_lb = smd_rma$ci.lb,
                   ci_ub = smd_rma$ci.ub,
                   pval = smd_rma$pval)
      meta_outp <- rbind(meta_outp, smd_meta)
    }else{
      message("\nWarning: ",i," is not suitable.")
    }
  }
  return(meta_outp)
  cat(paste0("== estimate > 0, ==> ", group_pair[1]," ==\n== estimate < 0, ==> ",group_pair[2], " =="))
}

alpha_meta_forest <- function(dt, sample_map,
                              ID="Sample_ID",
                              Group="Disease",
                              Project="Project_samplesite",
                              index=c("shannon","obs"),
                              output_dir="meta_alpha_forest",
                              title="Meta"
){
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  inter <- intersect(colnames(dt), sample_map[[ID]])
  if(length(inter)==0) stop("No overlapping samples between dt and sample_map")
  dt <- dt[, inter, drop=FALSE]
  sample_map <- sample_map[sample_map[[ID]] %in% inter, ]
  
  alpha_df <- data.frame(
    Sample  = colnames(dt),
    shannon = vegan::diversity(t(dt), index = "shannon"),
    obs     = colSums(dt > 0)
  )
  write.csv(alpha_df, file=file.path(output_dir,"alpha_df.csv"), row.names=FALSE)
  
  meta_input <- left_join(
    sample_map[, c(ID, Group, Project)],
    alpha_df,
    by = setNames("Sample", ID)
  )
  
  meta_res <- meta_metafor(dt=meta_input, sample=ID, group=Group,
                           group_pair=c("Disease","Control"), proj=Project,
                           measure="SMD", method="REML")
  
  write.csv(meta_res, file=file.path(output_dir,"meta_alpha_results.csv"), row.names=FALSE)
  
  forest_df <- meta_res %>%
    dplyr::select(feature, proj, yi, ci_lb, ci_ub, estimate) %>%
    arrange(feature, estimate)
  
  p <- ggplot(forest_df, aes(y=feature, x=yi)) +
    geom_point(aes(x=yi), color="grey40", shape=1, size=3) +          # 每个项目效应值
    geom_segment(aes(x=ci_lb, xend=ci_ub, y=feature, yend=feature),
                 color="grey40") +
    geom_point(aes(x=estimate), color="black", shape=18, size=3) +    # 汇总效应量
    geom_vline(xintercept = 0, linetype="dashed", color="red") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text=element_text(size=12)) +
    labs(x="Hedges' g (SMD)", y="", title=title)
  
  ggsave(file.path(output_dir,"meta_alpha_forest.pdf"), p, width=6, height=6)
  
  return(list(meta_res=meta_res, plot=p))
}

p.virus <- alpha_meta_forest(
  dt = vir.infection,
  sample_map = meta.infection,
  ID="Sample_ID",
  Group="Disease",
  Project="Project_samplesite",
  index=c("shannon","obs"),
  title="Virus",
  output_dir="alpha_meta_analysis/virus/"
)

#figure 6c----------------------------------------------------------------------
align_dt_sample <- function(dt, sample_map, ID=NA){
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dt = dt[,sample_map[,ID]] %>% filter(rowSums(.) !=0)
  list(dt=dt, sample_map=sample_map)
}

zy_pcoa <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NULL,
                    ado_method="bray", pca_method="bray",shape_by = NULL,  #点的形状映射
                    levels=0.95, star_plot=F, ellipse_plot=T,
                    cut_rate = NA, cut_num=NA,adjust=NULL,
                    title="PCoA", x=1, y=2, ados=T, mydist=NULL){
  sample_map <- filter(sample_map,!is.na(!!sym(group)))
  if(typeof(mydist) != "NULL"){
    fmt_profile = align_dist_sample(mydist, sample_map, ID=ID)
    mydist = fmt_profile$dist
  }else{
    fmt_profile = align_dt_sample(dt, sample_map, ID=ID)
    dt = fmt_profile$dt
  }
  sample_map = fmt_profile$sample_map
  
  
  if(is.finite(cut_rate) || is.finite(cut_num)){
    sample_map = filter_group(sample_map, group=group, cut_rate = cut_rate, cut_num=cut_num)
  }
  
  ## colors 
  if ( typeof(sample.color) == "NULL" ){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(new_label=paste(!!sym(group), " (", count, ")", sep=""))
  
  new_label <- structure(group_summ$new_label,names=as.character(unlist(group_summ[,group])))
  
  message(paste(length(unique(sample_map[,group])), "of groups to plot"))
  
  if(typeof(mydist) == "NULL"){
    mydist = vegdist(t(dt), method = pca_method)
  }else{
    mydist = as.dist(mydist)
  }
  
  
  ado_r2 = ado_p = NA
  if (isTRUE(ados)){
    if(length(unique(sample_map[,group])) > 1){
      ## adonis
      ado = adonis2(mydist ~ sample_map[,group])
      ado_r2 = round(ado$R2[1]*100, digits = 1)  
      ado_p = ado$`Pr(>F)`[1]
      
      substitle <- paste0("'R'^2~'='~'", ado_r2, "%'~~italic('p')~'='~'", ado_p, "'") %>% 
        as.formula() %>% 
        eval()
    }
  }
  
  if(!is.null(adjust)){
    ## adonis-adjust
    ado.adj = adonis2(mydist ~ sample_map[,adjust]+sample_map[,group], by = "terms",parallel=40)
    if ("sample_map[, group]" %in% rownames(ado.adj)){
      ado_r2.adj = round(ado.adj$R2[1]*100, digits = 1) 
      ado_p.adj = ado.adj$`Pr(>F)`[2]
      substitle <- paste0(
        "R² = ", ado_r2, "%, p = ", ado_p, "\n",
        "After adjusted: R² = ", ado_r2.adj, "%, p = ", ado_p.adj)
    }
  }
  
  ## PCoA
  pcoa = cmdscale(mydist, k=min(10, dim(mydist)[1]-1), eig=T)
  eigs = signif(pcoa$eig/sum(pcoa$eig), 4)*100
  point = pcoa$points
  
  colnames(point) = paste("pcoa.", 1:ncol(point),sep="")
  
  xlab = paste("PCoA", x, " (",eigs[x],"%)", sep="")
  ylab = paste("PCoA", y, " (",eigs[y],"%)", sep="")

  group = ifelse(length(unique(sample_map[,group])) == 1, sample.color[1], group)
  dm = merge(point, sample_map, by.x='row.names', by.y=ID)
  
  if (!is.null(shape_by)) {
    stopifnot(shape_by %in% colnames(dm))
  }
  
 
  if (is.null(shape_by)) {
    
    p1 <- ggscatter(
      data = dm,
      x = paste("pcoa.", x, sep = ""),
      y = paste("pcoa.", y, sep = ""),
      color = group,
      star.plot = star_plot,
      ellipse.level = levels,
      ellipse = ellipse_plot
    )
    
  } else {
    shape_levels <- unique(dm[[shape_by]])
    shape_values <- c(0,1,2,3,4)[seq_along(shape_levels)]
    
    p1 <- ggscatter(
      data = dm,
      x = paste("pcoa.", x, sep = ""),
      y = paste("pcoa.", y, sep = ""),
      color = group,
      fill  = NA,
      shape = shape_by,      
      star.plot = star_plot,
      ellipse.level = levels,
      ellipse = ellipse_plot
    )+
      scale_shape_manual(
        values = setNames(shape_values, shape_levels)
      ) +
      guides(
        shape = guide_legend(override.aes = list(size = 4))
      ) 
  }
  
  p1 <- p1+theme_bw()+
    geom_vline(xintercept=0, color="gray", linetype="dashed")+
    geom_hline(yintercept=0, color="gray", linetype="dashed")+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", linewidth=0.25),
          panel.border = element_rect(colour="black", linewidth=0.25))+
    scale_fill_manual(values=sample.color, guide="none")+
    scale_color_manual(values=sample.color, labels=new_label)+
    labs(x=xlab, y=ylab, title=title, subtitle = substitle)
 
  list(plot=p1, new_label=new_label)
}

pcoa.infection <- zy_pcoa(dt=vir.infection, sample_map=meta.infection, group="Disease", ID="Sample_ID", sample.color=c("#ee9d9e","#9ee4d8"),
                          ado_method="bray", pca_method="bray",shape_by="NCBI_BioProject_ID",
                          levels=0.95, star_plot=F, ellipse_plot=T,
                          cut_rate = NA, cut_num=NA,
                          title="PCoA of respiratory infection", x=1, y=2, ados=T, mydist=dist.vir.infection)
#Figure 6d----------------------------------------------------------------------
crd <- c("COPD","Cystic fibrosis")
meta.crd <- filter(meta.disease,Disease_name %in% crd)
vir.crd <- vir.disease[,meta.crd$Sample_ID]
dist.vir.crd <- dist.vir.disease[meta.crd$Sample_ID,meta.crd$Sample_ID]

meta.crd$group <- paste(meta.crd$Disease_name,meta.crd$Disease,meta.crd$NCBI_BioProject_ID,sep="_")

pcoa.crd <- zy_pcoa(dt=vir.crd, sample_map=meta.crd, group="group", ID="Sample_ID", sample.color=c("#4e728f","#e28c68","#b3a6c6","#d5dc8f"),
                    ado_method="bray", pca_method="bray",shape_by=NULL,
                    levels=0.95, star_plot=F, ellipse_plot=T,
                    cut_rate = NA, cut_num=NA,
                    title="PCoA of chronic respiratory disease", x=1, y=2, ados=T, mydist=dist.vir.crd)
#Figure 6e----------------------------------------------------------------------
host_plot <- host.stat %>%
  mutate(host.genus = ifelse(count == 1, "others", host.genus)) %>%
  group_by(host.genus) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(count) %>%
  mutate(host.genus = factor(host.genus, levels = c("others",setdiff(host.genus,"others"))))  # 保持排序

p <- ggplot(host_plot, aes(x = host.genus, y = count)) +
  geom_col(fill = "#9ee4d8") +
  coord_flip() +
  labs(x = "Host genus", y = "Number of vOTUs") +
  theme_bw(base_size = 14)

#Figure 6f----------------------------------------------------------------------
calc_auc <- function(dt, pred=NA, true=NA, group=NULL,acc=F,levels=NA,
                     boot_n=2000){
  roc.list = list()
  
  if(typeof(group) == "NULL"){
    if(sum(is.na(levels)) == 1){
      levels = unique(dt[,true])
    }
    roc.list['AUC'] = list( roc(dt[,true], dt[,pred], levels=levels))#计算roc曲线并储存为list
    grps = "AUC" 
    
  }else{
    if(length(group) > 1){ 
      old_group = group
      dt = my_join(dt, group)
      group = "zy_tmp_group"
    }
    grps = unique(dt[,group])
    if(length(grps) == 1){ 
      if(sum(is.na(levels)) == 1){
        levels = unique(dt[,true])
      }
      roc.list['AUC'] = list(roc(dt[,true], dt[,pred], levels=levels))
    }else{
      for(g in grps){
        temp_dt = dt[dt[,group]==g,]
        if(sum(is.na(levels)) == 1){
          levels = unique(temp_dt[,true])
        }
        roc.list[as.character(g)] = list(  roc(temp_dt[,true], temp_dt[,pred], levels=levels))
      }
    }
  }
  names_ = names(roc.list)
  result_auc = matrix(NA, ncol=6,nrow=length(grps),
                      dimnames=list(names_, c("low","auc","high","low_acc","acc","high_acc")))
  for(i in names_){
    b = sprintf("%0.4f",ci(roc.list[[i]], of="auc")*100) 
    result_auc[i,c(1,2,3)] = as.numeric(b)
    if(isTRUE(acc)){ 
      ac = ci.coords(roc.list[[i]], x="best", ret="accuracy", transpose=F) 
      ac = sprintf("%0.4f",unlist(ac)*100) 
      ac[2] = sprintf("%0.4f",coords(roc.list[[i]], x="best", ret="accuracy", transpose=F)*100)
      result_auc[i,c(4,5,6)] = ac
    }
  }
  
  result_auc = as.data.frame(result_auc)
  if(exists("old_group")){
    result_auc$zy_tmp_group = rownames(result_auc)
    result_auc = my_split(result_auc, group, old_group)
    result_auc$zy_tmp_group <- NULL
    rownames(result_auc) = NULL
  } else if( typeof(group) != "NULL" ){
    result_auc[,group] = rownames(result_auc)
    rownames(result_auc) = NULL
  }
  list(table=result_auc)
}

get_plot <- function(dt, title,levels_order = NULL){
  x = calc_auc(dt, pred="control", true="group", group=c("train","predict"))
  x$table
  
  pdt <- dcast(x$table, train ~ predict, value.var='auc')
  rownames(pdt) = pdt[,1]
  pdt = pdt[,-1]
  #pdt = pdt[c(colnames(pdt), "LevelOneStudyOut"),]
  
  if(!is.null(levels_order)){
    full_order <- c(levels_order, "LevelOneStudyOut")
    pdt <- pdt[full_order, levels_order]
  } else {
    pdt <- pdt[c(colnames(pdt), "LevelOneStudyOut"),]
  }
  
  heat_matrix = pdt
  
  mycols = colorRamp2(
    breaks = c(0, 50, 60, 70, 80, 90, 100),
    colors = c(
      "white",
      "#deebf7",
      "#c6dbef",
      "#9ecae1",
      "#6baed6",
      "#3182bd",
      "#08519c"
    )
  )
  
  sample_order <- colnames(heat_matrix)
  ann_df <- meta.infection %>%
    select(Project_samplesite, Disease_name) %>%
    distinct() %>%
    slice(match(sample_order, Project_samplesite)) %>%
    column_to_rownames("Project_samplesite")
  ann_df$Disease_name <- gsub("_"," ",ann_df$Disease_name)
  color.infection <- color.disease[unique(ann_df$Disease_name)]
  col_ha <- HeatmapAnnotation(df = ann_df, col = list(Disease_name = color.infection),
                              show_annotation_name = F)
  
  row_df <- rbind(ann_df, "LevelOneStudyOut" = c(Disease_name = "LevelOneStudyOut"))
  color.infection.row <- color.infection
  color.infection.row["LevelOneStudyOut"] <- "white"
  row_ha <- rowAnnotation(df = row_df, col = list(Disease_name = color.infection.row),
                          show_annotation_name = F)
  
  ht <- Heatmap(heat_matrix, 
                column_title = title,
                col = mycols,        
                border=T,
                cluster_rows = F,
                cluster_columns = F,
                row_names_gp = gpar(fontsize=10),column_names_gp = gpar(fontsize=10),
                #rect_gp = gpar(col = "white", lwd = 5), 
                rect_gp = gpar(col = NA),
                row_gap=unit(0,'mm'), column_gap = unit(0,'mm'),
                bottom_annotation = col_ha,
                right_annotation = row_ha,
                cell_fun = function(j, i, x, y, width, height, fill){
                  if(heat_matrix[i,j] > 0){
                    grid.text(sprintf("%0.2f", heat_matrix[i,j]), x, y, gp=gpar(fontsize=10))
                  }
                }
  )
  return(ht)
}

rf.result <- read.table("res.all.prediction.tsv", sep="\t", header=T, check.names=F)
ht <- get_plot(dt=rf.result, title="Random forest using 85 features",
               levels_order=c("BWPR_this_study_Oropharynx","BWPR_this_study_Sputum",
                              "PRJNA413615_Oropharynx","INFPN_unpublished_Sputum",
                              "VP_unpublished_Sputum","PRJNA781460_Nasopharynx"))

#Figure 6g----------------------------------------------------------------------
zy_alpha = function(dt=NA, sample_map=NA, group="Group", ID="Sample", # 必须参数
                    index="shannon", 
                    sample.color=NA, 
                    box_width=0.5, 
                    title="alpha diversity",
                    violin = F
){
  ## colors 
  if (any(is.na(sample.color))){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  message(paste(length(sample.color), "of groups to plot"))
  
  ## align dt and group
  inter <- intersect(colnames(dt),sample_map[,ID])
  if (length(inter) == 0) {  
    stop("Error: No overlapping columns found between dt and sample_map.")  
  } 
  dt <- dt[,colnames(dt) %in% inter]
  dt = dt[rowSums(dt)!=0,]
  sample_map <- filter(sample_map, sample_map[,ID] %in% inter)
  sample_map <- filter(sample_map, !is.na(!!sym(group)))
  
  #alpha
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }else{
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  dm = merge(alpha,sample_map, by.x='row.names', by.y=ID)
  comp = combn(as.character(unique(dm[,group])),2,list)
  
  p = ggplot(dm, aes(x=.data[[group]], y=alpha,fill=.data[[group]]))
  if(isTRUE(violin)){
    p <- p+
      geom_violin()+
      geom_boxplot(width=box_width, fill="white",
                   position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.colour = NA)
  }else{
    p <- p+ 
      geom_boxplot(position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
  }
  
  ylabs = structure(c("Number of OTUs","Shannon index", "1 - Simpson index", "Invsimpson index"),
                    names=c("obs", "shannon", "simpson","invsimpson"))
  ylab = ylabs[tolower(index)]
  
  
  p <- p+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values=sample.color)+
    geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=sigFunc)+
    labs(title=title, y = ylab, x=NULL)
  return(p)
}

op.shannon <- zy_alpha(dt=op.vir, sample_map=op.map, group="Group", ID="Sample", 
                          index="shannon", 
                          sample.color=color.op, 
                          box_width=0.5, 
                          title= "Shannon index of URTI groups (op samples)", 
                          violin = F)

op.obs <- zy_alpha(dt=op.vir, sample_map=op.map, group="Group", ID="Sample", 
                      index="obs", 
                      sample.color=color.op, 
                      box_width=0.5, 
                      title= "Observed vOTUs of URTI groups (op samples)", 
                      violin = F)


#Figure 6h----------------------------------------------------------------------
zy_pcoa_se_sd <- function(dt=NA, sample_map=NA, group=NA, ado_group=NA,ID=NA, sample.color=NULL,
                          ado_method="bray", pca_method="bray",
                          err_bar_type="se",
                          cut_rate = NA, cut_num=NA,
                          err_width=0.01,
                          title="PCoA", x=1, y=2){
  
  fmt_profile = align_dt_sample(dt, sample_map, ID=ID)
  dt = fmt_profile$dt
  sample_map = fmt_profile$sample_map
  
  if(is.finite(cut_rate) || is.finite(cut_num)){
    sample_map = filter_group(sample_map, group=group, cut_rate = cut_rate, cut_num=cut_num)
  }
  
  ## colors 
  if ( typeof(sample.color) == "NULL" ){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(new_label=paste(!!sym(group), " (", count, ")", sep=""))
  
  new_label <- structure(group_summ$new_label,names=as.character(unlist(group_summ[,group])))
  
  message(paste(length(unique(sample_map[,group])), "of groups to plot"))
  
  otu.dist = vegdist(t(dt), method = pca_method) # calc dist matrix
  
  if(length(unique(sample_map[,ado_group])) > 1){
    ## adonis
    ado = adonis2(otu.dist~sample_map[,ado_group], method = ado_method,na.action = na.omit)
    ado_r2 = round(ado$R2[1], digits = 4)
    ado_p = ado$`Pr(>F)`[1]
  }else{
    ado_r2 = NA
    ado_p = NA
  }
  
  ## PCoA
  pcoa = cmdscale(otu.dist, k=10, eig=T)
  eigs = signif(pcoa$eig/sum(pcoa$eig), 4)*100
  point = pcoa$points
  
  colnames(point) = paste("pcoa.", 1:ncol(point),sep="")
  
  xlab = paste("PCoA", x, " (",eigs[x],"%)", sep="")
  ylab = paste("PCoA", y, " (",eigs[y],"%)", sep="")
  
  substitle <- paste0("'R'^2~'='~'", ado_r2, "'~~italic('p')~'='~'", ado_p, "'") %>% 
    as.formula() %>% 
    eval()
  group = ifelse(length(unique(sample_map[,group])) == 1, sample.color[1], group)
  point = point[sample_map[,ID], ]
  dm = merge(point, sample_map, by.x='row.names', by.y=ID)
  tmp_mean = aggregate(point, by=list(c(sample_map[,group])), mean)
  tmp_sd  = aggregate(point, by=list(c(sample_map[,group])), sd)
  tmp_sd  = aggregate(point, by=list(c(sample_map[,group])), sd)
  
  
  plot_data = do.call(data.frame, aggregate(point, by=list(c(sample_map[,group])), 
                                            FUN = function(x){
                                              c(mean = mean(x), sd = sd(x), se = sd(x)/sqrt(length(x)))
                                            }))
  
  plot_x = paste("pcoa.", x, ".mean", sep="")
  plot_y = paste("pcoa.", y, ".mean", sep="")
  
  x_offset = paste("pcoa.", x, ".", err_bar_type, sep="")
  y_offset = paste("pcoa.", y, ".", err_bar_type, sep="")
  
  p1 <- ggplot(data=plot_data, aes(x = .data[[plot_x]], y = .data[[plot_y]], color = Group.1))+
    geom_point(size=10) +
    geom_errorbar(aes( xmin = .data[[plot_x]] - .data[[x_offset]],
                       xmax = .data[[plot_x]] + .data[[x_offset]]),
                  width=err_width
    )+
    geom_errorbar(aes( ymin = .data[[plot_y]] - .data[[y_offset]],
                       ymax = .data[[plot_y]] + .data[[y_offset]]),
                  width=err_width
    )+
    theme_bw()+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", linewidth=0.25),
          panel.border = element_rect(colour="black", linewidth=0.25))+
    scale_fill_manual(values=sample.color, guide="none")+
    scale_color_manual(values=sample.color, labels=new_label)+
    labs(x=xlab, y=ylab, title=title, subtitle = substitle)
  list(plot=p1, new_label=new_label)
}

op.pcoa.sd <- zy_pcoa_se_sd(dt=op.vir, sample_map=op.map, group="Group", ado_group="Group",ID="Sample", sample.color=color.op,
                               ado_method="bray", pca_method="bray",
                               err_bar_type="se",
                               cut_rate = NA, cut_num=NA,
                               err_width=0.01,
                               title="PCoA of URTI groups (op samples)", x=1, y=2)

#Figure 6i----------------------------------------------------------------------
group.color <- c("#80bdc1","#e193b3","#d8ae9c","#99c06e")
group.color <- setNames(group.color,unique(pathogen.meta$Infection.group))

p.votu.shannon <- zy_alpha(dt=dtvir.pathogen, sample_map=pathogen.meta, group="Infection.group", ID="Sample",
                           sample.color=group.color, 
                           index="shannon", 
                           box_width=0.5, 
                           title="Virome shannon", 
                           violin = F
)

p.votu.obs <- zy_alpha(dt=dtvir.pathogen, sample_map=pathogen.meta, group="Infection.group", ID="Sample",
                       sample.color=group.color, 
                       index="obs", 
                       box_width=0.5, 
                       title="Virome obs", 
                       violin = F
)

#Figure 6j----------------------------------------------------------------------
zy_dbrda <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NULL,
                     ado_method="bray", pca_method="bray",
                     levels=0.95,ellipse_plot=F,star_plot=F,
                     title="dbRDA", x=1, y=2){
  fmt_profile = align_dt_sample(dt, sample_map, ID=ID)
  dt = fmt_profile$dt
  sample_map = fmt_profile$sample_map
  
  ## colors 
  
  if (typeof(sample.color) == "NULL"){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(new_label=paste(!!sym(group), " (", count, ")", sep=""))
  new_label <- structure(group_summ$new_label,names=unlist(group_summ[,group]))
  
  message(paste(length(sample.color), "of groups to plot"))
  
  otu.dist = vegdist(t(dt), method = pca_method)
  
  if(length(unique(sample_map[,group])) > 1){
    ## adonis
    ado = adonis2(otu.dist ~ sample_map[,group], method = ado_method)
    ado_r2 = round(ado$R2[1], digits = 4)
    ado_p = ado$`Pr(>F)`[1]
    
  }else{
    ado_r2 = NA
    ado_p = NA
  }
  
  ## dbrda
  db_rda = capscale(otu.dist ~ sample_map[,group])
  db_rda_score = scores(db_rda, choices = 1:10)
  point = db_rda_score$sites
  
  eigs = round(summary(eigenvals(db_rda))[2,]*100, digits=4)
  eigs_name = names(eigs)
  
  xlab = paste(eigs_name[x], " (", eigs[x], "% )", sep="")
  ylab = paste(eigs_name[y], " (", eigs[y], "% )", sep="")
  subtitle <- paste0("'R'^2~'='~'", ado_r2, "'~~italic('p')~'='~'", ado_p, "'") %>% 
    as.formula() %>% 
    eval()
  
  dm = merge(point, sample_map, by.x='row.names', by.y=ID)
  p1 <- ggscatter(data=dm, x=eigs_name[x],y=eigs_name[y],
                  color=group,
                  ellipse.level = levels,ellipse = ellipse_plot,
                  star.plot = star_plot
  )+
    theme_bw()+
    geom_vline(xintercept=0, color="gray", linetype="dashed")+
    geom_hline(yintercept=0, color="gray", linetype="dashed")+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", linewidth=0.25),
          panel.border = element_rect(colour="black", linewidth=0.25))+
    scale_fill_manual(values=sample.color, guide="none")+
    scale_color_manual(values=sample.color, labels=new_label)+
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle)
  list(plot=p1, new_label=new_label)
}

p.infection.db <- zy_dbrda(dt=dtvir.pathogen, sample_map=pathogen.meta, group="Infection.group", ID="Sample",
                           sample.color=group.color,
                           ado_method="bray", pca_method="bray",
                           levels=0.65,ellipse_plot=T,star_plot=F,
                           title="Virome dbRDA", x=1, y=2)
#Figure 6k----------------------------------------------------------------------
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

best.roc <- plot_roc(dt=preddt, pred="prob", true="Group", group="type", levels=c("Infection","Non-infection"),
                     fill=T, title="Infection",
                     cols = c("Viral infection" = "#f7999a", "Bacterial infection" = "#7ecebb"), conf_level=0.95, boot_n=2000)

