setwd('/share/data2/guorc/Project/AVGC/2025.Nov7.AVGC/data/data_generation/03.fun')

#抗菌基因的收集Endolysin、AMP etc================================================================================================================
pc = read.csv('votu.prot.clu.tsv',sep='\t',stringsAsFactors = F,header = F)
# pcf = pc[pc$cluster!='',]
# ff=max(as.numeric(sub('PC_','',pcf$cluster)))
# pcf = pc[pc$cluster=='',]
# pc$cluster[pc$cluster=='']=paste('PC_',(ff+1):(ff+1+nrow(pcf)),sep='')
pc$cluster = paste('PC_',as.numeric(as.factor(pc$V1)),sep='')


#refseq##
vref = read.csv('votu.prot.clu.faa.ncbi.anno',sep='\t',stringsAsFactors = F,header = F)
aa = vref[grepl('endolysin|lysozyme|peptidoglycan hydrolase|peptidoglycanhydrolase',vref$V13,ignore.case = T),]
#aa$pc = pc$cluster[match(aa$V1,pc$V2)]
ff = data.frame(id=unique(aa$V1),db='refseq')


#pfam##
pfam = read.csv('votu.prot.clu.faa.pfam',sep='',header=F,stringsAsFactors = F)
aa = pfam[grepl('Amidase_2|Amidase_3|Amidase_5|lysozyme',pfam$V7,ignore.case = T),]
#aa$pc = pc$cluster[match(aa$V1,pc$V2)]
ff = rbind(ff,data.frame(id=unique(aa$V1),db='pfam'))


#kegg##
kegg = read.csv('votu.kegg.btp.uniq',sep='',header=F,stringsAsFactors = F)
kegg$V2 = sub('.*\\|','',kegg$V2)
gg = c('K01185','K02395','K03791','K07273','K08259','K11331','K13381','K13915','K17733','K21471',
       'K21473','K21474','K21508','K21509','K25438','K25707','K26029','K26031','K26440','K26837',
       "K18920","K18921","K18919","K18922","K19048","K18862","K19148","K19149","K19150","K19151",
       "K19152","K19153","K19154","K19779","K07171","K18841","K19155","K06218","K19157","K19158",
       "K07154","K07339","K19160","K05839","K03973","K19163","K19092","K19093","K07341","K19166",
       "K07334","K13651","K18828","K19686","K16214","K21487","K21489","K21491","K21493","K07062",
       "K21496","K19167","K18837","K19168","K18839")
aa = kegg[kegg$V2%in%gg,]
#aa$pc = pc$cluster[match(aa$V1,pc$V2)]
ff = rbind(ff,data.frame(id=unique(aa$V1),db='kegg'))

# #AMP
# aa = amp
# aa$pc = plen$pc[match(aa$V1,plen$id)]
# aa = aa[!is.na(aa$pc),]
# #ff = unique(c(ff,aa$V1))
# ff = rbind(ff,data.frame(id=unique(aa$V1),db='AMP'))

#cazy
cazy = read.csv('votu.cazy.btp.uniq',sep='',header=F,stringsAsFactors = F)
aa = cazy[grepl(paste(c('GH73','GH19','GH23','GH108','GH25','GH24','GH104'),collapse='|'),cazy$V2),]
ff = rbind(ff,data.frame(id=unique(aa$V1),db='cazy'))

endo_db = plyr::count(ff$id) #有几个数据库注释成功
endo = pc[pc$cluster%in%pc$cluster[pc$V2%in%ff$id],]
nrow(endo)
length(unique(endo$cluster))

###绘图展示
host = read.csv('votu.host.taxonomy',sep='\t',stringsAsFactors = F,header = F)
host = unique(host[,c(1,4:8)])
colnames(host) = c("id","phylum","class","order","family","genus")
# dat$pc = endo$cluster[match(dat$V1,endo$protein_id)]
# dat$contigid = sub('_\\d+$','',dat$V2)
# dat$vc = vc$preVC[match(dat$contigid,vc$Genome)]
# dat$vc[is.na(dat$vc)] = 'Others'
# dat$ps = dat$V1
dat = endo
dat$contigid = sub('_\\d+$','',dat$V2)
dat$ps = dat$cluster

# for(x in c("phylum","class","order","family","genus")){
#   aa = unique(data.frame(id=host$id,host=host[,x]))
#   datf = merge(dat,aa,by.x='contigid',by.y='id',all.x = T)
#   datf$host[is.na(datf$host)] = 'Unassigned'
#   bb = plyr::count(unique(datf[datf$host!='Unassigned',c('ps','host')])$ps)
#   print(sum(bb$freq>1)/length(unique(datf$ps)))
#   print(sum(bb$freq==1)/length(unique(datf$ps)))
# }

#宿主的树文件
datf = merge(dat,host,by.x='contigid',by.y='id',all.x = T)
datf[is.na(datf)] = 'Unassigned'
datf = unique(datf)
datf$order = sub('_[A-Z]$','',datf$order)
datf$class[datf$class=='c__Actinomycetes'] = 'c__Actinomycetia'
aa = unique(datf[,c("phylum","class","order","family","genus")])
aa$root = 'root'
aa = aa[,c('root',"phylum","class","order","family","genus")]
library(taxonomizr)
nwk <- makeNewick(taxa = aa, quote = NULL,terminator = ";")
#write.table(nwk,'endo.tree.linage.list.pc.nwk',sep='\t',quote = F,col.names = F,row.names = F)

xx=unique(datf[,5:10])
x1=xx[xx$phylum!='Unassigned',]
x2=xx[!xx$ps%in%x1$ps,]
xx=rbind(x1,x2)
#write.table(xx,'xx',se='\t',row.names = F,quote = F)


#树注释文件
cc = c()
for(x in c("phylum","class","order","family","genus")){
  lab=x
  aa = unique(datf[,c('ps',lab)])
  colnames(aa) = c('ps','host')
  bb = aa[aa$host!='Unassigned',]
  bb = plyr::count(bb$ps)
  aa = aa[!(aa$host=='Unassigned' & aa$ps%in%bb$x),]
  aa$type[aa$ps%in%c(bb$x[bb$freq==1])] = 'uniq'
  aa$type[aa$ps%in%c(bb$x[bb$freq!=1])] = 'share'
  aa$type[is.na(aa$type)] = 'uniq'
  xx = plyr::count(unique(aa[aa$host!="Unassigned",c('ps','type')])$type)
  xx$rate = xx$freq/sum(xx$freq)
  print(x)
  print(xx)
  aa = plyr::count(aa[,c('host','type')])
  aa = dcast(host~type,data=aa)
  aa[is.na(aa)] = 0
  aa$lab=x
  cc = rbind(cc,aa)
}

cc = rbind(cc[cc$lab=='genus',],cc[cc$lab!='genus' & cc$host!='Unassigned',])
anno = paste(cc$host,ifelse(cc$lab=='genus','1','1'),sqrt(cc$share+cc$uniq),cc$share,cc$uniq,sep=',')
#write.table(anno,'endo.tree.linage.list.pc.nwk.lab',sep='\t',quote = F,col.names = F,row.names = F)
cc[cc$lab=='genus','host']
unique(str_sub(cc$host,1,3))



#挑选每个ps中的代表序列，优选选择宿主类型最少的那一条作为代表序列
tab = c()
for(x in unique(datf$pc)){
  aa = datf[datf$pc==x,]
  aa$new_id = paste(sub('PC_0*','s',aa$pc),'_',as.numeric(as.factor(aa$ps)),sep='')
  aa$species = sub('_[A-Z]_','_',aa$species)
  aa$species = sub('_[A-Z]$','',aa$species)
  aa = unique(aa[,c('new_id','V2','ps','species')])
  
  #挑选代表序列
  bb = plyr::count(aa[,c('new_id','ps','V2')])
  cc = bb[bb$ps==bb$V2 & bb$freq==1,]
  bb = bb[!bb$new_id%in%cc$new_id,]
  bb = bb[order(bb$new_id,bb$freq),]
  cc = rbind(cc,bb[!duplicated(bb$new_id),])
  
  #标记可能宿主
  bb = aa[aa$species!='Unassigned',]
  if(nrow(bb)>0){
    bb = aggregate(bb$species,list(bb$new_id),function(x){paste(sort(unique(x)),collapse = ',')})
    cc$host = bb$x[match(cc$new_id,bb$Group.1)]
    cc$host[is.na(cc$host)] = 'Unassigned'
  }else{
    cc$host = 'Unassigned'
  }
  
  
  #统计基因总数
  bb = plyr::count(unique(aa[,c('new_id','V2')])$new_id)
  cc$geneNum = bb$freq[match(cc$new_id,bb$x)]
  
  #序列直接比对成功的数据库数量
  cc$map_db = endo_db$freq[match(cc$V2,endo_db$x)]
  cc$map_db[is.na(cc$map_db)] = 0
  
  tab = rbind(tab,data.frame(new_id=cc$new_id,geneNum=cc$geneNum,gene_id=cc$V2,map_db=cc$map_db,Host=cc$host))
}
