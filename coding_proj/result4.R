
#----result 3 #-----------


#---library R  package#--------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)

#---result1 analyses#-------------
result_path <- "result_proj"
dir.create(result_path)
#--creat path save result2
res1path <- paste(result_path,"/result4",sep = "")
dir.create(res1path)


#--inpurt data#-----------


#--micro root
ps <- readRDS("./dat_proj/Micro_seq_data/ps.rds")

map = as.data.frame(sample_data(ps))
map$Group = gsub("CK","Control",map$Group)
sample_data(ps) = map

ps_sub <- phyloseq::subset_samples(ps,Group %in% c("Control_insect","CF_insect","BOF_insect"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
otu_table = as.data.frame(t(vegan_otu(ps_sub)))
#-- change name
colnames(otu_table) = gsub("B","",gsub("_insect_","",colnames(otu_table)))
map = sample_data(ps_sub)
map$ID = gsub("B","",gsub("_insect_","",as.character(map$ID)))
row.names(map) = gsub("B","",gsub("_insect_","",as.character(map$ID)))
map$Group <- gsub("B","",map$Group2)
map$Group = gsub("CK","Control",map$Group)

# rebuilding phylsoeq abject for unified sample name
ps_insect <- phyloseq(otu_table(otu_table,taxa_are_rows = T),
                    tax_table(ps_sub),
                    sample_data(map)
)
ps_insect


#-- set size of picture #-------

length(unique(map$Group))
num <- map$Group %>% 
  unique() %>%
  length()

num
width = num
height = num

#--set theme of plot#-----

#----set color or fill#------------
library(RColorBrewer)#调色板调用包
#调用所有这个包中的调色板
display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
display.brewer.pal(9,"Set1")
# colset1 <- brewer.pal(9,"Set1")[c(6,5,4)]
colset1 <- brewer.pal(9,"Set1")[c(5,6,4)]

colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))



#  order of x axis #------------
axis_order = c("Control","CF","OF")



##--root microbiol community analyse #-----------
subpath <- paste(res1path,"/result_insect_micro",sep = "")
dir.create(subpath)


#--alpha diversity#----------
alppath = paste(subpath,"/alpha/",sep = "")
dir.create(alppath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")

# 使用分面出图，得益于我开发的R包EasyStat
alp = alpha(ps = ps_insect,inde="Shannon",group = "Group",Plot = TRUE )
index= alp[[3]]
head(index)
#--从这里发现bof第六个样本存在问题：BOF_insect_6
sel = c(7,12,9)
data = cbind(data.frame(ID = index$ID,group = index$Group),index[sel])
head(data)
#
result = MuiKwWlx(data = data,num = c(3:5))
result

FileName <- paste(alppath,"/alpha_diversity_all_abc.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(alppath,"/alpha_diversity_Sample.csv", sep = "")
write.csv(index,FileName,sep = "")

result1 = FacetMuiPlotresultBox(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]] + scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1


res = FacetMuiPlotresultBar(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1)
p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) + theme_bw()  + 
  mytheme1+ guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = FacetMuiPlotReBoxBar(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 1)
p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + theme_bw()  + mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(alppath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 5, height = 8)

FileName <- paste(alppath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = 5, height = 8)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = 5, height = 8)


#--beta PCOA ordination and plot#------------
betapath = paste(subpath,"/beta/",sep = "")
dir.create(betapath)

source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

result = BetaDiv(ps = ps_insect, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)
p2_1 = result[[1]] + scale_fill_manual(values = colset1,guide = F)+
  scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.2,0.2))

p2_1
#带标签图形出图
p2_2 = result[[3]] + scale_fill_manual(values = colset1,guide = F)+
  scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.2,0.2))

p2_2
#提取总体比较
TResult =result[[5]]
head(TResult)

#---------精修图
plotdata =result[[2]]
head(plotdata)
# 求均值
cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)

# p2$layers[[2]] = NULL
# library(ggcor)
library(ggsci)
p2_3 = p2_1 +geom_segment(data = segs,
                          mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
  scale_fill_manual(values = colset1,guide = F)+
  scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.2,0.2))


p2_3
# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_adonis.csv", sep = "")
write.csv(pair,FileName)

FileName <- paste(betapath,"/a2_bray_PCOA.pdf", sep = "")
ggsave(FileName, p2_1, width = 4, height = 4)
FileName1 <- paste(betapath,"/a2_bray_PCOA.jpg", sep = "")
ggsave(FileName1 , p2_1, width = 4, height = 4)

FileName <- paste(betapath,"/a2_bray_PCOA_label.pdf", sep = "")
ggsave(FileName, p2_2, width = 8, height = 8)
FileName1 <- paste(betapath,"/a2_bray_PCOA_label.jpg", sep = "")
ggsave(FileName1 , p2_2, width = 8, height = 8)

FileName <- paste(betapath,"/a2_bray_PCOA_adjust.pdf", sep = "")
ggsave(FileName, p2_3, width = 6, height = 6)
FileName1 <- paste(betapath,"/a2_bray_PCOA_adjust.jpg", sep = "")
ggsave(FileName1 , p2_3, width = 6, height = 6)

#------phylum starck bar plot #----------


# BiocManager::install("ggalluvial")
# devtools::install_github("taowenmicro/EasyStat",force = TRUE)
library(EasyMicrobiome)

#-------门类水平展示
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(subpath,"/barplot/",sep = "")
dir.create(barpath)
j = "Phylum"
result = barMainplot(ps = ps_insect,j = "Phylum",rep = 6,axis_ord = NULL,label = FALSE ,sd = FALSE,Top = 10)

p3_1 <- result[[1]] + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p3_1
p3_2  <- result[[3]] + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p3_2



FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
ggsave(FileName1, p3_1, width = 6, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
ggsave(FileName2, p3_1, width = 6, height =8 )

FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
ggsave(FileName1, p3_2, width = 6, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
ggsave(FileName2, p3_2, width = 6, height =8 )

#---flower plot #-------

flowpath = paste(subpath,"/flowplot/",sep = "")
dir.create(flowpath)
# source("G:\\Shared_Folder\\Function_local\\R_function\\micro/flowerplot.R")
# flowerplot(ps = ps_insect,rep = 6,path =flowpath )
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")
# #---调整花瓣效果
p0_2 <- ggflower(ps = ps_insect,
                 rep = 3,
                 group = "Group",
                 start = 1, # 风车效果
                 m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                 a = 0.2, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "green"
)
p0_2 = p0_2 +  scale_fill_manual(values = colset1,guide = F)
p0_2
FileName1 <- paste(flowpath,"/flower",".pdf", sep = "")
ggsave(FileName1, p0_2, width = 6, height =6 )
FileName2 <- paste(barpath,"/flower",".jpg", sep = "")
ggsave(FileName2, p0_2, width = 6, height =6 )


#-- network #-----------

netpath = paste(subpath,"/network/",sep = "")
dir.create(netpath)

library(igraph)
library(sna)

result = ggClusterNet::network(ps = ps_insect,N = 0.0001,r.threshold=0.6,p.threshold=0.05,label = FALSE,path = netpath ,zipi = F)
# ?network
# 全部样本的网络比对
p4_1 = result[[1]] + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p4_1
# 全部样本网络参数比对
data = result[[2]]
data
plotname1 = paste(netpath,"/network_all.jpg",sep = "")
ggsave(plotname1, p4_1,width = 10,height = 8)
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p4_1,width = 10,height = 8)

tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)

#---conbine plot of result1 #-----------

library(patchwork)

layout <- "
AAAABBBBBCCCC
AAAABBBBBCCCC
AAAADDDDDCCCC
AAAADDDDDCCCC
EEEEEEEEEEEEE
EEEEEEEEEEEEE
EEEEEEEEEEEEE
"
p_result4_1 <- p1_1 + p2_3 + p3_2 + p0_2 + p4_1 +
  plot_layout(design = layout)

# sp1 <- p1_1 |(p2_3/p0_2)|p3_1
# sp1/p4_1

plotname = paste(res1path,"/result4_1.pdf",sep = "")
ggsave(plotname, p_result4_1,width = 14,height = 12)
# p_result1_1 <- p1_1 |(p2_3)|p3_1|p4_1
# p_result1_1 
# plotname = paste(subpath,"/result1_1.pdf",sep = "")
# ggsave(plotname, p_result1_1,width = 30,height = 10)


#----pls-pm--model find the the most important object#----------

# 代码写好了，等待后面在一起整合，吧，因为这部分想要用R出图，似乎没有什么比较好的方案


#-source tracker 用于寻找微生物群落来源

#------------虫子微生物群落来源#-----------

#------------虫子微生物群落来源#-----------
stpath = paste(subpath,"/source_tracker/",sep = "")
dir.create(stpath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro\\FEAST.R",encoding = "UTF-8")
# 取子集
ps = readRDS("./Seq_cab_66/data/ps.rds")

map = as.data.frame(sample_data(ps))
map
ps_sub <- phyloseq::subset_samples(ps, Group2 %in% c("BOF"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 200, TRUE);ps_sub #筛选序列数量大于1的
# ps_sub = transform_sample_counts(ps_sub, function(x) x / sum(x) )

# #清空内存
# rm(list=ls()) 



result = FEAST(ps = ps_sub,group = "Group",sinkG = "BOF_insect",sourceG = c("BOF_soil","BOF_root","BOF_leaf"))
result


p5_1 <- Plot_FEAST(data = result)
p5_1
p5_2 <- MuiPlot_FEAST(data = result)
p5_2
ps_sub <- phyloseq::subset_samples(ps, Group2 %in% c("CF"));ps_sub
# ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 200, TRUE);ps_sub #筛选序列数量大于1的

result = FEAST(ps = ps_sub,group = "Group",sinkG = "CF_insect",sourceG = c("CF_soil","CF_root","CF_leaf"))
result

p5_2 <- Plot_FEAST(data = result)
p5_2

MuiPlot_FEAST(data = result)

# 取子集
ps = readRDS("./Seq_cab_66/data/ps.rds")

ps_sub <- phyloseq::subset_samples(ps, Group2 %in% c("CK"));ps_sub
# ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 200, TRUE);ps_sub #筛选序列数量大于1的

result = FEAST(ps = ps_sub,group = "Group",sinkG = "CK_insect",sourceG = c("CK_root","CK_leaf"))
result

# -例子
p5_3 <- Plot_FEAST(data = result)
p5_3
MuiPlot_FEAST(data = result)

# MuiPlot_FEAST(data = result)

library(patchwork)
p5 <- p5_1|p5_2|p5_3
p5
plotname1 = paste(stpath,"/source_tracker.jpg",sep = "")
ggsave(plotname1, p5,width = 10,height = 4)
plotname1 = paste(stpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p5,width = 10,height = 4)


#------- R index--from soil--to --insect#-------

ind = c(0.421,0.275,0.819,0.221)

dat = data.frame(ID = c("soil","root","leaf","insect"),Adonis_text_R2 = ind)
dat

subpath = paste(res1path,"/R-pcoa/",sep = "")
dir.create(subpath )

p <- ggplot(dat) + geom_bar(aes(x = ID,y = Adonis_text_R2,fill= ID),stat = "identity") 
p
p = p +  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset2) + scale_x_discrete(limits = c("soil","root","leaf","insect")) + mytheme1
p
plotname = paste(subpath,"/R0-soil-insect.pdf",sep = "")
ggsave(plotname, p,width = 6,height = 4)



#---distance#-----------

dispath = paste(res1path,"/distance/",sep = "")
dir.create(dispath)

map = as.data.frame(sample_data(ps))
head(map)


tre1 = c("BOF_soil","BOF_root","BOF_leaf","BOF_insect")
tre2 = c("CF_soil","CF_root","CF_leaf","CF_insect")

# tre3 = c("CK","C","AC","CD","CD","ABC","BCD","ABCD")
# tre4 = c("CK","D","AD","BD","CD","ACD","BCD","ABCD")

# datab  = data.frame(tre1 = tre1,tre2 = tre2,tre= tre3,tre3 = tre4)
datab  = data.frame(tre1 = tre1,tre2 = tre2)

# refer = "BOF_soil"
refer = c("BOF_soil","CF_soil")

plot = list()
i = 1 
for (i in 1:length(refer)) {
  tre = as.character(datab[[i]])
  map = as.data.frame(sample_data(ps))
  maps<- dplyr::filter(as.tibble(map),Group %in% tre)
  row.names(maps) = maps$ID
  ps_sub = ps
  sample_data( ps_sub ) = as.data.frame(maps);ps_sub 
  
  
  ps_rela  = transform_sample_counts(ps_sub , function(x) x / sum(x) );ps_rela 
  
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  otu= as.data.frame(t(vegan_otu(ps_rela)))
  
  head(otu)
  #-添加第一个时间的均值作为标准，求取均值列加入其后
  mapsub = as.data.frame(sample_data(ps_rela))
  #挑选time为1的样本
  ID = mapsub$ID[mapsub$Group == refer[i]]
  #挑选样本
  head(otu)
  ID = as.character(ID)
  index = otu[,ID] # 筛选并计算均值
  head(index,n = 5)
  
  index_mean= rowMeans(index)
  otu$index_mean=index_mean # 均值添加到OTU表
  head(otu)
  # 计算包括终点均值的所有样品bray距离
  bray_curtis = vegan::vegdist(t(otu), method = "bray")
  bray_curtis= as.matrix(bray_curtis)
  
  #--------------整理数据-
  # 计算各组内部差异程度
  # 建立一个存储组内差异的数据框
  dat = t(as.data.frame(c("sampleA","sampleB","0.15","group","genosite")))
  colnames(dat) = c("sampleA","sampleB","distance","group","type")
  rownames(dat) = c("test")
  
  
  # 每个样品与final对应的距离
  # ii = 1
  for (ii in tre){
    group = rownames(mapsub[mapsub$Group %in% ii,])
    # m = 1
    for (m in 1:(length(group))) {
      x = c(group[m],"index_mean",bray_curtis[group[m],"index_mean"],ii,tre)
      dat=rbind(dat,x)
    }
  }
  dat = as.data.frame(dat[-1,], stringsAsFactors=F) # 删除首行框架数据
  
  # dat = dat[dat$distance != 0,] # 
  # 距离转换为3位数值
  dat$distance=round(as.numeric(as.character(dat$distance)), digits=3)
  # 分组添加levels，分组有时间顺序
  dat$group = factor(dat$group, levels=unique(dat$group))
  
  all = dat
  # 开始尝试散点图拟合
  require(splines)
  require(MASS)
  library(ggplot2)
  
  all1 = all
  all1$group = as.numeric(all1$group)
  
  head(all1)
  library(dplyr)
  iris_groups<- group_by(all, group)
  all3<- dplyr::summarise(iris_groups, mean(distance))
  colnames(all3) = c("group","mean")
  
  datawt = data.frame(ID  =all$sampleA,group = all$group,distance = all$distance)
  library(EasyStat)
  result= aovMcomper (data = datawt, i= 3,method_Mc = "Tukey")
  # 提取多重比较结果
  result[[1]]
  
  PlotresultBox = aovMuiBoxP(data = datawt, i= 3,sig_show ="abc",result = result[[1]])
  #提取图片
  p = PlotresultBox[[1]]
  p = p + geom_line(aes(x=group, y=mean,group = 1),data = all3) +geom_point(aes(x=group, y=mean),all3,alpha=1,pch = 21,fill="red") +
    labs(y = "Bray curies distance")
  p  = p + theme_bw() + mytheme1 
  
  plot[[i]] = p
  
}

library(ggpubr)
p  = ggarrange(plotlist = plot, common.legend = TRUE, legend="right")
p

FileName <- paste(dispath,"/distance.pdf", sep = "")
ggsave(FileName, p, width = 12, height = 6)
