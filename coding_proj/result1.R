


# 备忘： 数据是一个有机的组合，代码也是一个有机的组合，proj是一个总览

cod_path <- "coding_proj"
dir.create(cod_path)

dat_path <- "dat_proj"
dir.create(dat_path)

result_path <- "result_proj"
dir.create(result_path)

#---library R  package#--------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)

#---result1 analyses#-------------

#--creat path save result1
res1path <- paste(result_path,"/result1",sep = "")
dir.create(res1path)




#--inpurt data#-----------

#==GC
psG <- readRDS("./dat_proj/GC_data/GC_soil/ps_soil.rds")
psG
re = otu_table(psG)
#-保留小数点位数达到取整目的
for (i in 1:dim(re)[2]) {
  re[,i] = round(re[,i],0)
}
#--将空缺值填充为0
re[is.na(re)] <- 0
otu_table(psG) <- re

psG  = transform_sample_counts(psG, function(x) x / sum(x) )
gctab_soil = as.data.frame(vegan_otu(psG))
map = sample_data(psG)
map$ID = gsub("B","",as.character(map$sample_ID))

row.names(map) = gsub("B","",as.character(map$sample_ID))
#--change sample name
rownames(gctab_soil) = gsub("B","",as.character(map$sample_ID))
map$sample_ID = NULL
map$Group = as.factor(gsub("B","",as.character(map$Group)))
# rebuilding phylsoeq abject for unified sample name
library(ggClusterNet)
taxG = as.data.frame(vegan_tax(psG))

for (i in 1:length(taxG$group)) {
  if (taxG$group[i] == "short chain carbon organic acids") {
    taxG$group[i] = "small-molecule organic acids"
  }
}
table(taxG$group)

tax_table(psG) = as.matrix(taxG)

psG <- phyloseq(otu_table(gctab_soil,taxa_are_rows = F),
                tax_table(psG),
                sample_data(map)
                )
#--micro

ps <- readRDS("./dat_proj/Micro_seq_data/ps.rds")

map = as.data.frame(sample_data(ps))
map$Group = gsub("CK","Contral",map$Group)
sample_data(ps) = map

ps_sub <- phyloseq::subset_samples(ps,Group %in% c("CF_soil","BOF_soil"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
otu_table = as.data.frame(t(vegan_otu(ps_sub)))
#-- change name
colnames(otu_table) = gsub("B","",gsub("_soil_","",colnames(otu_table)))
# rebuilding phylsoeq abject for unified sample name
ps <- phyloseq(otu_table(otu_table,taxa_are_rows = T),
                tax_table(ps),
                sample_data(psG)
)

#-soil chemical propertity
env = read.csv("./dat_proj/soil_chem/soil_chem_CF_OF.csv")


#-- set size of picture #-------

length(unique(map$Group))
num <- map$Group %>% 
  unique() %>%
  length()

num
width = num
height = num



#----set color or fill#------------
library(RColorBrewer)#调色板调用包
#调用所有这个包中的调色板
display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
display.brewer.pal(9,"Set1")
colset1 <- brewer.pal(9,"Set1")[c(5,4)]
colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))



#  order of x axis #------------
axis_order = c("CF","OF")

#--alpha--diversity#----------
alppath = paste(res1path,"/alpha/",sep = "")
dir.create(alppath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")

# 使用分面出图，得益于我开发的R包EasyStat
alp = alpha(ps = ps,inde="Shannon",group = "Group",Plot = TRUE )
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
  scale_fill_manual(values = colset1) + guides(fill=F)
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
betapath = paste(res1path,"/beta/",sep = "")
dir.create(betapath)

source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

result = BetaDiv(ps = ps, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)
p2_1 = result[[1]] + scale_fill_manual(values = colset1,guide = F)+
  scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.2,0.2))

p2_1
#带标签图形出图
p2_2 = result[[3]] + scale_fill_manual(values = colset1,guide = F)+
  scale_color_manual(values = colset1,guide = F) + may_theme1 + theme(legend.position = c(0.2,0.2))

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
ggsave(FileName, p2_1, width = 6, height = 6)
FileName1 <- paste(betapath,"/a2_bray_PCOA.jpg", sep = "")
ggsave(FileName1 , p2_1, width = 6, height = 6)

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
barpath = paste(res1path,"/barplot/",sep = "")
dir.create(barpath)
j = "Phylum"
result = barMainplot(ps = ps,j = "Phylum",rep = 6,axis_ord = NULL,label = FALSE ,sd = FALSE,Top = 10)

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

flowpath = paste(res1path,"/flowplot/",sep = "")
dir.create(flowpath)
# source("G:\\Shared_Folder\\Function_local\\R_function\\micro/flowerplot.R")
# flowerplot(ps = ps,rep = 6,path =flowpath )
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")
p0_1 <- ggflower(ps = ps,
                 rep = 6,
                 group = "Group",
                 start = 1, # 风车效果
                 m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                 a = 0.2, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow"
)
# p + scale_fill_brewer(palette = "Paired")



FileName1 <- paste(flowpath,"ggflower.pdf", sep = "")
ggsave(FileName1, p0_1, width = 4, height = 4)
FileName2 <- paste(flowpath,"ggflower.jpg", sep = "")
ggsave(FileName2, p0_1, width = 4, height = 4 )

#-- network #-----------

netpath = paste(res1path,"/network/",sep = "")
dir.create(netpath)

library(igraph)
library(sna)
# ,yourmem = mytheme1
source("G:\\Shared_Folder\\Function_local\\R_function\\my_R_packages\\ggClusterNet\\R\\networkplot.R")
result = ggClusterNet::network(ps = ps,N = 0.001,r.threshold=0.6,
                               p.threshold=0.05,label = FALSE,path = netpath ,
                               zipi = TRUE
                               )

# 全部样本的网络比对
p4_1 = result[[1]] + scale_fill_brewer(palette = "Paired") +  scale_x_discrete(limits = axis_order) + mytheme1
p4_1
# 全部样本网络参数比对
data = result[[2]]
plotname1 = paste(netpath,"/network_all.jpg",sep = "")
ggsave(plotname1, p4_1,width = 10,height = 8)
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p4_1,width = 10,height = 8)

tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)

#---conbine plot of result1 #-----------

library(patchwork)
# 尝试了集中布局，不太如人意

layout <- "
AAAABBBBBCCCC
AAAABBBBBCCCC
AAAADDDDDCCCC
AAAADDDDDCCCC
EEEEEEEEEEEEE
EEEEEEEEEEEEE
EEEEEEEEEEEEE
"
p_result1_1 <- p1_1 + p2_3 + p3_2 + p0_1 + p4_1 +
  plot_layout(design = layout)

# sp1 <- p1_1 |(p2_3/p0_2)|p3_1
# sp1/p4_1

plotname = paste(res1path,"/result1_1.pdf",sep = "")
ggsave(plotname, p_result1_1,width = 14,height = 12)

#--------result1  soil GC analyses plot#------------
#-------GC strack bar #-----------
# source("G:\\Shared_Folder\\Function_local\\R_function\\my_R_packages\\EasyMicrobiome\\R\\barMainplot.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")


barpath = paste(res1path,"/GC_soil_tax_class_barplot/",sep = "")
dir.create(barpath)


psG1 <- psG %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    !Metabolite.name %in% c("Unknown","unknown"),
    # row.names(tax_table(ps_rela ))%in%c("SH010924.07FU_KF986690_reps_singleton","SH020983.07FU_JN235282_refs")
  )
psG1

tax = as.data.frame(vegan_tax(psG1))
table(tax$group)


ps_rela  = transform_sample_counts(psG1, function(x) x / sum(x) )


library(EasyMicrobiome)
tax_glom_wt <- function(ps = ps,ranks = "Phylum") {
  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))
  # building group
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)
  taxcon <- tax[1:match(ranks,colnames(tax))]
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]
  row.names(taxcon) <- taxcon[[ranks]]
  pscon <- phyloseq(
    otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    tax_table(as.matrix(taxcon)),
    sample_data(ps)
  )
  return(pscon)
}
psdata <- tax_glom_wt(ps = ps_rela ,ranks = "group")

otu = as.data.frame((vegan_otu(psdata)))
map = as.data.frame(sample_data(psdata))
map
mapp = map[,c(4,2)]
colnames(mapp) = c("ID","group")
data = cbind(mapp,otu)
result = MuiKwWlx(data = data,num = c(3:ncol(data)))

data = data
i =c(3:ncol(data)) 
result = result

res <- MuiPlotStackBar(data = data,i =c(3:ncol(data)) ,result = result,errbar = F)
#--提取图片
g1_1 = res[[1]] + scale_fill_manual(values = colset3) +  scale_x_discrete(limits = axis_order) + mytheme1 + labs (y = "Relative abundance (%)")
  
g1_1
#--提取数据
plotdata = res[[2]]
j = "GCsoil"

FileName1 <- paste(barpath,"/a2_",j,"_wlx_sbar",".pdf", sep = "")
ggsave(FileName1, g1_1, width = 8, height =8 )
FileName2 <- paste(barpath,"/a2_",j,"_wlx_sbar",".jpg", sep = "")
ggsave(FileName2, g1_1, width = 8, height =8 )

FileName2 <- paste(barpath,"/a2_",j,"_wlx_sbar_data",".csv", sep = "")
write.csv( plotdata,FileName2 )

result1 = FacetMuiPlotresultBox(data = data,num = c(3:ncol(data)),result = result,sig_show ="abc",ncol = 3 )
g1_2 <- result1[[1]] + scale_fill_manual(values = colset1) +  scale_x_discrete(limits = axis_order) + mytheme1
g1_2


FileName1 <- paste(barpath,"/a2_",j,"_wlx_sbar_facet",".pdf", sep = "")
ggsave(FileName1, g1_2, width = 10, height =10 )
FileName2 <- paste(barpath,"/a2_",j,"_wlx_sbar_facet",".jpg", sep = "")
ggsave(FileName2, g1_2, width = 10, height =10)

#--------beta ordanation plot #--------
betapath = paste(res1path,"/GCbeta/",sep = "")
dir.create(betapath)

source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

result = BetaDiv(ps = psG1, group = "Group", dist = "bray", method = "PCA", Micromet = "adonis", pvalue.cutoff = 0.05)

#不带标签默认出图
g2_1 = result[[1]]  + scale_fill_manual(values = colset1,guide = NULL) +  scale_x_discrete(limits = axis_order) + mytheme1 +
  theme(legend.position = c(0.6,0.2))
g2_1
#带标签图形出图
g2_2 = result[[3]] + scale_fill_manual(values = colset1,guide = NULL) +  scale_x_discrete(limits = axis_order) + mytheme1 +
  theme(legend.position = c(0.6,0.2))
g2_2


#--提取出图数据
plotdata =result[[2]]
head(plotdata)
#--提取两两比较结果
pairResult =result[[4]]
head(pairResult)

#提取总体比较
TResult =result[[5]]
head(TResult)



FileName <- paste(betapath,"/PCA_plotdata.csv", sep = "")
write.csv(plotdata,FileName)

FileName <- paste(betapath,"/pairResult.csv", sep = "")
write.csv(pairResult,FileName)

FileName <- paste(betapath,"/TResult.csv", sep = "")
write.csv(TResult,FileName)

FileName <- paste(betapath,"/a2_PCA.pdf", sep = "")
ggsave(FileName, g2_1, width = 6, height = 6)


FileName1 <- paste(betapath,"/a2_PCA_label.pdf", sep = "")
ggsave(FileName1 , g2_2, width = 6, height = 6)


plotdata =result[[2]]
head(plotdata)


# 求均值
cent <- aggregate(cbind(x,y) ~ Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)


g2_1$layers[[2]] = NULL
# library(ggcor)
library(ggsci)

g2_3 = g2_1 +geom_segment(data = segs,
                     mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group)) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 20,color = "blue",fill = "blue") +
  scale_fill_manual(values = colset1,guide = F)+
  scale_color_manual(values = colset1,guide = F) + mytheme1 + theme(legend.position = c(0.4,0.2))
g2_3



# 提取两两检测结果
pair = result[[4]]
FileName <- paste(betapath,"Pair_adonis.csv", sep = "")
write.csv(pair,FileName)


FileName <- paste(betapath,"/a2_bray_PCA_adjust.pdf", sep = "")
ggsave(FileName, g2_3, width = 6, height = 6)
FileName1 <- paste(betapath,"/a2_bray_PCA_adjust.jpg", sep = "")
ggsave(FileName1 , g2_3, width = 6, height = 6)


#----------GC different soil heatmap？ #------------


#-----人工指定分组信息
group1 = c("CF","OF")
# group2 = c("CK","CF")
b= data.frame(group1)

source("G:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")

# path = getwd()
wlxpath = paste(res1path,"/GCwlxtest/",sep = "")
dir.create(wlxpath)
ps_rela  = transform_sample_counts(psG1, function(x) x / sum(x) )
result = statSuper(ps = ps_rela,group  = "Group",artGroup = b,method = "wilcox")
head(result)



fileName =  paste(wlxpath,"/wlxtest_arti.csv",sep = "")
write.csv(result,fileName)

#----筛选差异大的分泌物可视化
#----筛选差异大的分泌物可视化
head(result)
difres <- result %>%  filter(CF_OF_fdr < 0.05 & Metabolite.name != "Unknown"& Metabolite.name != "unknown" )
colnames(difres)
difres$Metabolite.name
fileName=  paste(wlxpath,"/wlxtest_arti_fdr05RmUnknown.csv",sep = "")
write.csv(difres,fileName)

#---差异分泌物可视化--这里的是标准化到100的丰度，注意不是到1.
#---
data <- difres[,c(sample_names(psG1))]
head(data)
data$ID = paste(difres$id,difres$Metabolite.name,sep = "_")

row.names(data) = data$ID
data$ID = NULL

# data <-sqrt (data)
data[data > 0.6]<-0.6
# wt2<-sqrt(wt2)


#,fontsize_col = 10 修改行大小
#fontsize_row =  10 修改列自己大小
library(pheatmap)
map = as.data.frame(sample_data(ps))
map
# annotation_row = data.frame(map$Group)
# rownames(annotation_row) = row.names(map)



annotation_col = data.frame(map$Group)
rownames(annotation_col) =  row.names(map)



g3_1 = pheatmap(data,fontsize=6,cellwidth = 20, cellheight = 10,cluster_rows = TRUE,
             color = colorRampPalette(brewer.pal(11,"Spectral"))(60),
             display_numbers = F,fontsize_col = 10,fontsize_row =  10,
             annotation_col = annotation_col)
g3_1

fileName =  paste(wlxpath,"/pheartmap_diff.pdf",sep = "")
ggsave(fileName, g3_1, width = 12, height =15 )

#--气泡图

library(reshape2)
data <- difres[,c(sample_names(ps))]
head(data)
data$ID = paste(difres$id,difres$Metabolite.name,sep = "_")

row.names(data) = data$ID
data$ID = NULL
head(data)
data$id = row.names(data)
pcm = melt(data, id = c("id"))
head(pcm)

# pcm$ID <- factor(pcm$ID,levels=unique(pcm$ID))

colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
#----样本在y轴上
g3_2 = ggplot(pcm, aes(y = id, x = variable)) + 
  geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
  labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
  # scale_fill_manual(values = colours, guide = FALSE) + 
  scale_x_discrete(limits = rev(levels(pcm$variable))) + mytheme1 

g3_2


fileName =  paste(wlxpath,"/bubbleplot_diff.pdf",sep = "")
ggsave(fileName,g3_2, width = 12, height =15 )

#------- network of soil GC MS result#-------------
library(tidyverse)



netpath = paste(res1path,"/GCnetwork/",sep = "")
dir.create(netpath)
library(ggClusterNet)
library(igraph)
library(ggpubr)
result = network(ps = ps_rela,N = 0,r.threshold=0.6,p.threshold=0.05,label = FALSE,path = netpath ,zipi = TRUE,fill = "group"
                 )


# 全部样本的网络比对
g4_1 = result[[1]] + mytheme1
g4_1
# 全部样本网络参数比对
data = result[[2]]
plotname1 = paste(netpath,"/network_all.jpg",sep = "")
ggsave(plotname1, g4_1,width = 10,height = 6)
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, g4_1,width = 10,height = 6)

tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)

library(patchwork)
require(ggplotify)
g3_1 = as.ggplot(g3_1)
p_result1_2 <-(g1_1 + g2_1 + g1_1 )/( g4_1)



layout <- "
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
AAAAAACCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
BBBBBBCCCCCCC
DDDDDDDDDDDDD
DDDDDDDDDDDDD
DDDDDDDDDDDDD
"
p_result1_2 <- g2_1 + g1_1 + g3_1 + g4_1 +
  plot_layout(design = layout)
p_result1_2

plotname = paste(res1path,"/result1_2.pdf",sep = "")
ggsave(plotname, p_result1_2,width = 15,height = 10)


#-------env plot #-----------
envpath = paste(res1path,"/env/",sep = "")
dir.create(envpath)

data <- env
data$group = map$Group

dataenv <- data %>% select(ID,group, everything())

head(dataenv)
#
result = MuiKwWlx(data = dataenv,num = c(3:11))
result

FileName <- paste(envpath,"/alpha_diversity_all_abc.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(envpath,"/alpha_diversity_Sample.csv", sep = "")
write.csv(index,FileName,sep = "")

result1 = FacetMuiPlotresultBox(data = dataenv,num = c(3:11),result = result,sig_show ="abc",ncol = 5 )
p1_1 = result1[[1]] + scale_x_discrete(limits = axis_order) + 
  theme_bw() + mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1


res = FacetMuiPlotresultBar(data = dataenv,num = c(3:11),result = result,sig_show ="abc",ncol = 5)
p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) + theme_bw()  + 
  mytheme1+ guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = FacetMuiPlotReBoxBar(data = dataenv,num = c(3:11),result = result,sig_show ="abc",ncol = 5)
p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + theme_bw()  + mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(envpath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 10, height = 6)

FileName <- paste(envpath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = 10, height = 6)

FileName <- paste(envpath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = 10, height = 6)



#--env and microbial community RDA or CCA #--------
envpath = paste(res1path,"/env_micro/",sep = "")
dir.create(envpath)

library(EasyStat)
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/rda-cca.R")

row.names(env) = env$ID
env$ID= NULL


result = RDA_CCA(ps = ps,env = env,path = envpath)


#提取图片
p2_1 = result[[1]] + scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) + mytheme1

# 提取作图数据
dataplot = result[[2]]
# 提取带有标记的图片
p2_2 = result[[3]]  + scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) + mytheme1

#提取理化提供群落关系的渐检验结果
aov = result[[4]]

##保存
plotnamea = paste( envpath,"/RDA_env.pdf",sep = "")
ggsave(plotnamea, p2_1, width = 8, height = 6)
plotnamea4 = paste( envpath,"/RDA_env.jpg",sep = "")
ggsave(plotnamea4, p2_1, width = 8, height = 6)


filenamea = paste( envpath,"dataplot.txt",sep = "")
write.table(dataplot ,file=filenamea,sep="\t",col.names=NA)

filenamea = paste( envpath,"aov.txt",sep = "")
write.table(aov,file=filenamea,sep="\t",col.names=NA)

plotnamea = paste( envpath,"/RDA_envlabel.pdf",sep = "")
ggsave(plotnamea, p2_2, width = 8, height = 6)
plotnamea4 = paste( envpath,"/RDA_envlabel.jpg",sep = "")
ggsave(plotnamea4, p2_2, width = 8, height = 6)


#--变量相对丰度标准化编号：ps1_rela
ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 

otu = as.data.frame(t(vegan_otu(ps_rela)))
mapping = as.data.frame( sample_data(ps_rela))

env.dat = env
# match env and fg datasets
samp.fg = colnames(otu)
env.st = decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)#
samp.env= rownames(env.st)
my.env = match(samp.fg, samp.env)
env.st2 = na.omit(env.st[my.env, ])  # omit the NA rows if without fg data
samp.env= rownames(env.st2)
my.fg = match(samp.env, samp.fg)
otu = otu[, my.fg]

# for CCA calculation
otu = t(otu)
C.whole = rda(otu, env.st2)
C.whole

#----------------------------------------选择环境变量-----------------------------
inf_factor = vif.cca(C.whole)
inf_factor

# delete varable with max inflation factor
na_env = which(is.na(inf_factor))
if(isTRUE(length(na_env) > "0") ){
  inf_factor = inf_factor[-na_env]
}

max_env = which(inf_factor == max(inf_factor))
env.st4 = env.st2
while ( inf_factor[max_env] > 20){
  env.st4 = env.st4[,-max_env]
  C.reduced = cca(otu, env.st4)
  inf_factor = vif.cca(C.reduced)
  max_env = which(inf_factor == max(inf_factor))
}
output2 = inf_factor ;output2

C.whole = rda(otu, env.st4)  ##rda(otu, env.st3)
total.chi = C.whole$tot.chi
ind.p = array(0,dim=c(1,ncol(env.st4)))
for(j in 1:ncol(env.st4)){
  # j = 1
  ind.par = rda(otu, env.st4[,j], env.st4[,-j])
  ind.chi = ind.par$CCA$tot.chi
  ind.per = ind.chi/total.chi
  ind.p[j] = ind.per
}
ind.p

rowname = colnames(env.st4);rowname
out = matrix(data=NA,ncol=length(colnames(env.st4)),nrow=1);out
out = ind.p
rownames(out) = "percent"
colnames(out) = rowname
out
plotname = paste( envpath,"/vpa_out.csv",sep = "")
write.csv(out,plotname )
#------提取解释比例
total.chi = C.whole$tot.chi;total.chi
total.constrained = C.whole$CCA$tot.chi ; total.constrained
# 解释的比例
explained.percent = (total.constrained) / total.chi;explained.percent
# 未解释的比例
unexplained.percent = (total.chi - total.constrained) / total.chi;unexplained.percent

exp = data.frame(ID = c("explained.percent","unexplained.percent"),count = c(explained.percent,unexplained.percent))
exp
plotname = paste( envpath,"/vpa_explain.csv",sep = "")
write.csv(exp,plotname)




#----------network  for microbiol community, GC ms soil and env #------------

netpath = paste(res1path,"/three_netwotk_bios/",sep = "")
dir.create(netpath)
row.names(env) <- env$ID 
env$ID = NULL
data1  = as.data.frame(t(env))
data1$ID = row.names(data1)
data1 = data1 %>% select(ID,everything())
# data1 = data1 [1:7]
Gru = data.frame(SampleID = row.names(data1),Group = "env")


# 导入细菌ps，通过相对丰度value1来过滤otu表格
# ps_sub <- phyloseq::subset_samples(ps,Group %in% c("OF"));ps_sub
ps16 = ps
#导入真菌ps，通过相对丰度value1来过滤otu表格
ps_sub <- phyloseq::subset_samples(psG1,Group %in% c("OF"));ps_sub
psIT = psG1
ps.merge <- merge16S_ITS (ps16,psIT,N16s = 0.002,NITS = 0.002,
                          onlygroup = TRUE,
                          dat1.lab = "bac",
                          dat2.lab = "GC"
)
ps.merge

sample_data(ps.merge)

library(sna)
library(igraph)
library(ggpubr)


# result <- corBiostripe(data = data1,group = Gru,ps = ps.merge)
result <- corBionetwork(ps = ps.merge,
                        N = 0.001,
                        r.threshold = 0.90, # 相关阈值
                        p.threshold = 0.05,
                        label = FALSE,
                        group = "Group",
                        env = data1[-1], # 环境指标表格
                        envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        path = netpath,# 结果文件存储路径
                        fill = "group", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = TRUE, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = FALSE, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        height = 10
)

p = result[[1]]
p
# 全部样本网络参数比对
data = result[[2]]
plotname1 = paste(netpath,"/network_all.jpg",sep = "")
ggsave(plotname1, p,width = 12,height = 12)
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 12,height = 12)
tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)

library(ggrepel)
result <- corBionetwork(ps = ps.merge,
                        N = 0.001,
                        r.threshold = 0.90, # 相关阈值
                        p.threshold = 0.05,
                        label = T,
                        group = "Group",
                        env = data1[-1], # 环境指标表格
                        envGroup = Gru,# 环境因子分组文件表格
                        # layout = "fruchtermanreingold",
                        path = netpath,# 结果文件存储路径
                        fill = "elements", # 出图点填充颜色用什么值
                        size = "igraph.degree", # 出图点大小用什么数据
                        scale = TRUE, # 是否要进行相对丰度标准化
                        bio = TRUE, # 是否做二分网络
                        zipi = FALSE, # 是否计算ZIPI
                        step = 100, # 随机网络抽样的次数
                        width = 12,
                        height = 10
)



p = result[[1]]
p

plotnamea = paste( netpath,"/net_label.pdf",sep = "")
ggsave(plotnamea, p, width = 50, height = 40,limitsize = FALSE)
plotnamea4 = paste( netpath,"/net_label.jpg",sep = "")
ggsave(plotnamea4, p, width = 50, height = 40,limitsize = FALSE)
