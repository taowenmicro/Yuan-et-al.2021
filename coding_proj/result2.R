
#---library R  package#--------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)

#---result1 analyses#-------------
result_path <- "result_proj"
dir.create(result_path)
#--creat path save result2
res1path <- paste(result_path,"/result2",sep = "")
dir.create(res1path)


#--inpurt data#-----------

#==GC
psG <- readRDS("./dat_proj/GC_data/GC_leaf//ps_leaf.rds")

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
library(ggClusterNet)
taxG = as.data.frame(vegan_tax(psG))

for (i in 1:length(taxG$group)) {
  if (taxG$group[i] == "short chain carbon organic acids") {
    taxG$group[i] = "small-molecule organic acids"
  }
}
table(taxG$group)

tax_table(psG) = as.matrix(taxG)


# rebuilding phylsoeq abject for unified sample name
psG <- phyloseq(otu_table(gctab_soil,taxa_are_rows = F),
                tax_table(psG),
                sample_data(map)
);psG

map = as.data.frame(sample_data(psG))
map$Group = gsub("CK","Control",map$Group)
sample_data(psG) = map

#--micro root
ps <- readRDS("./dat_proj/Micro_seq_data/ps.rds")

map = as.data.frame(sample_data(psG))
map$Group = gsub("CK","Control",map$Group)
# sample_data(ps) = map
ps_sub <- phyloseq::subset_samples(ps,Group %in% c("CK_root","CF_root","BOF_root"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
otu_table = as.data.frame(t(vegan_otu(ps_sub)))
#-- change name
colnames(otu_table) = gsub("B","",gsub("_root_","",colnames(otu_table)))
# rebuilding phylsoeq abject for unified sample name
ps_root <- phyloseq(otu_table(otu_table,taxa_are_rows = T),
               tax_table(ps),
               sample_data(map)
)
ps_root

#--micro leaf
ps <- readRDS("./dat_proj/Micro_seq_data/ps.rds")

ps_sub <- phyloseq::subset_samples(ps,Group %in% c("CK_leaf","CF_leaf","BOF_leaf"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
otu_table = as.data.frame(t(vegan_otu(ps_sub)))
#-- change name
colnames(otu_table) = gsub("B","",gsub("_leaf_","",colnames(otu_table)))
# rebuilding phylsoeq abject for unified sample name
ps_leaf <- phyloseq(otu_table(otu_table,taxa_are_rows = T),
                    tax_table(ps),
                    sample_data(psG)
)
ps_leaf


#- leaf chemical propertity
env = read.csv("./dat_proj/leaf_chem/leaf_index.csv")
head(env)
env$ID = gsub("B","",env$ID)


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
colset1 <- brewer.pal(9,"Set1")[c(5,6,4)]

colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))

map$Group


#  order of x axis #------------
axis_order = c("Control","CF","OF")

##--root microbiol community analyse #-----------
subpath <- paste(res1path,"/result_root_micro",sep = "")
dir.create(subpath)


#--alpha diversity#----------
alppath = paste(subpath,"/alpha/",sep = "")
dir.create(alppath)


source("G:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")

# 使用分面出图，得益于我开发的R包EasyStat
alp = alpha(ps = ps_root,inde="Shannon",group = "Group",Plot = TRUE )
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

result = BetaDiv(ps = ps_root, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)
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
result = barMainplot(ps = ps_root,j = "Phylum",rep = 6,axis_ord = NULL,label = FALSE ,sd = FALSE,Top = 10)

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
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")


# #---调整花瓣效果
p0_1 <- ggflower(ps = ps_root,
              rep = 3,
              group = "Group",
              start = 1, # 风车效果
              m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
              a = 0.2, # 花瓣胖瘦
              b = 1, # 花瓣距离花心的距离
              lab.leaf = 1, # 花瓣标签到圆心的距离
              col.cir = "green"
)
p0_2 = p0_1 +  scale_fill_manual(values = colset1,guide = F)

FileName1 <- paste(flowpath,"/flower",".pdf", sep = "")
ggsave(FileName1, p0_2, width = 6, height =6 )
FileName2 <- paste(barpath,"/flower",".jpg", sep = "")
ggsave(FileName2, p0_2, width = 6, height =6 )

#-- network #-----------

netpath = paste(subpath,"/network/",sep = "")
dir.create(netpath)



library(igraph)
library(sna)

?network
result = ggClusterNet::network(ps = ps_root,N = 0.001,r.threshold=0.6,p.threshold=0.05,label = FALSE,path = netpath ,zipi = TRUE)
# ?network
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

# library(patchwork)
# p_result1_1 <- p1_1 |(p2_3)|p3_1|p4_1
# p_result1_1 


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
EEEEEEEEEEEEE
"
p_result2_1 <- p1_1 + p2_3 + p3_2 + p0_2 + p4_1 +
  plot_layout(design = layout)

# sp1 <- p1_1 |(p2_3/p0_2)|p3_1
# sp1/p4_1

plotname = paste(res1path,"/result2.pdf",sep = "")
ggsave(plotname, p_result2_1 ,width = 18,height = 12)
plotname = paste(res1path,"/result2_1.pdf",sep = "")
ggsave(plotname, p_result2_1 ,width = 18,height = 15)

