

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
ggsave(plotname1, p,width = 12,height = 5)
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

otu = NULL
tax = NULL
map = NULL
ps = ps.merge
N = 0.001
r.threshold = 0.90 # 相关阈值
p.threshold = 0.05
label = FALSE
group = "Group"
env = data1[-1] # 环境指标表格
envGroup = Gru# 环境因子分组文件表格
# layout = "fruchtermanreingold",
path = netpath# 结果文件存储路径
fill = "group"# 出图点填充颜色用什么值
size = "igraph.degree"# 出图点大小用什么数据
scale = TRUE # 是否要进行相对丰度标准化
bio = TRUE# 是否做二分网络
zipi = FALSE# 是否计算ZIPI
step = 100# 随机网络抽样的次数
width = 12
height = 10
