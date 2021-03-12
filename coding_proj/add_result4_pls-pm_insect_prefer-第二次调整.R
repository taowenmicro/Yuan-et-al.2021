
#--pls-pm模型的调试将花费很长的时间，尤其是元素很多的情况下
# --开始构造代码，自动识别数据范围。


#--设计用于pls-pm 提取数据列的函数#-------------

coldsata <- function(a){
  d = cumsum(a)
  b <- cumsum(a) - a + 1
  A = list()
  for ( i  in 1:length(a)){
    A[[i]]  <-  seq(c(b[i],d[i])[1],c(b[i],d[i])[2], 1)
  }
  return(A)
}





library(EasyMicrobiome)
library(phyloseq)
#----不同指标之间的关系#-------
res1path = "./result_proj/result4/"
plspath = paste(res1path,"/pls-pm_insect_prefer/",sep = "")
dir.create(plspath)


## 数据准备#——------------
# 数据准备-本次主题需要使用的分析全部都准备到这里，方便随时使用，
#这在分析不公指标之间关系的时候显得非常方便。

#---土壤代谢组#----
library("vegan")
ps_GCsoil = readRDS("./GCMS_cab_36/小菜蛾实验叶土代谢组测样/土壤/ps_soil.rds")
re = otu_table(ps_GCsoil)
#-保留小数点位数达到取整目的
for (i in 1:dim(re)[2]) {
  re[,i] = round(re[,i],0)
}
#--将空缺值填充为0
re[is.na(re)] <- 0
otu_table(ps_GCsoil) <- re
ps_GCsoil = transform_sample_counts(ps_GCsoil, function(x) x / sum(x) )
# ps_GCsoil = filter_taxa(ps_GCsoil, function(x) sum(x ) > 0.01 , TRUE);ps_GCsoil #筛选序列数量大于1的

gctab_soil = as.data.frame(vegan_otu(ps_GCsoil))
map = sample_data(ps_GCsoil)
map
#--改名
rownames(gctab_soil) = as.character(map$sample_ID)
gctab_soil[1:6,1:6]
addtab = as.data.frame(matrix(0,ncol = ncol(gctab_soil),nrow = 6))
row.names(addtab) = paste("CK",1:6,sep = "")
colnames(addtab) = colnames(gctab_soil)
addtab[1:6,1:6]
gctab_soil = rbind(gctab_soil,addtab)
gctab_soil[1:18,1:6]


#--土壤化学性质#---------
envs = read.csv("./dat_proj/soil_chem/soil_chem_CF_OF.csv")
head(envs)
row.names(envs) = envs$ID
envs$ID = NULL
addtab = as.data.frame(matrix(0,ncol = ncol(envs),nrow = 6))
row.names(addtab) = paste("CK",1:6,sep = "")
colnames(addtab) = colnames(envs)
addtab[1:6,1:6]
envs = rbind(envs,addtab)
envs[1:18,1:6]
row.names(envs) = gsub("OF","BOF",row.names(envs))
envs = envs[match(row.names(gctab_soil),row.names(envs)),]


#-----叶片代谢组
library("vegan")
#下面开始做两样本相关性检测用otu表格和分泌物数据做
ps_GCleaf = readRDS("./GCMS_cab_36/小菜蛾实验叶土代谢组测样/植物//ps_leaf.rds")
re = otu_table(ps_GCleaf)
#-保留小数点位数达到取整目的
for (i in 1:dim(re)[2]) {
  re[,i] = round(re[,i],0)
}
#--将空缺值填充为0
re[is.na(re)] <- 0
otu_table(ps_GCleaf) <- re

ps_GCleaf  = transform_sample_counts(ps_GCleaf, function(x) x / sum(x) )
ps_sub <- phyloseq::subset_samples(ps_GCleaf,Group %in% c("CK","CF","BOF"));ps_sub
# ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE)#筛选序列数量大于1的
gctab_leaf = as.data.frame((vegan_otu(ps_sub)))
map = sample_data(ps_sub)
map
gctab_leaf = gctab_leaf[map$ID,]
#--改名比较麻烦
rownames(gctab_leaf) = as.character(map$sample_ID)
gctab_leaf[1:18,1:6]

# ## 使用pcoa排序结果代替叶片代谢组结果
# result = BetaDiv(ps = ps_sub, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)
# p1 = result[[1]]
# p1
# gctab_leaf_xy <- result[[2]][,1:2]
# rownames(gctab_leaf_xy) = as.character(map$sample_ID)
# colnames(gctab_leaf_xy) = paste("GCleaf",colnames(gctab_leaf_xy),sep = "")
# gctab_leaf = gctab_leaf_xy
# 
# gctab_leaf = gctab_leaf[match(row.names(gctab_soil),row.names(gctab_leaf)),]


#--挑选差异的分泌物
# 
id = read.csv("./result_proj/result2/result_leaf_GC/wlxtest/wlxtest_arti_fdr05RmUnknown.csv")
head(id)
# id <- c("RE_1282",
#   "RE_2194",
#   "RE_824",
#   "RE_2251",
#   "RE_435",
#   "RE_1270",
#   "RE_866",
#   "RE_1208",
#   "RE_668",
#   "RE_1338",
#   "RE_985",
#   "RE_2581",
#   "RE_1106",
#   "RE_2287",
#   "RE_817",
#   "RE_1300",
#   "RE_1424",
#   "RE_803",
#   "RE_747",
#   "RE_1371")


ps_sub <- ps_sub %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    # Species %in%c("Fusarium_oxysporum","Fusarium_keratoplasticum") 
    row.names(tax_table(ps_sub ))%in%c(id$id)
  )
ps_sub
# ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE)#筛选序列数量大于1的
gctab_leaf = as.data.frame((vegan_otu(ps_sub)))
map = sample_data(ps_sub)
map
gctab_leaf = gctab_leaf[map$ID,]
#--改名比较麻烦
rownames(gctab_leaf) = as.character(map$sample_ID)
gctab_leaf[1:18,1:6]
gctab_leaf = gctab_leaf[match(row.names(gctab_soil),row.names(gctab_leaf)),]




#---土壤微生物组数据
ps_micro = readRDS("./Seq_cab_66/data/ps.rds")
map = sample_data(ps_micro)
map
ps_sub <- phyloseq::subset_samples(ps_micro,Group %in% c("CF_soil","BOF_soil"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
ps_rela  = transform_sample_counts(ps_sub, function(x) x / sum(x) );ps_rela 
# ps_rela = filter_taxa(ps_rela, function(x) sum(x ) > 0.01 , TRUE);ps_rela #筛选序列数量大于1的

otu_table = as.data.frame(t(vegan_otu(ps_rela)))
map = sample_data(ps_rela)
map
#--改名
colnames(otu_table) = gsub("_soil_","",colnames(otu_table))

addtab = as.data.frame(matrix(0,ncol = 6,nrow = nrow(otu_table)))
colnames(addtab) = paste("CK",1:6,sep = "")
row.names(addtab) = row.names(otu_table)
addtab[1:6,1:6]
otu_table = cbind(otu_table,addtab)
otu_table = otu_table[,match(rownames(gctab_soil),colnames(otu_table))]
otu_table[1:6,1:18]


#--根系微生物
ps_sub <- phyloseq::subset_samples(ps_micro,Group %in% c("CK_root","CF_root","BOF_root"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
ps_rela  = transform_sample_counts(ps_sub, function(x) x / sum(x) );ps_rela
# ps_rela = filter_taxa(ps_rela, function(x) sum(x ) > 0.01 , TRUE);ps_rela #筛选序列数量大于1的

otu_table_root = as.data.frame(t(vegan_otu(ps_rela)))
colnames(otu_table_root)
colnames(otu_table_root) = gsub("_root_","",colnames(otu_table_root))
otu_table_root = otu_table_root[,match(rownames(gctab_soil),colnames(otu_table_root))]

otu_table_root[1:6,1:18]

#--叶片微生物
ps_sub <- phyloseq::subset_samples(ps_micro,Group %in% c("CK_leaf","CF_leaf","BOF_leaf"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
ps_rela  = transform_sample_counts(ps_sub, function(x) x / sum(x) );ps_rela 
# ps_rela = filter_taxa(ps_rela, function(x) sum(x ) > 0.01 , TRUE);ps_rela #筛选序列数量大于1的
otu_table_leaf = as.data.frame(t(vegan_otu(ps_rela)))
colnames(otu_table_leaf)
colnames(otu_table_leaf) = gsub("_leaf_","",colnames(otu_table_leaf))
otu_table_leaf = otu_table_leaf[,match(rownames(gctab_soil),colnames(otu_table_leaf))]
otu_table_leaf[1:6,1:18]

# # #--代表
# result = BetaDiv(ps = ps_rela, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)
# p1 = result[[1]]
# p1
# otu_table_leaf_xy <- result[[2]][,1:2]
# rownames(otu_table_leaf_xy) = gsub("_leaf_","",rownames(otu_table_leaf_xy))
# otu_table_leaf_xy = otu_table_leaf_xy[match(rownames(gctab_soil),rownames(otu_table_leaf_xy)),]
# otu_table_leaf_xy
# colnames(otu_table_leaf_xy) = paste("leaf",colnames(otu_table_leaf_xy),sep = "")
# otu_table_leaf <- as.data.frame(t(otu_table_leaf_xy))




#--虫子微生物
ps_sub <- phyloseq::subset_samples(ps_micro,Group %in% c("CK_insect","CF_insect","BOF_insect"));ps_sub
ps_sub = filter_taxa(ps_sub, function(x) sum(x ) > 1 , TRUE)#筛选序列数量大于1的
ps_rela  = transform_sample_counts(ps_sub, function(x) x / sum(x) );ps_rela 
# ps_rela = filter_taxa(ps_rela, function(x) sum(x ) > 0.01 , TRUE);ps_rela #筛选序列数量大于1的
otu_table_insect = as.data.frame(t(vegan_otu(ps_rela)))
colnames(otu_table_insect)
colnames(otu_table_insect) = gsub("_insect_","",colnames(otu_table_insect))
otu_table_insect = otu_table_insect[,match(rownames(gctab_soil),colnames(otu_table_insect))]
otu_table_insect[1:6,1:18]

# #--代表
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("G:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")
result = BetaDiv(ps = ps_rela, group = "Group", dist = "bray", method = "PCoA", Micromet = "adonis", pvalue.cutoff = 0.05)
p1 = result[[1]]
p1
otu_table_insect_xy <- result[[2]][,1:2]
rownames(otu_table_insect_xy) = gsub("_insect_","",rownames(otu_table_insect_xy))
otu_table_insect_xy = otu_table_insect_xy[match(rownames(gctab_soil),rownames(otu_table_insect_xy)),]
otu_table_insect_xy
colnames(otu_table_insect_xy) = paste("insect",colnames(otu_table_insect_xy),sep = "")
otu_table_insect <- as.data.frame(t(otu_table_insect_xy))



#---叶片营养数据（env）
library(tidyverse)
env = read.csv("./other_index/leaf_index.csv")
head(env)
colnames(env) = c("ID","weight","Salicylic_acid","Total_free_amino_acids","Nitrite","Nitrate","Soluble_protein","Soluble_sugar")
env = env %>% filter(!ID %in% as.character(env$ID[1:6]))
row.names(env) = env$ID
env$ID = NULL
decostand(env, MARGIN=1, "normalize")
env = decostand(env, MARGIN=1, "hellinger")
addtab = as.data.frame(matrix(0,ncol = ncol(env),nrow = 6))
row.names(addtab) = paste("CK",1:6,sep = "")
colnames(addtab) = colnames(env)
addtab[1:6,1:6]
env = rbind(env,addtab)
env[1:18,1:6]
env = env[match(row.names(gctab_soil),row.names(env)),]

#-----------plsPM--路径模型绘制#------
library("plspm")
pathplspm = read.csv("./plsPM路径_insect_prefer.csv",row.names = 1)
head(pathplspm)
pathplspm = as.matrix(pathplspm)
# 2、描绘路径
innerplot(pathplspm)
# 3、设定模块 即设定潜变量包含的显变量数据，这里我们一共有五个潜变量，也就是说这些潜变量是不可以直接通过测定的，只能有显变量来推断
df = cbind(
  t(otu_table),
  gctab_soil,
  gctab_leaf,
  env[2:7],
  t(otu_table_root),
  t(otu_table_leaf),
  t(otu_table_insect),
  envs,
  env[1]
)
df = as.matrix(df)
df[is.na(df)] = 0
df = as.data.frame(df)
dim(df)

a_dt <- c(dim(t(otu_table))[2],
          dim(gctab_soil)[2],
          dim(gctab_leaf)[2],
          dim(env[2:7])[2],
          dim(t(otu_table_root))[2],
          dim(t(otu_table_leaf))[2],
          dim(t(otu_table_insect))[2],
          dim(envs)[2],
          dim(env[1])[2]
)
sum(a_dt)
id = c("OTUsoil","GCsoil","GCleaf","env","OTUroot","OTUleaf","OTUinsect","env_soil","eat_prefer")
sat_blocks <- coldsata(a_dt)
names(sat_blocks) = id

# 4、设置模块数
# A表示是反映型的变量，隐变量是显变量的原因，箭头指向显变量
# B表示是影响型的变量，显变量是隐变量的原因，箭头指向隐变量
sat_mod <- rep("A", 9)


satpls <- plspm(df, pathplspm,sat_blocks, modes = sat_mod,
                scaled = T,boot.val = F)
#-------------------------------------------------------------
#-下面这个错误往往会在进行boot.val为  T的时候出现，-------
# Error in Path[k1, k2] <- path_lm$coef[-1, 1] :    
# 被替换的项目不是替换值长度的倍数
#-------------------------------------------------------------


#-----------模型细节查看
#路径效应指数
satpls$path_coefs
#显变量对隐变量的解释，这里可以查看weigt，比loading好一些
satpls$outer_model
#这里是R方
aa= summary(satpls)
aa$inner_summary
satpls$unidim
satpls$inner_summary[, "R2", drop = FALSE]

#显著性
satpls$inner_model

#提取模型拟合度
satpls$gof

innerplot(satpls)
#提取直接影响和间接影响
satpls$effects

#--并没有做boot检验
satpls$boot


#----结构虽然不是特别好，但是也不坏，这里探索一下可视化的方法，因为默认可视化实在是太过分了。

# 转化为相关边文件

#--这是R值矩阵
satpls$path_coefs


#--
i = 1

for (i in 1:length(names(satpls$inner_model))) {
  names(satpls$inner_model)[i]
  data = data.frame(from = row.names(as.data.frame(satpls$inner_model[[i]]))[-1],to = names(satpls$inner_model)[i],p = satpls$inner_model[[i]][-1,4])
  head(data)
  if (i == 1) {
    oridata = data
  }
  if (i != 1) {
    oridata = rbind( oridata,data)
  }
}
dim(oridata)
head(oridata)


#---------提取效应值矩阵
satpls$path_coefs


#----提取相关矩阵和p值矩阵
corr <- as.matrix(satpls$path_coefs)
# p.value <- as.matrix(occor.p )

corr

# 构造边文件
edges <- tibble::tibble(to = rep(row.names(corr), ncol(corr)),
                        from = rep(colnames(corr), each = nrow(corr)),
                        r = as.vector(corr))
as.data.frame(edges)
edges <- dplyr::filter(edges, as.vector(lower.tri(corr)))
#-----------合并数据框

head(edges)

edges <- edges %>% 
  inner_join(oridata,by = c("from","to")) 

as.data.frame(edges)
#--去除不显著的边
r.threshold=0
p.threshold=0.15
edges <- dplyr::filter(edges,p < p.threshold)
edges$r
library(tidyverse)
colnames(edges)[3] = "weight"
edges$weight = round(edges$weight,4)




#---设置边的正负
E.color <- edges$weight
edges$direction <- ifelse(edges$weight>0, "pp",ifelse(edges$weight<0, "np","ns"))
#----节点的属性--名称---门类七个--总体丰度信息，分组丰度信息--差异信息--设置自己需要标注信息都可以加入其中。
node = tibble::tibble(name = row.names(satpls$inner_summary))
row.names(node) = node$name
# node = data.frame(name = c(node$name[1:2],"GCleaf",node$name[3:7]),row.names = c(node$name[1:2],"GCleaf",node$name[3:7]))
node 

#--将R值加入到node中
a = satpls$inner_summary[, "R2", drop = FALSE]
node$R2 = paste("R2: ",round(a$R2,3),sep = "")
dim(
  node
)

tbl_graph = tidygraph::tbl_graph(nodes = node, edges = edges, directed = TRUE)
tbl_graph

library(ggraph)
p = ggraph(tbl_graph , layout = 'linear', circular = TRUE) +
  # geom_edge_arc(aes(edge_colour = direction,
  #                   edge_width = weight,
  #                   label = weight) +
  # geom_label(aes(x = x/2,y = y/2,label = name))+
  # geom_edge_fan() +
  geom_edge_link(aes(edge_colour = direction,
                     edge_width = weight,
                     label = weight,
                     start_cap = label_rect(from),
                     end_cap = label_rect(to)),
                 arrow = arrow(length = unit(8, 'mm')))+
  geom_node_point(aes())+
  geom_node_label(aes(label = name),vjust = -1) +
  geom_node_label(aes(label = R2),vjust = 1)  +
  theme_void()

p



filename = paste(plspath,"/PLSPM-change-R.pdf",sep = "")

ggsave(filename,p,width = 10,height = 8)

ggraph(tbl_graph , layout = 'linear', circular = TRUE) +
  geom_edge_arc(aes(edge_colour = direction,
                    # edge_width = weight,
                    label = weight,
                    start_cap = label_rect(from),
                    end_cap = label_rect(to)),
                arrow = arrow(length = unit(8, 'mm')),vjust = -1) +
  geom_node_point(aes())+
  geom_node_label(aes(label = name),vjust = -1,hjust = -1) +
  geom_node_label(aes(label = R2),vjust = 1)  +
  theme_void()


#--换一种圈--类似上一个，所以不做
# ggraph(tbl_graph , layout = 'linear', circular = TRUE) +
#   geom_edge_arc(aes(edge_colour = direction,
#                     # edge_width = weight,
#                     label = weight,
#                     start_cap = label_rect(from),
#                     end_cap = label_rect(to)),
#                 arrow = arrow(length = unit(8, 'mm')),vjust = -1) +
#   geom_node_point(aes())+
#   geom_node_label(aes(label = name),vjust = -1,hjust = -1) +
#   theme_void()
set.seed(2)
p = ggraph(tbl_graph , layout = 'graphopt') +
  geom_edge_link(aes(edge_colour = direction,
                     edge_width = weight,
                     label = weight,
                     start_cap = label_rect(from),
                     end_cap = label_rect(to)),
                 arrow = arrow(length = unit(8, 'mm')),vjust = -1) +
  geom_node_point(aes())+
  geom_node_label(aes(label = name),hjust = 0.1,vjust = -0.1) +
  geom_node_label(aes(label = R2),vjust = 1)  +
  theme_void()
p

filename = paste(interpath,"/PLSPM_graphopt_line.pdf",sep = "")
ggsave(filename,p,width = 10,height = 8)

filename = paste(interpath,"/PLSPM_edge.csv",sep = "")
write.csv(edges,filename)


# #-----------使用PC坐标轴做path分析失败#-----
# df = cbind(env,gctab_soil_xy,gctab_leaf_xy,otu_table_xy,otu_table_root_xy,otu_table_leaf_xy,otu_table_insect_xy)
# head(df)
# df = as.matrix(df)
# df[is.na(df)] = 0
# df = df ^2
# df = as.data.frame(df)
# colnames(path)
# 
# sat_blocks <-list(12:13,8:9,10:11,2:7,14:15,16:17,18:19,1)
# # example_scaling = list(c("NUM", "NUM"),
# #                        c("NUM", "NUM"),
# #                        c("NUM", "NUM"),
# #                        c("NUM", "NUM", "NUM", "NUM", "NUM", "NUM"),
# #                        c("NUM", "NUM"),
# #                        c("NUM", "NUM"),
# #                        c("NUM", "NUM"),
# #                        c("NUM")
# #                        )
# id = c("OTUsoil","GCsoil","GCleaf","env","OTUroot","OTUleaf","OTUinsect","eat_prefer")
# names(sat_blocks) = id
# sat_mod <- rep("A", 8)
# 
# satpls <- plspm(df, path,sat_blocks, modes = sat_mod,
#                 scaled = example_scaling,boot.val = F)
