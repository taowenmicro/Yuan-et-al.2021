


library(phyloseq)
library(ggClusterNet)
#--random forest for microbial communities#--------
ps <- readRDS("./dat_proj/Micro_seq_data/ps.rds")

map = as.data.frame(sample_data(ps))
map$Group = gsub("CK","Control",map$Group)
sample_data(ps) = map

library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
library(randomForest)
library(tidyverse)

set.seed(5)


result = MicroRF(ps = ps,group  = "Group3",optimal = 300,rfcv = F,nrfcvnum = 5,min = -1,max = 5)


datimp <- result[[5]]
datimp$id
#--按照丰度大小分配
ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela
ps_sub <- ps_rela %>%
  subset_taxa(
    row.names(tax_table(ps_rela ))%in%c(as.character(datimp$id))
  )
ps_sub

otu = as.data.frame((vegan_otu(ps_sub)))
head(otu)
otu$ID = row.names(otu)

datsum <- otu %>% 
  inner_join(as.data.frame(sample_data(ps_rela)))  %>%
  group_by(Group3) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

head(datsum)

i = 2
A = c()
for (i in 2:ncol(datsum)) {
  A[i] <- datsum$Group3[match(max(datsum[,i]),datsum[[i]])]
  
}

table(A)

otusel = data.frame(ID = colnames(datsum)[-1],group = A[-1])
head(otusel)

tax = as.data.frame(vegan_tax(ps_sub))
head(tax)
tax$ID = row.names(tax)

otu = as.data.frame(t(vegan_otu(ps_sub)))
head(otu)
otu$ID = row.names(otu)
otufia <- otusel %>%
  inner_join(otu) %>% 
  inner_join(tax)

id <- otufia %>% filter(group %in% c("insect"))
id$ID

#--特征OTU的可视化
ps_sub <- ps_rela%>%
  subset_taxa(
    row.names(tax_table(ps_rela))%in% id$ID
  )

otu = as.data.frame(t(otu_table(ps_sub)))
map = as.data.frame(sample_data(ps_sub))
map = map[,2]
otuGroup = merge(otu,map,by= "row.names",all = F)
head(otuGroup)
colnames(otuGroup)[1] = "ID"
library("reshape2")
datap = melt(otuGroup, id=c("ID","Group"))
head(datap)
library(ggplot2)
p = ggplot() + geom_boxplot(data = datap,aes(x = Group,y = value,color = Group)) +
  facet_wrap(~ variable,scales = "free",ncol = 2) 
p 

p = p + theme_bw()+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    
    plot.title = element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
    axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
    axis.text = element_text(size = 20,face = "bold"),
    axis.text.x = element_text(colour = "black",size = 14,angle = -60,hjust = -0.1),
    axis.text.y = element_text(colour = "black",size = 14),
    legend.text = element_text(size = 15,face = "bold")
  ) 
p




#---微生物的不同部位网络#-------
library(igraph)
library(network)
library(sna)
library(ggplot2)
library(ggrepel)
library(ggClusterNet)

otu = otufia[,3:68]
head(otu)
row.names(otu ) = otufia$ID
ps_m <- phyloseq(otu_table(otu,taxa_are_rows =T))
ps_m
#----------计算相关#----
result = corMicro (ps = ps_m,N = 0,r.threshold=0.8,p.threshold=0.05,method = "pearson")

#--提取相关矩阵
cor = result[[1]]
# head(cor)
ps_net = result[[3]]


table = as.data.frame(vegan_otu(ps_net))

netClu = data.frame(ID = otufia$ID,group = otufia$group)
# netClu$group = as.factor(sapply(strsplit(netClu$group , "_"), `[`, 2))
netClu$group = as.factor(netClu$group )
head(netClu)

##-----人工置顶半径大小和圆心位置

#--这里我设置r都是相同的，也可以设置不同，然后包装成一个向量就可以了#-------
xs = as.data.frame(table(netClu$group))
r = rep(15,length(xs$Freq))
r
#----准备圆心坐标，往往与你的设计有关#---------
# 有多少个模块就提供多少个坐标
#--指定坐标吮顺序按照一下指定
xs$Var1
#-------人工准备坐标
ax1 = c(120,0)
ax2 = c(130,-30)
ax3 = c(140,-70)
ax4 = c(130,-110)
# ax5 = c(120,-140)
# ax6 = c(12,0)
# ax7 = c(13,-30)
# ax8 = c(14,-70)
# ax9 = c(13,-110)
# ax10 = c(12,-140)
# ax11 = c(60,-140)
# da = rbind(ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11)
da = rbind(ax1,ax2,ax3,ax4)

#-------计算网络布局-得到节点坐标=node#---------
result2 = ArtifCluster(cor = cor,nodeGroup =netClu,r = r,da =da)

node = result2[[1]]
head(node)


otu_table = as.data.frame(t(vegan_otu(ps_net)))


#-----计算边#--------
edge = edgeBuild(cor = cor,plotcord = node)
node$ID = row.names(node)
nodes <- node %>% inner_join(netClu)
head(nodes)
head(edge)

#载入所需R包#---------


pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                         data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group),pch = 21,size = 5, data = nodes) + scale_fill_manual(values = brewer.pal(12,"Paired")) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

#---处理分开#——---------
netClu = data.frame(ID = otufia$ID,group = otufia$group)
# netClu$group = as.factor(sapply(strsplit(netClu$group , "_"), `[`, 2))
netClu$group = as.factor(netClu$group )
head(netClu)

##-----人工置顶半径大小和圆心位置

#--这里我设置r都是相同的，也可以设置不同，然后包装成一个向量就可以了#-------
xs = as.data.frame(table(netClu$group))
r = rep(15,length(xs$Freq))
r
#----准备圆心坐标，往往与你的设计有关#---------
# 有多少个模块就提供多少个坐标
#--指定坐标吮顺序按照一下指定
xs$Var1
#-------人工准备坐标
ax1 = c(120,0)
ax2 = c(130,-30)
ax3 = c(140,-70)
ax4 = c(130,-110)
ax5 = c(120,-140)
ax6 = c(12,0)
ax7 = c(13,-30)
ax8 = c(14,-70)
ax9 = c(13,-110)
ax10 = c(12,-140)
ax11 = c(60,-140)
da = rbind(ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11)
# da = rbind(ax1,ax2,ax3,ax4)

#-------计算网络布局-得到节点坐标=node#---------
result2 = ArtifCluster(cor = cor,nodeGroup =netClu,r = r,da =da)

node = result2[[1]]
head(node)


otu_table = as.data.frame(t(vegan_otu(ps_net)))


#-----计算边#--------
edge = edgeBuild(cor = cor,plotcord = node)
node$ID = row.names(node)
nodes <- node %>% inner_join(netClu)
head(nodes)
head(edge)

#载入所需R包#---------


pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group),pch = 21,size = 5, data = nodes) + scale_fill_manual(values = brewer.pal(12,"Paired")) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

