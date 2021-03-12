

result_path <- "result_proj"
dir.create(result_path)

#---library R  package#--------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)

#---result1 analyses#-------------

#--creat path save result1
res1path <- paste(result_path,"/result6",sep = "")
dir.create(res1path)

# result 4 #---------
width = num
height = num

#--set theme of plot#-----
#--设置图形的样式和配色


#----set color or fill#------------
library(RColorBrewer)#调色板调用包
#调用所有这个包中的调色板
display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
display.brewer.pal(9,"Set1")
colset1 <- brewer.pal(9,"Set1")[c(6,5,4)]
colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
colset1 <- brewer.pal(9,"Set1")[c(5,6,4)]

#- leaf chemical propertity
env = read.csv("./dat_proj/leaf_chem/leaf_index.csv")
head(env)
env$ID = gsub("B","",env$ID)



#  order of x axis #------------
axis_order = c("Control","CF","OF")

#--1  picture of eaters eat leaf #-------


envpath = paste(res1path,"/env/",sep = "")
dir.create(envpath)

#---叶片营养数据（env）
library(tidyverse)
env = read.csv("./other_index/leaf_index.csv")
head(env)
colnames(env) = c("ID","weight","Salicylic_acid","Total_free_amino_acids","Nitrite","Nitrate","Soluble_protein","Soluble_sugar")
# env = env %>% filter(!ID %in% as.character(env$ID[1:6]))
row.names(env) = env$ID

env$group = c(rep("Control",6),rep("CF",6),rep("OF",6))
library(tidyverse)
dataenv <- env %>% 
  dplyr::select(ID,group, everything())

head(dataenv)
#
result = MuiKwWlx(data = dataenv,num = c(3))
result

FileName <- paste(envpath,"/alpha_diversity_all_abc.csv", sep = "")
write.csv(result,FileName,sep = "")
# FileName <- paste(envpath,"/alpha_diversity_Sample.csv", sep = "")
# write.csv(index,FileName,sep = "")

result1 = FacetMuiPlotresultBox(data = dataenv,num = c(3),result = result,sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]] + scale_x_discrete(limits = axis_order) + 
   mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1

FileName <- paste(envpath,"weight", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 5, height = 5)

## 水杨酸 #-----------
result = MuiKwWlx(data = dataenv,num = c(4))
result

FileName <- paste(envpath,"/水杨酸_alpha_diversity_all_abc.csv", sep = "")
write.csv(result,FileName,sep = "")
# FileName <- paste(envpath,"/alpha_diversity_Sample.csv", sep = "")
# write.csv(index,FileName,sep = "")

result1 = FacetMuiPlotresultBox(data = dataenv,num = c(4),result = result,sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]] + scale_x_discrete(limits = axis_order) + 
   mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
p1_1

FileName <- paste(envpath,"水杨酸", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 5, height = 5)



#--3  experiment of eater eat leaf #-----------
library(tidyverse)
library(EasyStat)
addpath = paste(res1path,"/add/",sep = "")
dir.create(addpath)
#  等待后期整合进来
dat0 <- read_csv("./dat_proj\\add-expriment/leaf_GC_expriment.csv")
dat1 = dat0[1:6]
dat1$group = as.factor(dat1$group)
result = MuiKwWlx(data = dat1,num = c(3:6))
result
result1 = FacetMuiPlotresultBox(data = dat1,num = c(3:6),result = result,sig_show ="abc",ncol = 2 )

# colset1 <- brewer.pal(9,"Set1")[c(5,6,4)]
p1 <- result1[[1]] + scale_x_discrete(limits = c("CK","Enriched","Depleted")) + theme_bw()  + mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = brewer.pal(9,"Set1")[c(6,5,4)])
p1

dat2 = dat0[c(1,7:9)]
colnames(dat2)[2] = "group"
result = MuiKwWlx(data = dat2,num = c(3:4))
result

result1 = FacetMuiPlotresultBox(data = dat2,num = c(3:4),result = result,sig_show ="abc",ncol = 1 )
p2 <- result1[[1]] + scale_x_discrete(limits = c("CK","CF","OF")) + theme_bw()  + mytheme1 + guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2

library(patchwork)

p0 <-  p2 + p1 +  plot_layout(widths = c(1, 2))

plotname = paste(addpath,"/add_expriment.pdf",sep = "")
ggsave(plotname, p0,width = 16,height = 8)
