

#--GC--soil#---------
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
psG <- phyloseq(otu_table(gctab_soil,taxa_are_rows = F),
                tax_table(psG),
                sample_data(map)
)

psGsoil <- psG %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    !Metabolite.name %in% c("Unknown","unknown"),
    # row.names(tax_table(ps_rela ))%in%c("SH010924.07FU_KF986690_reps_singleton","SH020983.07FU_JN235282_refs")
  )
psGsoil

#--GC-leaf#-----------------
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
# rebuilding phylsoeq abject for unified sample name
psG <- phyloseq(otu_table(gctab_soil,taxa_are_rows = F),
                tax_table(psG),
                sample_data(map)
);psG

map = as.data.frame(sample_data(psG))
map$Group = gsub("CK","Contral",map$Group)
sample_data(psG) = map

psGleaf <- psG %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    !Metabolite.name %in% c("Unknown","unknown"),
    # row.names(tax_table(ps_rela ))%in%c("SH010924.07FU_KF986690_reps_singleton","SH020983.07FU_JN235282_refs")
  )
psGleaf

#--合并
taxs = as.data.frame(vegan_tax(psGsoil))
taxs$Metabolite.name
length(taxs$Metabolite.name)
length(unique(taxs$Metabolite.name))
otu = as.data.frame(t(vegan_otu(psGsoil)))
head(otu)
otu$ID = row.names(otu)

taxssub = data.frame(ID = row.names(taxs),name = taxs$Metabolite.name,group = taxs$group) %>% 
  distinct(name,.keep_all = T)
otutaxs = otu %>% inner_join(taxssub)
head(otutaxs)
otutaxs$ID = NULL
colnames(otutaxs)[1:12] = paste(colnames(otutaxs)[1:12],"_soil",sep = "")



taxl = as.data.frame(vegan_tax(psGleaf))
length(taxl$Metabolite.name)
length(unique(taxl$Metabolite.name))
otu = as.data.frame(t(vegan_otu(psGleaf)))
head(otu)
otu$ID = row.names(otu)

taxlsub = data.frame(ID = row.names(taxl),name = taxl$Metabolite.name,group = taxl$group) %>% 
  distinct(name,.keep_all = T)
otutaxl = otu %>% inner_join(taxlsub)
head(otutaxl)
colnames(otutaxl)[1:18] = paste(colnames(otutaxl)[1:18],"_leaf",sep = "")
otutaxl$ID = NULL

# intersect(taxl$Metabolite.name,taxs$Metabolite.name)
otutax <- otutaxs %>% full_join(otutaxl,by = c("name","group"))
head(otutax)
#--合并为phyloseq对象#------
length(otutax$name)
length(unique(otutax$name))

otu = otutax[c(1:12,15:32)]
row.names(otu) = paste("GC_",1:nrow(otu),sep = "")
head(otu)
otu = as.matrix(otu)
otu[is.na(otu)] = 0

tax = data.frame(row.names = row.names(otu),name = otutax$name,class = otutax$group)
id = sapply(strsplit(colnames(otu) , "_"), `[`, 1)
group2 = substring(id, 1, 2)

zone = sapply(strsplit(colnames(otu) , "_"), `[`, 2)



map = data.frame(row.names = colnames(otu),ID = colnames(otu),Group =  paste(group2,zone,sep = "_"),Group2 = id,Group3 = zone )

ps= phyloseq(otu_table(otu,taxa_are_rows= T),
             tax_table(as.matrix(tax)),
             sample_data(map)
)


set.seed(5)
result = MicroRF(ps = ps,group  = "Group",optimal = 100,rfcv = F,nrfcvnum = 5,min = -1,max = 5)
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
  group_by(Group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

head(datsum)

i = 2
A = c()
for (i in 2:ncol(datsum)) {
  A[i] <- datsum$Group[match(max(datsum[,i]),datsum[[i]])]
  
}

table(A)

otusel = data.frame(ID = colnames(datsum)[-1],group = A[-1])
head(otusel)

tax = as.data.frame(vegan_tax(ps_sub))
tax$ID = row.names(tax)
head(tax)
colnames(tax)[2] = "class"
otu = as.data.frame(t(vegan_otu(ps_sub)))
head(otu)
otu$ID = row.names(otu)
otufia2 <- otusel %>%
  inner_join(otu) %>% 
  inner_join(tax)

head(otufia2)
head(otufia)
colnames(otufia)

#---准备两组环境数据#-------

envs = read.csv("./dat_proj/soil_chem/soil_chem_CF_OF.csv")
head(envs)
envl = read.csv("./dat_proj/leaf_chem/leaf_index.csv")
head(envl)
envl$ID = gsub("B","",envl$ID)
