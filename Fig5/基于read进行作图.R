rm(list=ls())
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(ggcor)
library(reshape2)
#基于read richness 是否符合ldg ddr
#根据基因的表，来计算
setwd('/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释')
data1<- read.table('all-archaea-ncyc.txt',header = T,row.names = 1,check.names = F)
data11<- read.table('all-archaea-mcyc.txt',header = T,row.names = 1,check.names = F)
data111 <- read.table('all-archaea-scyc.txt',header = T,row.names = 1,check.names = F)
metadata <- readxl::read_xlsx("sample.xlsx")
location1 <- read.table('Location',header = T,row.names = 1)

#处理下metadata的列名问题
metadata1 <- metadata[,-1]
row.names(metadata) <- metadata$id

#计算richness 
Richness <- specnumber(data1,MARGIN = 2)
data2 <- data.frame(Richness)
Richness1 <- specnumber(data11,MARGIN = 2)
data22 <- data.frame(Richness1)
Richness11 <- specnumber(data111,MARGIN = 2)
data222 <- data.frame(Richness11)
data_all <- cbind(data2,data22,data222)
color1 <- brewer.pal(12,"Paired")

df3 <- merge(location1,data_all,by="row.names",all = T)
colnames(df3)[3:5] <-c("Ncyc","Mcyc","Scyc")
df3$Location <- factor(df3$Location, levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang","Qingdao","Dongying","Dandong"))
 
ggplot(df3,aes(Location,Ncyc,fill=Location))+
  geom_violin()+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)
                                                                                                     
pdf("richness-ncyc.pdf",width = 10,height = 6)
ggplot(df3,aes(Location,Ncyc,fill=Location))+
  geom_violin()+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)
dev.off()

pdf("richness-mcyc.pdf",width = 10,height = 6)
ggplot(df3,aes(Location,Mcyc,fill=Location))+
  geom_violin()+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)
dev.off()

pdf("richness-scyc.pdf",width = 10,height = 6)
ggplot(df3,aes(Location,Scyc,fill=Location))+
  geom_violin()+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)
dev.off()


#补充下LDG和DDR的图
#LDG
#读入统计的数据，N,S,CH4
library(ggpubr)
library(ggpmisc)
library(reshape2)
location2 <- read.table('location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
colnames(data_all_ldg)[2:5] <-c("Ncyc","Mcyc","Scyc","Latitude")
data_all_ldg <- data_all_ldg[,-1]

data_new <- melt(data_all_ldg,id.vars = c("Latitude","Location"))
ggplot(data_new,aes(Latitude,log10(value),color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_cor(aes(color=variable))+
  theme_classic()+
  ylab("Richness")

#DDR
library(reshape2)
library(geosphere)
#计算similirty
#ncyc
bray_distance_data_all <- vegdist(t(data1),method = "bray")
data_ncyc <- as.matrix(bray_distance_data_all,row.names=1)
data_ncyc[upper.tri(data_ncyc)]=NA
data_ncyc <- melt(data_ncyc,na.rm = T)
data_ncyc$group <- "ncyc"
#mcyc
bray_distance_data_all <- vegdist(t(data11),method = "bray")
data_mcyc <- as.matrix(bray_distance_data_all,row.names=1)
data_mcyc[upper.tri(data_mcyc)]=NA
data_mcyc <- melt(data_mcyc,na.rm = T)
data_mcyc$group <- "mcyc"
#scyc
bray_distance_data_all <- vegdist(t(data111),method = "bray")
data_scyc <- as.matrix(bray_distance_data_all,row.names=1)
data_scyc[upper.tri(data_scyc)]=NA
data_scyc <- melt(data_scyc,na.rm = T)
data_scyc$group <- "scyc"

#合并三个模块的因子
data2 <- rbind(data_scyc,data_mcyc,data_ncyc)
write.table(data2,"bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#计算distance
location3 <- read.table('location-ddr',header = T,row.names = 1,check.names = F)
distance1 <-distm(location3,fun = distCosine)
rownames(distance1)<- rownames(location3)
colnames(distance1) <- rownames(location3)
data_distance1 <- as.matrix(distance1,row.names=1)
data_distance1[upper.tri(data_distance1)]=NA
data_distance2 <- melt(data_distance1,na.rm = T)
write.table(data_distance2,"bray_similirty_distance1.txt",quote = FALSE, row.names = FALSE,sep = "\t")


geo <- vegdist(location3,method = "euclidean")
geo<-pcnm(geo)
geo<-as.data.frame(geo$vectors)


#经过python脚本处理后，绘制图片
data1 <- read.table("bray_similarity_distance.txt",header = T)
library(ggplot2)
ggplot(data1,aes(x=log10(distance),y=log10(similarity),color=group))+
  geom_point(size=1)+
  geom_smooth(method = "lm",linewidth=2)+theme_classic()+xlim(5.25,6.5)+
  stat_cor(aes())+theme_bw()
  




#  stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = '~~~~')),label.y = "bottom", label.x = "left")






library(reshape2)
library(dplyr) 



#基于read的分类注释
setwd("G://work//宏基因组分析//基于read宏基因组注释//taxonomy//all")
all_files <- list.files("G://work//宏基因组分析//基于read宏基因组注释//taxonomy//all")
location <- read.table("G://work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0
data_all1 <-data_all1[-18,]

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"

dd <-merge(t(dd1),location,
          by = "row.names",)

color2 <- c("#FFFFB3","#8DD3C7",  "#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462" ,"#B3DE69", "#FCCDE5", 
            "#D9D9D9" ,"#CCEBC5","#BC80BD" )

c1 <- c("Others","CandidatusAenigmarchaeota","CandidatusWoesearchaeota","CandidatusHeimdallarchaeota","CandidatusThorarchaeota",
        "CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota", "Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
hh$Phylum = factor(hh$Phylum,order = T,levels =c1)

hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




plot_read <- ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))



plot_16s+plot_read






pdf("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//read-分类.pdf",width = 10,height = 6)
ggplot(hh)+
  geom_bar(aes(x =Location,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.5)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 65,colour = "black",family = "Times",size =15))
dev.off()






#基于分类的LDG作图
#补充下LDG和DDR的图
#LDG
#读入统计的数据，分类数据
library(ggpubr)
library(ggpmisc)
library(reshape2)

setwd("~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/all-species-number")
all_files <- list.files("~/Desktop//article1/宏基因组分析/基于read宏基因组注释//taxonomy//all-species-number/")
location <- read.table("~/Desktop//article1/宏基因组分析/基于read宏基因组注释/Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}
data_all1[is.na(data_all1)] <- 0

#进行数据抽平
#数据抽平
colSums(data_all1)
#使用该代码抽平
data_all_Flattening = as.data.frame(t(rrarefy(t(data_all1), min(colSums(data_all1)))))
tasv_data_all_Flatteningg <- t(data_all_Flattening)
Richness11 <- specnumber(data_all_Flattening,MARGIN = 2)
data_all <- data.frame(Richness11)
location2 <- read.table('/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释/location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
colnames(data_all_ldg)[3] <-c("Latitude")
data_all_ldg <- data_all_ldg[,-1]

data_new <- melt(data_all_ldg,id.vars = c("Latitude","Location"))
ggplot(data_new,aes(Latitude,value,))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_cor(aes())+theme_bw()+
  ylab("Richness")


#ddr作图

setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//")
bray_distance_tasv <- vegdist(tasv_data_all_Flatteningg,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
write.table(data2,"taxonomy_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#计算distance,前面算好了，直接借用
#读入python，处理好的文件
data1 <- read.table("~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/bray_similarity_distance.txt",header = T)
library(ggplot2)
ggplot(data1,aes(x=log10(distance),y=log10(similarity)))+
  geom_point(color="#8DEEEE")+
  geom_smooth(method = "lm",color="black")+theme_classic()+xlim(5.25,6.5)+stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )

#计算alpha 多样性

Richness11 <- specnumber(data_all_Flattening,MARGIN = 2)
data_all <- data.frame(Richness11)
data_all <- data.frame(Richness11)
location2 <- read.table('C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
data_all_ldg$Location <- factor(data_all_ldg$Location, levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                                          "Qingdao","Dongying","Dandong"))
ggplot(data_all_ldg,aes(Location,Richness11,fill=Location))+
  geom_violin(width=1.5)+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)



#B多样性分析,使用上面抽平过的数据
asv1 <- tasv_data_all_Flatteningg
distance <- vegdist(asv1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(asv1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location,by = "row.names")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(asv1 ~ Location,data=location,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)





p1<- ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))

ggsave(filename = "~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/read_pcoa.pdf",plot=p1,width=10,height =8 )




#N循环的途径，分类作图
library(reshape2)
library(dplyr) 


setwd("~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy//N-taxonomy/")
all_files <- list.files("~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy//N-taxonomy/",pattern="-taxonomy-number",full.names=TRUE)
location <- read.table("~/Desktop//article1/宏基因组分析//基于read宏基因组注释/Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"

dd <-merge(t(dd1),location,
           by = "row.names",)

color2 <- c("#FFFFB3","#8DD3C7",  "#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462" ,"#B3DE69", "#FCCDE5", 
            "#D9D9D9" ,"#CCEBC5","#BC80BD" )

c1<- c("Others","CandidatusMicrarchaeota","CandidatusWoesearchaeota","CandidatusHeimdallarchaeota",
       "CandidatusThorarchaeota","CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota","Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                                  "Qingdao","Dongying","Dandong"))
hh$Phylum = factor(hh$Phylum,order = T,levels = c1)



pdf("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//N-taxonomy//read-分类.pdf",width = 10,height = 6)
ggplot(hh)+
  geom_bar(aes(x =Location,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.5)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 65,colour = "black",family = "Times",size =15))
dev.off()

hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




plot_ncyc <- ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))





#beta 多样性
#beta 多样性
#这里采用的asv 应该是 基因和地点的表格作为第一数据

asv1 <- read.table("~/Desktop//article1/宏基因组分析//基于read宏基因组注释/t-all-archaea-ncyc.txt")
distance <- vegdist(asv1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(asv1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location,by = "row.names")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(asv1 ~ Location,data=location,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)

p1<-ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))


ggsave(filename = "~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/read_ncyc.pdf",plot=p1,width=10,height =8 )




#MCyc 分类作图


setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//MCyc-taxonomy/")
all_files <- list.files("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//MCyc-taxonomy/", pattern="-taxonomy-number",full.names=TRUE)
location <- read.table("C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"

dd <-merge(t(dd1),location,
           by = "row.names",)

color2 <- c("#FFFFB3","#8DD3C7",  "#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462" ,"#B3DE69", "#FCCDE5", 
            "#D9D9D9" ,"#CCEBC5","#BC80BD" )

c1<- c("Others","CandidatusWoesearchaeota","CandidatusAenigmarchaeota","CandidatusHeimdallarchaeota",
       "CandidatusThorarchaeota","CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota","Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
hh$Phylum = factor(hh$Phylum,order = T,levels = c1)



pdf("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//MCyc-taxonomy//read-分类.pdf",width = 10,height = 6)
ggplot(hh)+
  geom_bar(aes(x =Location,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.5)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 65,colour = "black",family = "Times",size =15))
dev.off()


hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




plot_mcyc <- ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))










#beta

#B多样性分析,基因的相关数据

asv1 <- read.table("~/Desktop//article1/宏基因组分析//基于read宏基因组注释/t-all-archaea-mcyc.txt")
distance <- vegdist(asv1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(asv1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)

data <- merge(plot_data,location,by = "row.names")

#基于bray-curtis距离进行计算
asv1.div <- adonis2(asv1 ~ Location,data=location,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)

p1<-ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))


ggsave(filename = "~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/read_mcyc.pdf",plot=p1,width=10,height =8 )













#SCyc 分类作图


setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//SCyc-taxonomy/")
all_files <- list.files("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//SCyc-taxonomy/", pattern="-taxonomy-number",full.names=TRUE)
location <- read.table("C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"

dd <-merge(t(dd1),location,
           by = "row.names",)

color2 <- c("#FFFFB3","#8DD3C7",  "#BEBADA" ,"#FB8072" ,"#80B1D3", "#FDB462" ,"#B3DE69", "#FCCDE5", 
            "#D9D9D9" ,"#CCEBC5","#BC80BD" )

c1<- c("Others", "CandidatusKorarchaeota","CandidatusAenigmarchaeota" ,"CandidatusHeimdallarchaeota","CandidatusThorarchaeota","CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota","Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
hh$Phylum = factor(hh$Phylum,order = T,levels = c1)



pdf("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//SCyc-taxonomy//read-分类.pdf",width = 10,height = 6)
ggplot(hh)+
  geom_bar(aes(x =Location,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.5)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 65,colour = "black",family = "Times",size =15))
dev.off()



hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




plot_scyc <- ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))








#beta 多样性
asv1 <- read.table("~/Desktop//article1/宏基因组分析//基于read宏基因组注释/t-all-archaea-scyc.txt")
distance <- vegdist(asv1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(asv1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)

data <- merge(plot_data,location,by = "row.names")

#基于bray-curtis距离进行计算
asv1.div <- adonis2(asv1 ~ Location,data=location,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)

p1<-ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))


ggsave(filename = "~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/read_scyc.pdf",plot=p1,width=10,height =8 )





plot_ncyc+plot_mcyc+plot_scyc

















#kegg作图
setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//kegg/")
all_files <- list.files("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//kegg/",pattern="-kegg-final",full.names=TRUE)
location <- read.table("C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {  
  df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0

#计算alpha 多样性
Richness <- specnumber(data_all1,MARGIN = 2)
Richness <- data.frame(Richness)
data_all <- data.frame(Richness)
data_all <- data.frame(Richness)
location2 <- read.table('C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
data_all_ldg$Location <- factor(data_all_ldg$Location, levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                                                "Qingdao","Dongying","Dandong"))
ggplot(data_all_ldg,aes(Location,Richness,fill=Location))+
  geom_violin(width=1.5)+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)



#kegg ldr 
library(vegan)
library(ggpubr)
library(ggpmisc)
library(reshape2)

Richness1 <- specnumber(data_all1,MARGIN = 2)
data1 <- data.frame(Richness1)

data_all_ldg <- merge(data1,location2,by = "row.names")
colnames(data_all_ldg)[3] <-c("Latitude")
data_all_ldg <- data_all_ldg[,-1]
data_new <- melt(data_all_ldg,id.vars = c("Latitude","Location"))
ggplot(data_new,aes(Latitude,value))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = '~~~~')),label.y = "bottom", label.x = "left")+
  theme_bw()+
  ylab("Richness")




#kegg ddr 计算
#DDR
library(reshape2)
library(geosphere)
setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//kegg/")
#计算similirty
tdata_all1 <- t(data_all1)
bray_distance_tdata_all1 <- vegdist(tdata_all1,method = "bray")
data1 <- as.matrix(bray_distance_tdata_all1,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
#写入相似性，后续python 处理
write.table(data2,"bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#计算distance
location3 <- read.table('C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\location-ddr',header = T,row.names = 1,check.names = F)
distance1 <-distm(location3,fun = distCosine)
rownames(distance1)<- rownames(location3)
colnames(distance1) <- rownames(location3)
data_distance1 <- as.matrix(distance1,row.names=1)
data_distance1[upper.tri(data_distance1)]=NA
data_distance2 <- melt(data_distance1,na.rm = T)
#写入距离，后续python 处理，脚本（read_bray_curtis）,还需要把QD改成CXD在sublime text 里面
write.table(data_distance2,"bray_similirty_distance1.txt",quote = FALSE, row.names = FALSE,sep = "\t")

#作图
#进入文本，写入ID
data1 <- read.table("bray_similarity_distance.txt",header = T)
library(ggplot2)
ggplot(data1,aes(x=log10(distance),y=log10(similarity)))+
  geom_point(color="#8DEEEE")+
  geom_smooth(method = "lm",color="black")+theme_classic()+xlim(5.25,6.5)+stat_cor(aes())


#beta 多样性
asv1<- read.table('~/Desktop//article1//宏基因组分析//基于read宏基因组注释//kegg/t-all-kegg.txt',header = T,row.names = 1,check.names = F)
distance <- vegdist(asv1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(asv1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)

data <- merge(plot_data,location,by = "row.names")

#基于bray-curtis距离进行计算
asv1.div <- adonis2(asv1 ~ Location,data=location,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)

p1<-ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))


ggsave(filename = "~/Desktop//article1/宏基因组分析//基于read宏基因组注释//taxonomy/read_kegg.pdf",plot=p1,width=10,height =8 )













#测试代码
data1<-read.table(file ="BH10_FDSW220017130-1r",header = T,row.names = 1,check.names = F)
data2<-read.table(file ="BH12_FDSW220017132-1r",header = T,row.names = 1,check.names = F)
data3<-read.table(file ="BH14_FDSW220017134-1r",header = T,row.names = 1,check.names = F)
data4<-read.table(file ="BH15_FDSW220017135-1r",header = T,row.names = 1,check.names = F)
data5<-read.table(file ="BH2_FDSW220017122-1r",header = T,row.names = 1,check.names = F)
data6<-read.table(file ="BH4_FDSW220017124-1r",header = T,row.names = 1,check.names = F)
data7<-read.table(file ="BH6_FDSW220017126-1r",header = T,row.names = 1,check.names = F)
data8<-read.table(file ="BH8_FDSW220017128-1r",header = T,row.names = 1,check.names = F)
data_all1 <- merge(data1,data2,by="row.names",all = T)
rownames(data_all1) <- data_all1$Row.names
data_all1 <- data_all1[,-1]
v <- c(3:8)
for (i in v) {
  df.now <- get(paste0("data",i))
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  
}

data_all1[is.na(data_all1)] <- 0





#载入包，读入数据
setwd('C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释')
data1<- read.table('all-archaea-ncyc.txt',header = T,row.names = 1,check.names = F)
data11<- read.table('all-archaea-mcyc.txt',header = T,row.names = 1,check.names = F)
data111 <- read.table('all-archaea-scyc.txt',header = T,row.names = 1,check.names = F)
metadata <- readxl::read_xlsx("sample.xlsx")
location1 <- read.table('Location',header = T,row.names = 1)



#N循坏
setwd('C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释')
data1<- read.table('all-archaea-ncyc.txt',header = T,row.names = 1,check.names = F)
#B多样性分析
tdata1 <- t(data1)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(df1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
group1 <-metadata['Location']

data <- merge(plot_data,metadata,by = "id")

ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))



#M-B多样性分析
tdata1 <- t(data11)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(df1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
group1 <-metadata['Location']

data <- merge(plot_data,metadata,by = "id")

ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))



#s-B多样性分析
tdata1 <- t(data111)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(df1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
group1 <-metadata['Location']

data <- merge(plot_data,metadata,by = "id")

ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))









#统计N基因数量
setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//all")
location <- read.table("C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")
data_all1<- read.table('C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\all-archaea-ncyc.txt',header = T,row.names = 1,check.names = F)

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"


dd <-merge(t(dd1),location,
           by = "row.names",)

color2 <- brewer.pal(11,"Spectral")


c1 <- c("Others","CandidatusAenigmarchaeota","CandidatusWoesearchaeota","CandidatusHeimdallarchaeota","CandidatusThorarchaeota",
        "CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota", "Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))





hh$Phylum = factor(hh$Phylum,order = T,levels = rev(c(row.names(data_all2)[1:10],"Others")))

hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

area_ncyc <- ggplot(hh3,aes(x =id,y = Abundance,fill = Phylum))+
  geom_area(position = position_fill())+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))



#统计M(CH4)
#统计M基因数量
setwd("C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//taxonomy//all")
location <- read.table("C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")
data_all1<- read.table('C:\\Users\\DELL\\Desktop\\work\\宏基因组分析\\基于read宏基因组注释\\all-archaea-mcyc.txt',header = T,row.names = 1,check.names = F)

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"


dd <-merge(t(dd1),location,
           by = "row.names",)

color2 <- brewer.pal(11,"Spectral")


c1 <- c("Others","CandidatusAenigmarchaeota","CandidatusWoesearchaeota","CandidatusHeimdallarchaeota","CandidatusThorarchaeota",
        "CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota", "Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))





hh$Phylum = factor(hh$Phylum,order = T,levels = rev(c(row.names(data_all2)[1:10],"Others")))

hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

area_mcyc <- ggplot(hh3,aes(x =id,y = Abundance,fill = Phylum))+
  geom_area(position = position_fill())+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))


#统计S基因数量
setwd("G://work//宏基因组分析//基于read宏基因组注释//taxonomy//all")
location <- read.table("G:\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
color1<-brewer.pal(12,"Set3")
data_all1<- read.table('G:\\work\\宏基因组分析\\基于read宏基因组注释\\all-archaea-scyc.txt',header = T,row.names = 1,check.names = F)

#按照数量大小排序：
order<-sort(rowSums(data_all1[,2:ncol(data_all1)]),index.return=TRUE,decreasing=T) 
data_all2<-data_all1[order$ix,]

dd1<-rbind(colSums(data_all2[11:as.numeric(length(rownames(data_all2))),]),data_all2[10:1,])
rownames(dd1)[1]<-"Others"


dd <-merge(t(dd1),location,
           by = "row.names",)

color2 <- brewer.pal(11,"Spectral")


c1 <- c("Others","CandidatusAenigmarchaeota","CandidatusWoesearchaeota","CandidatusHeimdallarchaeota","CandidatusThorarchaeota",
        "CandidatusLokiarchaeota","Crenarchaeota","CandidatusThermoplasmatota","CandidatusBathyarchaeota", "Euryarchaeota","Thaumarchaeota")
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')

kk<-tidyr::gather(dd,Phylum,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))





hh$Phylum = factor(hh$Phylum,order = T,levels = rev(c(row.names(data_all2)[1:10],"Others")))

hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")




ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.9)+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

area_scyc <- ggplot(hh3,aes(x =id,y = Abundance,fill = Phylum))+
  geom_area(position = position_fill())+
  scale_fill_manual(values = color2) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))




area_ncyc+area_mcyc+area_scyc


















#计算 kegg 的richness
setwd("~/Desktop//article1//宏基因组分析//基于read宏基因组注释//kegg/")
all_files <- list.files("~/Desktop//article1//宏基因组分析//基于read宏基因组注释/kegg/",pattern="-kegg-final",full.names=TRUE)
location <- read.table("~/Desktop//article1//宏基因组分析//基于read宏基因组注释/Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {  
    df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0

#nsc+kegg 统一绘图

setwd('~/Desktop/article1/宏基因组分析/基于read宏基因组注释/')
data1<- read.table('all-archaea-ncyc.txt',header = T,row.names = 1,check.names = F)
data11<- read.table('all-archaea-mcyc.txt',header = T,row.names = 1,check.names = F)
data111 <- read.table('all-archaea-scyc.txt',header = T,row.names = 1,check.names = F)
metadata <- readxl::read_xlsx("sample.xlsx")
location1 <- read.table('Location',header = T,row.names = 1)

#处理下metadata的列名问题
metadata1 <- metadata[,-1]
row.names(metadata) <- metadata$id

#计算richness 
Richness <- specnumber(data1,MARGIN = 2)
data2 <- data.frame(Richness)
Richness1 <- specnumber(data11,MARGIN = 2)
data22 <- data.frame(Richness1)
Richness11 <- specnumber(data111,MARGIN = 2)
data222 <- data.frame(Richness11)
Richness111 <- specnumber(data_all1,MARGIN = 2)
data2222 <- data.frame(Richness111)


data_all <- cbind(data2,data22,data222,data2222)
color1 <- brewer.pal(12,"Paired")

df3 <- merge(location1,data_all,by="row.names",all = T)
colnames(df3)[3:6] <-c("Ncyc","Mcyc","Scyc","kegg")
df3$Location <- factor(df3$Location, levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang","Qingdao","Dongying","Dandong"))

#补充下LDG和DDR的图
#LDG
#读入统计的数据，N,S,CH4
library(ggpubr)
library(ggpmisc)
library(reshape2)
location2 <- read.table('location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
colnames(data_all_ldg)[2:6] <-c("Ncyc","Mcyc","Scyc","kegg","Latitude")
data_all_ldg <- data_all_ldg[,-1]

data_new <- melt(data_all_ldg,id.vars = c("Latitude","Location"))
ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_cor(aes(color=variable))+theme_bw()+
  ylab("Richness")


ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+stat_poly_eq(formula = y ~ x,
             aes(label = paste(..eq.label.., ..rr.label..,p.value.label,sep ="~~~")),
             parse = TRUE)



upper <- ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_cor(aes(color=variable))+theme_classic()+
  ylab("")+
  xlab("")+
  coord_cartesian(ylim=c(3000,4000))+
  scale_y_continuous(breaks = c(3000, 4000, 250)) + # 以250为单位划分Y轴
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank())
  
  
down <- ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_cor(aes(color=variable))+theme_classic()+
  ylab("Richness")+
  coord_cartesian(ylim=c(0,200))


p1<- ggarrange(upper,down,ncol=1,align = "hv",common.legend = T)

ggsave("NCS_KEGG_ldg_2024.3.22.pdf",p1,width = 8,height = 6)


#DDR
library(reshape2)
library(geosphere)
#计算similirty
#ncyc
bray_distance_data_all <- vegdist(t(data1),method = "bray")
data_ncyc <- as.matrix(bray_distance_data_all,row.names=1)
data_ncyc[upper.tri(data_ncyc)]=NA
data_ncyc <- melt(data_ncyc,na.rm = T)
data_ncyc$group <- "ncyc"
#mcyc
bray_distance_data_all <- vegdist(t(data11),method = "bray")
data_mcyc <- as.matrix(bray_distance_data_all,row.names=1)
data_mcyc[upper.tri(data_mcyc)]=NA
data_mcyc <- melt(data_mcyc,na.rm = T)
data_mcyc$group <- "mcyc"
#scyc
bray_distance_data_all <- vegdist(t(data111),method = "bray")
data_scyc <- as.matrix(bray_distance_data_all,row.names=1)
data_scyc[upper.tri(data_scyc)]=NA
data_scyc <- melt(data_scyc,na.rm = T)
data_scyc$group <- "scyc"
#kegg
bray_distance_data_all <- vegdist(t(data_all1),method = "bray")
data_kegg <- as.matrix(bray_distance_data_all,row.names=1)
data_kegg[upper.tri(data_kegg)]=NA
data_kegg <- melt(data_kegg,na.rm = T)
data_kegg$group <- "kegg"






#合并三个模块的因子
data2 <- rbind(data_scyc,data_mcyc,data_ncyc,data_kegg)
write.table(data2,"ncs-kegg-bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#计算distance
location3 <- read.table('location-ddr',header = T,row.names = 1,check.names = F)
distance1 <-distm(location3,fun = distCosine)
rownames(distance1)<- rownames(location3)
colnames(distance1) <- rownames(location3)
data_distance1 <- as.matrix(distance1,row.names=1)
data_distance1[upper.tri(data_distance1)]=NA
data_distance2 <- melt(data_distance1,na.rm = T)
write.table(data_distance2,"bray_similirty_distance1.txt",quote = FALSE, row.names = FALSE,sep = "\t")

#经过python脚本处理后，绘制图片
data1 <- read.table("ncs-kegg-bray_similarity_distance.txt",header = T)
library(ggplot2)
p1<- ggplot(data1,aes(x=log10(distance),y=log10(similarity),color=group))+
  geom_point(size=1)+
  geom_smooth(method = "lm",linewidth=1.5)+theme_classic()+xlim(5.25,6.5)+
  stat_cor(aes())+
  scale_color_manual(values = c("#A07DB6","#78A72B","#EB736A","#1EB4B8"))+
  labs(x="Geographical distance (log10 scale)",y="Similarity (log10 scale)")

ggsave("nms_kegg_ddr_2024.3.22.pdf",p1,width = 8,height = 6)


#  stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = '~~~~')),label.y = "bottom", label.x = "left")



#16s+read作图

library(tidyverse)
library(ggprism)
library(vegan)
library(RColorBrewer)
library(reshape2)
#读入数据
setwd("G://work/16s分析/16s-all")
asv<- read.table("ASVs_counts_clean-rdp.txt",sep = '\t')
taxa1<- read.table("ASVs_taxonomy-clean-rdp.txt", sep="\t")
group<- readxl::read_xlsx("sample_df.xlsx")
rownames(group) <- group$id 

#进行抽平在进行分析
colSums(asv)
#使用该代码抽平
asv_Flattening = as.data.frame(t(rrarefy(t(asv), min(colSums(asv)))))

tasv_Flattening <- t(asv_Flattening)
#查看抽平后的每个样本的和
colSums(asv_Flattening)
#将抽平后的otu表保存到该工作目录，准备后面的多样性分析
asv <- asv_Flattening



Richness <- specnumber(asv,MARGIN = 2)
richness <- data.frame(Richness)
location2 <- read.table('G://work//宏基因组分析//基于read宏基因组注释//location-ldg.txt',header = T,row.names = 1,check.names = F)

data_all_ldg <- merge(richness,location2,by = "row.names")
colnames(data_all_ldg)[3] <-c("Latitude")
data_all_ldg <- data_all_ldg[,-1]








library(ggpubr)
library(ggpmisc)
library(reshape2)

setwd("~/Desktop/article1/宏基因组分析/基于read宏基因组注释//taxonomy//all-species-number/")
all_files <- list.files("~/Desktop/article1/宏基因组分析//基于read宏基因组注释//taxonomy//all-species-number/")
location <- read.table("~/Desktop/article1/宏基因组分析/基于read宏基因组注释/Location")
color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}
data_all1[is.na(data_all1)] <- 0

#进行数据抽平
#数据抽平
colSums(data_all1)
#使用该代码抽平
data_all_Flattening = as.data.frame(t(rrarefy(t(data_all1), min(colSums(data_all1)))))
#查看抽平结果
colSums(data_all_Flattening)


tasv_data_all_Flatteningg <- t(data_all_Flattening)

Richness11 <- specnumber(data_all_Flattening,MARGIN = 2)
data_all <- data.frame(Richness11)


Richness111 <- specnumber(data_all1,MARGIN = 2)
data_all_2 <- data.frame(Richness111)



location2 <- read.table('C://Users//DELL//Desktop//work//宏基因组分析//基于read宏基因组注释//location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
colnames(data_all_ldg)[3] <-c("Latitude")
data_all_ldg <- data_all_ldg[,-1]







#经过上述整合放在txt文件里
library(ggplot2)
read_16s <- read.table("~/Desktop/article1/fig2/16s_read_ldg-3",header = T,row.names = 1,check.names = F)
color1 <- c("#1f77b4","#f6b5d1")
p1<-ggplot(read_16s,aes(read_16s$Latitude,read_16s$Richness,color=read_16s$group))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )+theme_bw()+
  scale_color_manual(values = color1)


ggsave("~/Desktop/article1/16s分析/16s-all/16s_reads_ldg_2025.3.13.pdf",p1,width = 8,height = 6)

#DDR
#计算read ddr



bray_distance_tasv <- vegdist(tasv_data_all_Flatteningg,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
write.table(data2,"C://Users//DELL//Desktop//work//article1/fig2//taxonomy_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")

# 计算16s ddr

tasv <- t(asv)
bray_distance_tasv <- vegdist(tasv,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
data(varespec)
write.table(data2,"C://Users//DELL//Desktop//work//article1/fig2//16s_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")

#读入相应的文件，作ddr分析

data1 <- read.table("~/Desktop/article1/fig2/16s_bray_similarity_distance.txt",header = T)
data1$group <- "16s"
data2 <- read.table("~/Desktop/article1/fig2//read_bray_similarity_distance.txt",header = T)
data2$group <- "read"
data3 <- rbind(data1,data2)

color1 <- c("#1f77b4","#f6b5d1")

library(ggplot2)
p1<-ggplot(data3,aes(x=log10(distance),y=log10(similarity),color=group))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+xlim(5.25,6.5)+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label..,p.value.label,sep ="~~~")),parse = TRUE)+
  scale_color_manual(values = color1)

ggsave("/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释/16s_reads_ldg_25.3.13.pdf",p1,width = 8,height = 6)




#基于reads 的mantel-test 分析
#基于read的分类注释
library(dplyr)
setwd("~/Desktop/article1///宏基因组分析//基于read宏基因组注释//all-number/")
all_files <- list.files("~/Desktop/article1///宏基因组分析//基于read宏基因组注释/all-number/",pattern = "-number")
location <- read.table("~/Desktop/article1///宏基因组分析//基于read宏基因组注释//Location")
env1  <- read.table("~/Desktop/article1///宏基因组分析//基于read宏基因组注释/env.txt",header = T,row.names = 1,sep="\t")

color1<-brewer.pal(12,"Set3")

v<-length(all_files)
for (i in (1:v)){
  if (i<2){
    print(i)
    data1<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
  }
  else if(i==2){
    data2<-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
    print(data1)
    data_all1 <- merge(data1,data2,by="row.names",all = T)
    rownames(data_all1) <- data_all1$Row.names
    data_all1 <- data_all1[,-1]
  }
  else {   df.now <-read.table(file=all_files[i],header=T,row.names = 1,check.names = F)
  data_all1 <- merge(data_all1,df.now,by="row.names",all = T)
  rownames(data_all1) <- data_all1$Row.names
  data_all1 <- data_all1[,-1]
  }
}

data_all1[is.na(data_all1)] <- 0

#数据抽平
colSums(data_all1)



#使用该代码抽平
asv_Flattening = as.data.frame(t(rrarefy(t(data_all1), min(colSums(data_all1)))))

tasv_Flattening <- t(asv_Flattening)
#查看抽平后的每个样本的和
colSums(asv_Flattening)
#将抽平后的otu表保存到该工作目录，准备后面的多样性分析
asv_Flattening
tasv<- t(asv_Flattening)

asv1 <- tasv[match(rownames(env1),rownames(tasv)),]
#环境因子相关性分析及展示：
env <-data.matrix(env1)
qcorrplot(correlate(env), type = "lower",method = "spearman") +
  geom_square()+
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')


df_mantel <- mantel_test(asv1, env, mantel_fun = 'mantel',
                         spec_dist = 'bray', 
                         env_dist = 'euclidean')




#定义标签：
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.25, Inf),
                    labels = c("< 0.1", "0.1 - 0.25", ">= 0.25")),#定义Mantel的R值范围标签
         df_p = cut(p, breaks = c(-Inf,0.001,0.01, 0.05, Inf),
                    labels = c("< 0.001", "0.001<x<0.01", "0.01<x<0.05",">= 0.05")))#定义Mantel的P值范围标



qcorrplot(correlate(env),method = "spearman", type = "upper",diag = FALSE,cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(sig_thres = 0.05,only_mark= T,sig_level = c(0.05, 0.01, 0.001),mark = c("*", "**", "***"), size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),limits=c(-1,1),breaks=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)) + #颜色设置
  geom_couple(df_mantel, aes(color = df_p,
                             size = df_r),curvature = nice_curvature(by="from"))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("#4A9746","#645F88","#CCCCCC"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")


#环境因子相关性分析及展示：
env <-data.matrix(env1)
quickcor(env, type = "lower",method = "spearman") +
  geom_square()+
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')


df_mantel <- mantel_test(tasv, env,
                         mantel.fun = 'mantel',spec.dist.method = 'bray', 
                         env.dist.method = 'euclidean')

#定义标签：
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.25, Inf),
                    labels = c("< 0.1", "0.1 - 0.25", ">= 0.25")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf,0.01, 0.05, Inf),
                    labels = c("< 0.01","0.01<x<0.05",">= 0.05")))#定义Mantel的P值范围标



quickcor(env,method = "spearman", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = '#2166AC', mid = 'white',low = '#CC0000',midpoint = 0.25) + #颜色设置
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("#009100","#996699","#D95F02","#CCCCCC"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")



#基于kegg 可视化

BiocManager::install("clusterProfiler")
library(clusterProfiler)





#相关代谢途径
library(RColorBrewer)
#mcyc
setwd('~/Desktop/article1///宏基因组分析//基于read宏基因组注释/mcyc/')
location <- read.table("~/Desktop/article1///宏基因组分析//基于read宏基因组注释/Location")
data1 <- read.table('pathway-number-mcyc',header = T,row.names = 1,sep = '\t',check.names = F)
dd <-merge(t(data1),location,by = "row.names")
kk<-tidyr::gather(dd,pathway,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,pathway) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
color2 <- brewer.pal(10,"Set3")



hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")
pdf("G://work//宏基因组分析//基于read宏基因组注释//mcyc/mcyc.pdf",width = 10,height = 6)
ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = pathway),position="fill",stat = "identity",width = 1)+
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  scale_fill_manual(values = color2)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
dev.off()



#ncyc
setwd('G://work/宏基因组分析/基于read宏基因组注释/ncyc/')
location <- read.table("G:\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
data1 <- read.table('pathway_ncyc_number',header = T,row.names = 1,sep = '\t',check.names = F)
dd <-merge(t(data1),location,by = "row.names")
kk<-tidyr::gather(dd,pathway,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,pathway) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

color2 <- brewer.pal(9,"Set3")



hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")
pdf("G://work//宏基因组分析//基于read宏基因组注释//ncyc/ncyc.pdf",width = 10,height = 6)
ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = pathway),position="fill",stat = "identity",width = 1)+
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  scale_fill_manual(values = color2)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

dev.off()




#scyc

color2 <- brewer.pal(11,"Spectral")


setwd('G://work/宏基因组分析/基于read宏基因组注释/scyc/')
location <- read.table("G:\\work\\宏基因组分析\\基于read宏基因组注释\\Location")
data1 <- read.table('pathway-scys-number',header = T,row.names = 1,sep = '\t',check.names = F)
dd <-merge(t(data1),location,by = "row.names")
kk<-tidyr::gather(dd,pathway,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,pathway) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

color2 <- brewer.pal(11,"Spectral")



hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")



pdf("G://work//宏基因组分析//基于read宏基因组注释//scyc/scyc.pdf",width = 10,height = 6)

ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = pathway),position="fill",stat = "identity",width = 1)+
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  scale_fill_manual(values = color2)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
dev.off()




#pathway 统计图

#设置环境,构建scyc-pathway柱状堆积图
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释/scyc")
#读入数据
data1<-read.table("scyc-pathway-number",row.names = 1,header = T,check.names = F)

location <- read.table("/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释/Location")
dd <-merge(t(data1),location,by = "row.names")
kk<-tidyr::gather(dd,pathway,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,pathway) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

color2 <- brewer.pal(10,"Set3")



hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")

pdf("scyc.pdf",width = 10,height = 6)

ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = pathway),position="fill",stat = "identity",width = 1)+
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  scale_fill_manual(values = color2)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
dev.off()




#设置环境,构建mcyc-pathway柱状堆积图
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释/mcyc")
#读入数据
data1<-read.table("mcyc-pathway-number",row.names = 1,header = T,check.names = F)

location <- read.table("/Users/hahamark/Desktop/article1/宏基因组分析/基于read宏基因组注释/Location")
dd <-merge(t(data1),location,by = "row.names")
kk<-tidyr::gather(dd,pathway,Abundance,-c(Location,Row.names))
hh <- kk %>%group_by(Location,pathway) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))

color2 <- brewer.pal(11,"Spectral")



hh1<- c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
        "Qingdao","Dongying","Dandong")
id1=rep(1:12, each=1)

hh2 <- data.frame(hh1,id1)
colnames(hh2)<-(c("Location","id"))
hh3 <- merge(hh,hh2,by="Location")

pdf("scyc.pdf",width = 10,height = 6)

ggplot(hh3)+
  geom_bar(aes(x =id,y = Abundance,fill = pathway),position="fill",stat = "identity",width = 1)+
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  scale_fill_manual(values = color2)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_x_continuous(breaks=1:12,labels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
dev.off()




