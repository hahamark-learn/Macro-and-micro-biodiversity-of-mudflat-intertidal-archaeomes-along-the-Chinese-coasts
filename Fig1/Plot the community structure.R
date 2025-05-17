#16s柱状堆积图
library(tidyverse)
library(ggprism)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(ggtree)
library(ggplot2)
library(dplyr)
library(reshape2)
library(linkET)
library(vegan)
library(RColorBrewer)
#读入数据，进行抽平处理
setwd("~/Desktop/article1/16s分析/16s-all/")
asv<-read.table("ASVs_counts_clean-rdp.txt",row.names = 1,header = T)
taxa<- read.table("ASVs_taxonomy-clean-rdp.txt",row.names = 1,header = T)
metadata <- readxl::read_xlsx("sample_df.xlsx")
rownames(metadata)<-metadata$id
#抽平处理
#数据抽平
colSums(asv)
#使用该代码抽平
asv_Flattening = as.data.frame(t(rrarefy(t(asv), min(colSums(asv)))))

tasv_Flattening <- t(asv_Flattening)
#查看抽平后的每个样本的和
colSums(asv_Flattening)
#将抽平后的otu表保存到该工作目录，准备后面的多样性分析
asv_Flattening
write.table(asv_Flattening,"ASVs_counts_clean_Flattening-rdp.txt",quote=F,sep = "\t")

#重新读入排除总和为0的那一行
asv_clean <- read.table("ASVs_counts_clean_Flattening-rdp-new.txt",row.names = 1,header = T)

#这里比较关键，此时asv为抽凭过的asv
asv <- asv_clean
dat <- merge(x=asv,y=taxa,by='row.names')
dat=dplyr::rename(dat,OTUID=Row.names)
aa<-aggregate(dat[,2:145],by=list(dat$Phylum),FUN=sum)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
#删除为0的那一行
aa <- aa[-5,]
head(aa, n = 3)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,1:ncol(aa)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-aa[order$ix,]
head(cc, n = 3)

##只展示排名前10的物种，之后的算作Others(根据需求改数字)
#dd<-rbind(colSums(cc[11:as.numeric(length(rownames(cc))),]),cc[10:1,])
#head(dd, n = 3)
#rownames(dd)[1]<-"Others"
#head(dd, n = 3)
bb<-merge(t(cc),dplyr::select(metadata,id,Location),
          by.x = "row.names",by.y ="id")
kk<-tidyr::gather(bb,Phylum,Abundance,-c(Location,Row.names))
kk$Phylum<-ifelse(kk$Phylum=='Unassigned','Others',kk$Phylum)
hh <- kk %>%group_by(Location,Phylum) %>%dplyr :: summarise(Abundance=sum(Abundance))
hh$Phylum <- factor(hh$Phylum,order = T,levels =rev(c("Thaumarchaeota","Euryarchaeota","Crenarchaeota","Woesearchaeota","Aenigmarchaeota","Nanohaloarchaeota","Diapherotrites","unclassified")))


color1 <- brewer.pal(12,"Set3")
color1 <- rev(c("#BC80BD","#CCEBC5","#B3DE69","#BEBADA","#8DD3C7","#DADF00","#FFCC33","#FFFFB3"))


hh$Phylum = factor(hh$Phylum,order = T,levels = rev(row.names(cc1)))
#hh$Location <-factor(hh$Location,levels= c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng'))
#c('Lianyungang','Zhuhai','Sanya','Dandong','Dongying','Qingdao','Xiamen','Shantou','Ningbo','Wenzhou','Beihai','Yancheng')
hh$Location <-factor(hh$Location,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
ggplot(hh)+
  geom_bar(aes(x =Location,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.93)+
  scale_fill_manual(values = color1) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 90))

ggplot(hh)+
  geom_bar(aes(x =Location,y = Abundance,fill = Phylum),position="fill",stat = "identity",width = 0.93)+
  scale_fill_manual(values = color1) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 65,colour = "black",family = "Times",size =15))



#绘制144个位置的图

hh1 <- kk%>% 
  group_by(Row.names) %>%
  mutate(precent = Abundance/sum(Abundance)*100)


hh1$Phylum <- factor(hh1$Phylum,order = T,levels =rev(c("Thaumarchaeota","Euryarchaeota","Crenarchaeota","Woesearchaeota","Aenigmarchaeota","Nanohaloarchaeota","Diapherotrites","unclassified")))

color1 <- rev(c("#BC80BD","#CCEBC5","#B3DE69","#BEBADA","#8DD3C7","#DADF00"))
hh1$precent <- as.numeric(hh1$precent)
hh1$Row.names <- as.factor(hh1$Row.names)
hh1$Row.names <- factor(hh1$Row.names,levels = c("FDMP22H003646_1a_L1_SY1","FDMP22H003646_1a_L1_SY11","FDMP22H003646_1a_L1_SY12","FDMP22H003646_1a_L1_SY13",
                                                 "FDMP22H003646_1a_L1_SY14","FDMP22H003646_1a_L1_SY15","FDMP22H003646_1a_L1_SY2","FDMP22H003646_1a_L1_SY3",
                                                 "FDMP22H003646_1a_L1_SY5","FDMP22H003646_1a_L1_SY6","FDMP22H003646_1a_L1_SY7","FDMP22H003646_1a_L1_SY9",
                                                 "FDMP22H003646_1a_L1_BH1","FDMP22H003646_1a_L1_BH10","FDMP22H003646_1a_L1_BH11","FDMP22H003646_1a_L1_BH12",
                                                 "FDMP22H003646_1a_L1_BH13","FDMP22H003646_1a_L1_BH14","FDMP22H003646_1a_L1_BH15","FDMP22H003646_1a_L1_BH2",
                                                 "FDMP22H003646_1a_L1_BH3","FDMP22H003646_1a_L1_BH4","FDMP22H003646_1a_L1_BH6","FDMP22H003646_1a_L1_BH8",
                                                 "FDMP22H003645_1a_L1_ZH1","FDMP22H003645_1a_L1_ZH11","FDMP22H003645_1a_L1_ZH12","FDMP22H003645_1a_L1_ZH13",
                                                 "FDMP22H003645_1a_L1_ZH14","FDMP22H003645_1a_L1_ZH15","FDMP22H003645_1a_L1_ZH2","FDMP22H003645_1a_L1_ZH3",
                                                 "FDMP22H003645_1a_L1_ZH4","FDMP22H003645_1a_L1_ZH5","FDMP22H003645_1a_L1_ZH7","FDMP22H003645_1a_L1_ZH9",
                                                 "FDMP22H003647_1a_L1_ST1","FDMP22H003647_1a_L1_ST10","FDMP22H003647_1a_L1_ST11","FDMP22H003647_1a_L1_ST12",
                                                 "FDMP22H003647_1a_L1_ST13","FDMP22H003647_1a_L1_ST14","FDMP22H003647_1a_L1_ST15","FDMP22H003647_1a_L1_ST2",
                                                 "FDMP22H003647_1a_L1_ST3","FDMP22H003647_1a_L1_ST5","FDMP22H003647_1a_L1_ST7","FDMP22H003647_1a_L1_ST9",
                                                 "FDMP22H003647_1a_L1_XM1","FDMP22H003647_1a_L1_XM10","FDMP22H003647_1a_L1_XM11","FDMP22H003647_1a_L1_XM12",
                                                 "FDMP22H003647_1a_L1_XM13","FDMP22H003647_1a_L1_XM14","FDMP22H003647_1a_L1_XM15","FDMP22H003647_1a_L1_XM2",
                                                 "FDMP22H003647_1a_L1_XM3","FDMP22H003647_1a_L1_XM4","FDMP22H003647_1a_L1_XM6","FDMP22H003647_1a_L1_XM8",
                                                 "FDMP22H003647_1a_L1_WZ1","FDMP22H003647_1a_L1_WZ10","FDMP22H003647_1a_L1_WZ11","FDMP22H003647_1a_L1_WZ13",
                                                 "FDMP22H003647_1a_L1_WZ14","FDMP22H003647_1a_L1_WZ15","FDMP22H003647_1a_L1_WZ2","FDMP22H003647_1a_L1_WZ3",
                                                 "FDMP22H003647_1a_L1_WZ4","FDMP22H003647_1a_L1_WZ5","FDMP22H003647_1a_L1_WZ7","FDMP22H003647_1a_L1_WZ9",
                                                 "FDMP22H003646_1a_L1_NB1","FDMP22H003646_1a_L1_NB11","FDMP22H003646_1a_L1_NB12","FDMP22H003646_1a_L1_NB14",
                                                 "FDMP22H003646_1a_L1_NB15","FDMP22H003646_1a_L1_NB2","FDMP22H003646_1a_L1_NB3","FDMP22H003646_1a_L1_NB4",
                                                 "FDMP22H003646_1a_L1_NB5","FDMP22H003646_1a_L1_NB6","FDMP22H003646_1a_L1_NB7","FDMP22H003646_1a_L1_NB9",
                                                 "FDMP22H003647_1a_L1_YC1","FDMP22H003647_1a_L1_YC12","FDMP22H003647_1a_L1_YC13","FDMP22H003647_1a_L1_YC14",
                                                 "FDMP22H003647_1a_L1_YC15","FDMP22H003647_1a_L1_YC2","FDMP22H003647_1a_L1_YC3","FDMP22H003647_1a_L1_YC4",
                                                 "FDMP22H003647_1a_L1_YC5","FDMP22H003647_1a_L1_YC7","FDMP22H003647_1a_L1_YC8","FDMP22H003647_1a_L1_YC9",
                                                 "FDMP22H003646_1a_L1_LYG1","FDMP22H003646_1a_L1_LYG11","FDMP22H003646_1a_L1_LYG12","FDMP22H003646_1a_L1_LYG13",
                                                 "FDMP22H003646_1a_L1_LYG14","FDMP22H003646_1a_L1_LYG15","FDMP22H003646_1a_L1_LYG2","FDMP22H003646_1a_L1_LYG3",
                                                 "FDMP22H003646_1a_L1_LYG4","FDMP22H003646_1a_L1_LYG5","FDMP22H003646_1a_L1_LYG6","FDMP22H003646_1a_L1_LYG7",
                                                 "FDMP22H003451_1a_L1_CXD_1412","FDMP22H003451_1a_L1_CXD_1413","FDMP22H003451_1a_L1_CXD_1414",
                                                 "FDMP22H003451_1a_L1_CXD_1415","FDMP22H003451_1a_L1_CXD_1604","FDMP22H003451_1a_L1_CXD_1608",
                                                 "FDMP22H003451_1a_L1_CXD_1609","FDMP22H003451_1a_L1_CXD_1610","FDMP22H003451_1a_L1_CXD_1611",
                                                 "FDMP22H003451_1a_L1_CXD_1612","FDMP22H003451_1a_L1_CXD_1613","FDMP22H003451_1a_L1_CXD_1614",
                                                 "FDMP22H003645_1a_L1_DY1","FDMP22H003645_1a_L1_DY10","FDMP22H003645_1a_L1_DY11","FDMP22H003645_1a_L1_DY13",
                                                 "FDMP22H003645_1a_L1_DY14","FDMP22H003645_1a_L1_DY15","FDMP22H003645_1a_L1_DY2","FDMP22H003645_1a_L1_DY3",
                                                 "FDMP22H003645_1a_L1_DY4","FDMP22H003645_1a_L1_DY6","FDMP22H003645_1a_L1_DY7","FDMP22H003645_1a_L1_DY9",
                                                 "FDMP22H003645_1a_L1_DD1","FDMP22H003645_1a_L1_DD10","FDMP22H003645_1a_L1_DD11","FDMP22H003645_1a_L1_DD13",
                                                 "FDMP22H003645_1a_L1_DD14","FDMP22H003645_1a_L1_DD15","FDMP22H003645_1a_L1_DD2","FDMP22H003645_1a_L1_DD3",
                                                 "FDMP22H003645_1a_L1_DD4","FDMP22H003645_1a_L1_DD6","FDMP22H003645_1a_L1_DD7","FDMP22H003645_1a_L1_DD9"))


ggplot(hh1,aes(x =hh1$Row.names,y = precent,fill = Phylum))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  scale_fill_manual(values = color1) +
  labs(x='Group',y='Abundance(%)')+guides(fill=guide_legend(reverse = TRUE))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 90))

ggsave("fig_16s_rdp.pdf", width = 8.27, height = 11.69)



#alpha多样性
metadata1 <- metadata[,-1]
row.names(metadata1) <- metadata$id
#shannon 指数
Shannon<-diversity(asv, index = "shannon", MARGIN = 2, base = exp(1))
Shannon <- data.frame(Shannon)
id<- data.frame((row.names(Shannon)))
colnames(id)<-"id"
shannon1  <- cbind(id,Shannon)
shannon2 <- merge(shannon1,metadata,by="id",all=T)
shannon2$Location <- factor(shannon2$Location, levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",                                                         "Qingdao","Dongying","Dandong"))
ggplot(shannon2,aes(Location,Shannon,fill=Location))+
  geom_violin()+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)


#除去统计的行数
asv <- asv_clean[,1:144]

#Simpson 指数
Simpson<-diversity(asv, index = "simpson", MARGIN = 2, base = exp(1))
#Richness 指数
Richness <- specnumber(asv,MARGIN = 2)
richness <- data.frame(Richness)
richness$id<- row.names(richness)
richness2 <- merge(richness,metadata,by="id",all=T)

color1 <- brewer.pal(12,"Paired")
#地图地点纬度排行
richness2$Location <- factor(richness2$Location, levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                                          "Qingdao","Dongying","Dandong"))
ggplot(richness2,aes(Location,Richness,fill=Location))+
  geom_violin()+
  geom_boxplot(position = position_identity(),width=0.05,fill="white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1,colour = "black",family = "Times",size =10),
        axis.text.y = element_text(family = "Times",size = 10,colour = "black"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color1)



#Beta多样性
#B多样性分析
asv1 <- t(asv)
distance <- vegdist(asv1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(asv1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
group1 <-metadata['Location']
color1 <- brewer.pal(12,"Set3")
data <- merge(plot_data,metadata,by = "id")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(asv1 ~ Location,data=metadata,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)



p1<- ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))

ggsave(filename = "16s_pcoa.pdf",plot=p1,width=10,height =8 )

#LDG
#有richness 和纬度坐标,手动调整，参考本文件夹下的richness-ldg-cankao
library(ggpmisc)
library(ggpubr)
Richness <- specnumber(asv,MARGIN = 2)
richness <- data.frame(Richness)
write.table(richness,"richness",sep="\t", quote=F)
richness_all <-read.table("Richness-ldg",sep = "\t",header = T,row.names = 1)
color1 <- brewer.pal(12,"Paired")
my.formula <- y ~ x


richness_all <-read.table("richness-1",sep = "\t",header = T,row.names = 1)
ggplot(richness_all,aes(x=Latitude,y=Richness))+
  geom_point(aes(colour=richness$group))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2),color="black")+
  scale_color_manual(values = color1)+theme_classic()





p2 <- ggplot(richness_all,aes(x=Latitude,y=Richness))+
  geom_point(color="#8DEEEE")+
  geom_smooth(method = "lm",color="black")+
  scale_color_manual()+ 
  stat_cor(aes())+
  theme_classic()






p1<-ggplot(richness_all,aes(x=Latitude,y=Richness))+
  geom_point(color="#8DEEEE")+
  geom_smooth(method = "lm",color="black")+
  scale_color_manual()+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label..,p.value.label,sep ="~~~")),
               parse = TRUE) +
  theme_classic()


ggsave("16s_ldg_2024.3.21.pdf",p1,width = 8,height = 6)

ggsave("16s_ldg_2024.3.24.pdf",p2,width = 8,height = 6)

#DDR
library(reshape2)
library(geosphere)
#计算similirty，使用vegdist 计算的相异度而不是相似性，所以下方数据合并的时候要用1-

tasv <- t(asv)
bray_distance_tasv <- vegdist(tasv,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
data(varespec)
write.table(data2,"asv_clean_bray_disimilirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#计算distance
location_data <- readxl::read_xlsx("location.xlsx")
location_name <- location_data$station
rownames(location_data) <- location_data$station
location_new_data <- location_data[,2:3]
rownames(location_new_data) <- location_data$station
distance1 <-distm(location_new_data,fun = distCosine)
rownames(distance1)<- location_data$station
colnames(distance1) <- location_data$station
data_distance1 <- as.matrix(distance1,row.names=1)
data_distance1[upper.tri(data_distance1)]=NA
data_distance2 <- melt(data_distance1,na.rm = T)
write.table(data_distance2,"bray_similirty_distance1.txt",quote = FALSE, row.names = FALSE,sep = "\t")

#作图:这里已经用1减过了，不用再减
data1 <- read.table("asv_clean_similarity_distance.txt",header = T)
library(ggplot2)
p1<-ggplot(data1,aes(x=log10(distance),y=log10(similarity)))+
  geom_point(color="#8DEEEE")+
  geom_smooth(method = "lm",color="black")+theme_classic()+xlim(5.25,6.5)+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label..,p.value.label,sep ="~~~")),
               parse = TRUE) +
  labs(x="Geographical distance (log10 scale)",y="Similarity (log10 scale)")


p2<-ggplot(data1,aes(x=log10(distance),y=log10(similarity)))+
  geom_point(color="#8DEEEE")+
  geom_smooth(method = "lm",color="black")+theme_classic()+xlim(5.25,6.5)+
  stat_cor(aes())+
  labs(x="Geographical distance (log10 scale)",y="Similarity (log10 scale)")




ggsave("16s_ddr_2024.3.21.pdf",p1,width = 8,height = 6)
ggsave("16s_ddr_2024.3.24.pdf",p2,width = 8,height = 6)
#mantel test
#加载包
library(vegan)
library(dplyr)
library(ggplot2)
library(linkET)
#读入数据
#asv<- read.table("ASVs_counts_clean-rdp.txt",sep = '\t')
#taxa1<- read.table("ASVs_taxonomy-clean-rdp.txt", sep="\t")
group<- readxl::read_xlsx("sample_df.xlsx")
env <- read.table("environmental_factors.txt",header = T,row.names = 1,sep="\t")
#数据变形，这里的asv应该为抽频过的asv
tasv <- t(asv)
asv1 <- tasv[match(rownames(env),rownames(tasv)),]



#环境因子相关性分析及展示：
env <-data.matrix(env)
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
  scale_color_manual(values = c("#4A9746","#645F88","#A6AD40","#CCCCCC"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")


#high = 'orange', mid = 'white',low = 'navyblue'





#优势种相关性分析：
#mantel test
#加载包
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
#读入数据
setwd("~/Desktop/article1/16s分析//16s-all")
asv<- read.table("ASVs_counts_clean-rdp.txt",sep = '\t')
taxa1<- read.table("ASVs_taxonomy-clean-rdp.txt", sep="\t")
group<- readxl::read_xlsx("sample_df.xlsx")
env <- read.table("environmental_factors.txt",header = T,row.names = 1,sep="\t")

#先统计处排名前十的优势种，进行相关性的分析
rownames(group) <- group$id 
dat <- merge(x=asv,y=taxa1,by='row.names')
dat=dplyr::rename(dat,OTUID=Row.names)
num_a <- ncol(asv)+1
aa<-aggregate(dat[,2:num_a],by=list(dat$Phylum),FUN=sum)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
head(aa, n = 3)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,2:ncol(aa)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-aa[order$ix,]

#数据变形
env1 <- env[,1]
env1<- data.frame(env1)
rownames(env1) <- row.names(env)
cc1<-cc[-7,]


tcc1 <- t(cc1)
asv1 <- tcc1[match(rownames(env1),rownames(tcc1)),]




#环境因子相关性分析及展示：
env <-data.matrix(env)
quickcor(env, type = "lower",method = "spearman") +
  geom_square()+
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')





df_mantel <- mantel_test(asv1, env1,
                         spec.select = list(Thaumarchaeota=1,
                                            Euryarchaeota=2
                         ),
                         mantel.fun = 'mantel',spec.dist.method = 'bray', 
                         env.dist.method = 'euclidean')







#定义标签：
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.25, Inf),
                    labels = c("< 0.1", "0.1 - 0.25", ">= 0.25")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf,0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01<x<0.05",">= 0.05")))#定义Mantel的P值范围标





quickcor(env,method = "spearman", type = "upper", cor.test = T, cluster.type = "all",show.diag = FALSE) +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = "#08306B", mid = 'white',low = "#A50026",limits=c(-1,1),breaks=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)) + #颜色设置
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("#009100","#996699","#D95F02","#CCCCCC"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")


#high = 'orange', mid = 'white',low = 'navyblue'


#增加PD相关计算
install.packages('piacnte')

library(picante)

data(phylocom)
names(phylocom)

phy <- phylocom$phylo # "phylo" class
comm <- phylocom$sample # "matrix" class
pd.result <- pd(comm, phy, include.root = TRUE)



#16S 计算PD,需要进化树的支持
library(iCAMP)
library(ape)
library(vegan)

range(rowSums(asv))
comm=vegan::rrarefy(asv,sample = min(rowSums(asv)))
comm=comm[,colSums(comm)>0,drop=FALSE]
phy_tree <- read.tree(file = "16s.tree")
pd.result <- pd(commin, tree, include.root = TRUE)






#读取icamp 相关数据
icis<-load("C://Users/DELL/Desktop/work/16s分析/raw-data/Test.iCAMP.SES.RC.detail.rda")
icamp_out <- res$bNRIiRCa
#保存icamp 结果

write.table(icamp_out,"C://Users/DELL/Desktop/work/16s分析/raw-data/icamp_out",sep = "\t",quote	=F)


#计算后，进行再次修改
library(ggpubr)
#总体的值的数量。
mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0.027471548,0.541221775,0.194135052,0.01740642,0.219765205)
)


#确定后的值
#0.027358204	0.541429746	0.195770125	0.017249856	0.218192069
mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0,
          0,
          0.240574008,
          0.074883868,
          0.684542124
  )
)

mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0,
          0.867235143,
          0.046369821,
          0.003135258,
          0.083259778
  )
)

icamp_all <- ggdonutchart(data=mydata, 
             x="Value",  
             label = "Type",  
             fill="Type",   
             lab.pos="in",
             palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)


ggsave("~/Desktop/article1/fig5/icamp_all.pdf",icamp_all,width=4,height = 5 )

#优势门相关的bin的icamp 值

#bin152

mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0,0.621822334,0.047948852,0.030253892,0.299974922)
)

icamp_152 <- ggdonutchart(data=mydata, 
                          x="Value",  
                          label = "Type",  
                          fill="Type",   
                          lab.pos="in",
                          palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)

ggsave("~/Desktop/article1/fig5/icamp_152.pdf",icamp_152,width=4,height = 5 )

#bin154

mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0,0.783590412,0.069071784,0.013846419,0.133491385)
)

icamp_154 <- ggdonutchart(data=mydata, 
                          x="Value",  
                          label = "Type",  
                          fill="Type",   
                          lab.pos="in",
                          palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)

ggsave("~/Desktop/article1/fig5/icamp_154.pdf",icamp_154,width=4,height = 5 )

#bin94

mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0.325210724,0.411561555,0.187534082,0.0000564661584656553,0.0756371735)
)

icamp_94 <- ggdonutchart(data=mydata, 
                          x="Value",  
                          label = "Type",  
                          fill="Type",   
                          lab.pos="in",
                          palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)

ggsave("~/Desktop/article1/fig5/icamp_94.pdf",icamp_94,width=4,height = 5 )


#bin11

mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0,0.111074718,0.409465771,0.032454358,0.447005153)
)

icamp_11 <- ggdonutchart(data=mydata, 
                         x="Value",  
                         label = "Type",  
                         fill="Type",   
                         lab.pos="in",
                         palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)

ggsave("~/Desktop/article1/fig5/icamp_11.pdf",icamp_11,width=4,height = 5 )

#bin6
mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0,0.187936431,0.317612602,0.020306167,0.474144801)
)

icamp_6 <- ggdonutchart(data=mydata, 
                         x="Value",  
                         label = "Type",  
                         fill="Type",   
                         lab.pos="in",
                         palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)

ggsave("~/Desktop/article1/fig5/icamp_6.pdf",icamp_6,width=4,height = 5 )


#bin5
mydata<-data.frame(
  Type=c("Hes","Hos","DL","HD","DR"),
  Value=c(0.000709092,0.325059261,0.210531288,0.062457265,0.401243094)
)

icamp_5 <- ggdonutchart(data=mydata, 
                        x="Value",  
                        label = "Type",  
                        fill="Type",   
                        lab.pos="in",
                        palette = c("#EE7671","#9C9D3B","#37B883","#28B7DD","#AB76B3")
)

ggsave("~/Desktop/article1/fig5/icamp_5.pdf",icamp_5,width=4,height = 5 )









#读取bin的结果
icis_bin<-load("C://Users/DELL/Desktop/work/16s分析/raw-data/Test.iCAMP.Summary.rda")
library(iCAMP)
library(ape)
setwd("C://Users/DELL/Desktop/work/16s分析/raw-data/")

# the OTU table file (Tab delimited txt file)
com.file="16s-otu.csv"

# the phylogenetic tree file
tree.file="16s.tree"

# the classification (taxonomy) information
clas.file="16s-class.csv"

# the treatment informaiton table
treat.file="site-group-1.csv"
rand.time=1000

commin=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                    as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                    check.names = FALSE))
tree=read.tree(file = tree.file)
clas=read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
treat=read.table(treat.file, header = TRUE, sep = ",", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)



icbin=iCAMP::icamp.bins(icamp.detail = res$detail,treat = treat,
                        clas=clas,silent=FALSE, boot = TRUE,
                        rand.time = rand.time,between.group = F)

prefix="Test1"


palette.pals('Set3',9)



#在优势的门里找到优势的种
#找到Thaumarchaeota 优势的目

dat_Thaumarchaeota <- subset(dat,dat$Phylum == "Thaumarchaeota")

dat_Thaumarchaeota <- dplyr::rename(dat_Thaumarchaeota,OTUID=Row.names)

aa<-aggregate(dat_Thaumarchaeota[,2:145],by=list(dat_Thaumarchaeota$Class),FUN=sum)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)


#找到Euryarchaeota 优势的目
dat_Euryarchaeota <- subset(dat,dat$Phylum == "Euryarchaeota")

dat_Euryarchaeota <- dplyr::rename(dat_Euryarchaeota,OTUID=Row.names)

aa<-aggregate(dat_Euryarchaeota[,2:145],by=list(dat_Euryarchaeota$Class),FUN=sum)

row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
head(aa, n = 3)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,1:ncol(aa)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-aa[order$ix,]


#找到 Crenarchaeota 优势的目

dat_Crenarchaeota <- subset(dat,dat$Phylum == "Crenarchaeota")

dat_Crenarchaeota <- dplyr::rename(dat_Crenarchaeota,OTUID=Row.names)

aa<-aggregate(dat_Crenarchaeota[,2:145],by=list(dat_Crenarchaeota$Class),FUN=sum)

row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
head(aa, n = 3)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,2:ncol(aa)]),index.return=TRUE,decreasing=T)   
#根据列求和结果对表格排序
cc<-aa[order$ix,]



#NST 分析
install.packages("NST")
