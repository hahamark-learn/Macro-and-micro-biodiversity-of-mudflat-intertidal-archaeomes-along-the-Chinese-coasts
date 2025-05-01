#contig_functional_pathway N-M-S KEGG
#加载相关的包
rm(list=ls())
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(ggcor)
#读入总体文件夹
setwd("~/Desktop/article1/宏基因组分析/contig_taxnonmy/")
#读相关数据
data1<- read.table('./cd_hit_result/ncyc/all_contig_ncyc.txt',header = T,row.names = 1,check.names = F)
data11<- read.table('./cd_hit_result/mcyc/all_contig_mcyc.txt',header = T,row.names = 1,check.names = F)
data111 <- read.table('./cd_hit_result/scyc/all_contig_scyc.txt',header = T,row.names = 1,check.names = F)
data1111 <- read.table('./cd_hit_result/kegg/all_contig_kegg.txt',header = T,row.names = 1,check.names = F)
metadata <- readxl::read_xlsx('~/Desktop/article1/宏基因组分析/基于read宏基因组注释/sample.xlsx')
location1 <- read.table('Location',header = T,row.names = 1)


#计算richness 
Richness <- specnumber(data1,MARGIN = 2)
data2 <- data.frame(Richness)
Richness1 <- specnumber(data11,MARGIN = 2)
data22 <- data.frame(Richness1)
Richness11 <- specnumber(data111,MARGIN = 2)
data222 <- data.frame(Richness11)
Richness111 <- specnumber(data1111,MARGIN = 2)
data2222 <- data.frame(Richness111)
data_all <- cbind(data2,data22,data222,data2222)
color1 <- brewer.pal(12,"Set3")


#计算alpha多样性
df3 <- merge(location1,data_all,by="row.names",all = T)
colnames(df3)[3:6] <-c("Ncyc","Mcyc","Scyc","kegg")
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

pdf("richness-ncyc-25.2.27.pdf",width = 10,height = 6)
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

pdf("richness-mcyc-25.2.27.pdf",width = 10,height = 6)
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

pdf("richness-scyc-25.2.27.pdf",width = 10,height = 6)
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

pdf("richness-kegg-25.2.27.pdf",width = 10,height = 6)
ggplot(df3,aes(Location,kegg,fill=Location))+
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


#Beta 多样性
#ncyc beta多样性
#B多样性分析,使用上面抽平过的数据
tdata1 <- t(data1)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(tdata1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location1,by = "row.names")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(tdata1 ~ Location,data=location1,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)
pdf("pcoa-ncyc-25.2.27.pdf",width = 10,height = 6)
ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))
dev.off()

#mcyc beta多样性
#B多样性分析,使用上面抽平过的数据
tdata1 <- t(data11)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(tdata1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location1,by = "row.names")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(tdata1 ~ Location,data=location1,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)
pdf("pcoa-mcyc-25.2.27.pdf",width = 10,height = 6)
ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))
dev.off()


#scyc beta多样性
#B多样性分析,使用上面抽平过的数据
tdata1 <- t(data111)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(tdata1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location1,by = "row.names")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(tdata1 ~ Location,data=location1,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)
pdf("pcoa-scyc-25.2.27.pdf",width = 10,height = 6)
ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))
dev.off()

#kegg beta多样性
#B多样性分析,使用上面抽平过的数据
tdata1 <- t(data1111)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(tdata1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location1,by = "row.names")
#基于bray-curtis距离进行计算
asv1.div <- adonis2(tdata1 ~ Location,data=location1,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)
pdf("pcoa-kegg-25.2.27.pdf",width = 10,height = 6)
ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))
dev.off()

#functional LDG 模式绘画
#LDG
#读入统计的数据，N,S,CH4
library(ggpubr)
library(ggpmisc)
library(reshape2)
location2 <- read.table('location-ldg.txt',header = T,row.names = 1,check.names = F)
data_all_ldg <- merge(data_all,location2,by = "row.names")
colnames(data_all_ldg)[2:7] <-c("Ncyc","Mcyc","Scyc","kegg","Latitude","Location")
data_all_ldg <- data_all_ldg[,-1]

data_new <- melt(data_all_ldg,id.vars = c("Latitude","Location"))
#这里使用stat_cor 没有定义。是计算了相关性，并不是拟合曲线的值
upper <- ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )+theme_classic()+
  ylab("")+
  xlab("")+
  coord_cartesian(ylim=c(1000,3000))+
  scale_y_continuous(breaks = c(1000, 3000, 500)) + # 以250为单位划分Y轴
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank())

down <- ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )+theme_classic()+
  ylab("Richness")+
  coord_cartesian(ylim=c(0,200))

p1 <- ggarrange(upper,down,ncol=1,align = "hv",common.legend = T)

ggsave("NCS_KEGG_ldg_25.3.13.pdf",p1,width = 8,height = 6)










#从新绘制ldg的图，并放在放在一起比较

p1<- ggplot(data_new,aes(Latitude,value,color=variable))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )+
  theme_classic()


ggsave("NCS_KEGG_ldg_20-25。3.13.pdf",p1,width = 8,height = 6)




#拟合曲线采用，容易受到其他因素的干扰。

data_kegg <- data_all_ldg[,3:5]

# 拟合鲁棒线性回归模型
robust_model <- rlm(kegg ~ Latitude, data = data_kegg)

# 生成预测值
data_kegg$y_pred <- predict(robust_model)

ggplot(data_kegg, aes(x = Latitude, y = kegg)) +
  geom_point(size = 3, color = "blue") +  # 绘制散点图
  geom_line(aes(y = y_pred), color = "red", size = 1) +  # 绘制鲁棒回归拟合曲线
  geom_smooth(method = "lm",formula =y ~ x)+
  labs(
    title = "Robust Linear Regression Fit",
    x = "X",
    y = "Y"
  ) +
  theme_minimal()+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )



#DDR
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
bray_distance_data_all <- vegdist(t(data1111),method = "bray")
data_kegg <- as.matrix(bray_distance_data_all,row.names=1)
data_kegg[upper.tri(data_kegg)]=NA
data_kegg <- melt(data_kegg,na.rm = T)
data_kegg$group <- "kegg"






#合并4个模块的因子
data2 <- rbind(data_scyc,data_mcyc,data_ncyc,data_kegg)
write.table(data2,"ncs-kegg-bray_similirty-25.2.27.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#计算distance
location3 <- read.table('location-ddr',header = T,row.names = 1,check.names = F)
distance1 <-distm(location3,fun = distCosine)
rownames(distance1)<- rownames(location3)
colnames(distance1) <- rownames(location3)
data_distance1 <- as.matrix(distance1,row.names=1)
data_distance1[upper.tri(data_distance1)]=NA
data_distance2 <- melt(data_distance1,na.rm = T)
write.table(data_distance2,"bray_similirty_distance1-25.2.27.txt",quote = FALSE, row.names = FALSE,sep = "\t")

#经过python脚本处理后，绘制图片
data1 <- read.table("ncs-kegg-bray_similarity_distance-25.2.27.txt",header = T,check.names = F)
library(ggplot2)
p1<- ggplot(data1,aes(x=log10(distance),y=log10(similarity),color=group))+
  geom_point(size=1)+
  geom_smooth(method = "lm",linewidth=1.5)+theme_classic()+xlim(5.25,6.5)+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )+
  scale_color_manual(values = c("#A07DB6","#78A72B","#EB736A","#1EB4B8"))+
  labs(x="Geographical distance (log10 scale)",y="Similarity (log10 scale)")

ggsave("nms_kegg_ddr_2025.3.13.pdf",p1,width = 8,height = 6)



#mantel-test分析：


#gene family  mantel test
library(linkET)
#读入总体文件夹
setwd("~/Desktop/article1/宏基因组分析/contig_taxnonmy/")
#读相关数据
data1<- read.table('./cd_hit_result/ncyc/all_contig_ncyc.txt',header = T,row.names = 1,check.names = F)
data11<- read.table('./cd_hit_result/mcyc/all_contig_mcyc.txt',header = T,row.names = 1,check.names = F)
data111 <- read.table('./cd_hit_result/scyc/all_contig_scyc.txt',header = T,row.names = 1,check.names = F)
data1111 <- read.table('./cd_hit_result/kegg/all_contig_kegg.txt',header = T,row.names = 1,check.names = F)
ncyc_function <- data.frame(t(data1))
mcyc_function <- data.frame(t(data11))
scyc_function <- data.frame(t(data111))
kegg_function <- data.frame(t(data1111))


#合并所有的参所
all_function <- cbind(ncyc_function,mcyc_function,scyc_function,kegg_function)
#读取环境因子
env_factor <-read.table("env_function.txt",header = T, row.names = 1, check.names = F)

#生成对应的table表格
mantel <- mantel_test(all_function, env_factor, mantel_fun = "mantel",
                      spec_select = list(
                        NCYC=1:50,MCYC=51:278,SCYC=279:463,KEGG=464:7973
                      ))%>% mutate(r_value = cut(r,breaks=c(-Inf,0.1,0.25,Inf),labels=c('<0.1','0.1-0.25','>=0.25'),right=F),
                                   p_value = cut(p,breaks=c(-Inf,0.001,0.01,0.05,Inf),labels=c('<=0.001','0.001-0.01','0.01-0.05','>=0.05'),right=F))






#最终绘图
color1 <- colorRampPalette(brewer.pal(11,"RdBu"))(50)



#绘图
p1<- qcorrplot(correlate(env_factor),method = "spearman", type = "upper",diag = FALSE,cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(sig_thres = 0.05,only_mark= T,sig_level = c(0.05, 0.01, 0.001),mark = c("*", "**", "***"), size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),limits=c(-1,1),breaks=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)) + #颜色设置
  geom_couple(mantel, aes(color = p_value,
                          size = r_value),curvature = nice_curvature(by="from"))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("#645F88","#4A9746","#CCCCCC"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")

ggsave("nms_kegg_mantel_2025.2.27.pdf",p1,width = 8,height = 6)
write.table(mantel,"mantel_function_contig-2025.2.27.txt",quote = F,sep = "\t")



#随机森林预测


