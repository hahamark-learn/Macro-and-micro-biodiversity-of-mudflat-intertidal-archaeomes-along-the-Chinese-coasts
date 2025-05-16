#kegg_contig_绘图
#加载对应的包
library(tidyverse)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggpmisc)
library(RColorBrewer)
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/")
files <- list.files(pattern = "kegg_result")
combined_data <- data.frame()
for (i in files){
  data <- read.table(i,header = T,row.names = 1,check.names = F,sep = "\t")
  if (nrow(combined_data)==0){
    combined_data <- data
  }else{
    combined_data <- merge(combined_data,data,by = "row.names",all = T)
    rownames(combined_data) <- combined_data$Row.names
    combined_data <-combined_data[, -1]
  }
  
}



combined_data[is.na(combined_data)] <- 0
write.table(file = "all_kegg_count_reads.txt",combined_data,quote = F,sep = "\t")
#通过coverm mean 来排出异常值#k00059
kegg_data_integer<- read.table("all_kegg_count_reads_clean.txt",header = T,row.names = 1,check.names = F)

#进行抽平在进行分析
colSums(kegg_data_integer)
#使用该代码抽平
kegg_data_Flattening = as.data.frame(t(rrarefy(t(kegg_data_integer), min(colSums(kegg_data_integer)))))
rownames(kegg_data_Flattening) <- row.names(kegg_data_integer)
tkegg_data_Flattening <- t(kegg_data_Flattening)
write.table(file = "/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/all_contig_kegg_raw_flattening.txt",kegg_data_Flattening,quote = F,sep = "\t")
#进行抽平在进行分析
richness_kegg <- specnumber(kegg_data_Flattening,MARGIN = 2)
richness_kegg<- data.frame(richness_kegg)
colnames(richness_kegg) <- "richness"
location_kegg <- read.table('~/Desktop/article1/宏基因组分析/基于read宏基因组注释/location-ldg.txt',header = T,row.names = 1,check.names = F)
data_kegg_richness_location <- merge(richness_kegg,location_kegg,by = "row.names")
colnames(data_kegg_richness_location)[3] <-c("Latitude")
colnames(data_kegg_richness_location)[4] <-c("location")
data_kegg_richness_location$type <- "function"
#读入16s相关数据已经抽平了：
data_16s <- read.table("~/Desktop/article1/16s分析/16s-all/ASVs_counts_clean_Flattening-rdp-new.txt",header = T,row.names = 1,check.names = F)
richness_16s <- data.frame(specnumber(data_16s,MARGIN = 2))
tdata_16s_Flattening <- t(data_16s)
colnames(richness_16s) <- "richness"
location_16s <- read.table("~/Desktop/article1/16s分析/16s-all/location-ldg.txt",header = T,row.names = 1)
data_16s_richness_location <- merge(richness_16s,location_16s,by="row.names")
colnames(data_16s_richness_location)[3] <-c("Latitude")
colnames(data_16s_richness_location)[4] <-c("location")
data_16s_richness_location$type <- "taxonomy"
#整合数据合并16s和read数据
data_all_richness_location <- rbind(data_16s_richness_location,data_kegg_richness_location)
#赋予颜色
color1 <- c("#1f77b4","#f6b5d1")
p1 <- ggplot(data_all_richness_location,aes(Latitude,richness,color=type))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )+theme_bw()+
  scale_color_manual(values = color1)

ggsave("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/16s_kegg_ldg_2025.5.9.pdf",p1,width = 8,height = 6)
#计算ddr
bray_distance_tasv <- vegdist(tkegg_data_Flattening,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
write.table(data2,"/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/kegg_flattening_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#经过for_ddr，python脚本处理后读入
data_kegg_ddr <- read.table("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/kegg_flattening_bray_similirty_get_location.txt",header = T)

#计算16s_ddr的情况
bray_distance_tasv <- vegdist(tdata_16s_Flattening,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
write.table(data2,"~/Desktop/article1/宏基因组分析/contig_taxnonmy/cd_hit_result/kegg/16s_flattening_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#经过for_ddr，python脚本处理后读入
data_16s_ddr <- read.table("~/Desktop/article1/宏基因组分析/contig_taxnonmy/cd_hit_result/kegg/16s_flattening_bray_similirty_get_location.txt",header = T)
#合并数据
data_16s_ddr$type <- "taxonomy"
data_kegg_ddr$type <- "function"
data_all_ddr <- rbind(data_16s_ddr,data_kegg_ddr)
p2 <- ggplot(data_all_ddr,aes(x=log10(distance),y=log10(similarity),color=type))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+xlim(5.25,6.5)+
  geom_smooth(method = "lm",formula =y ~ x)+
  stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label..,p.value.label,sep ="~~~")),parse = TRUE)+
  theme_bw()+
  scale_color_manual(values = color1)
p2
ggsave("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/16s_kegg_ddr_2025.5.9.pdf",p2,width = 8,height = 6)


#pcoa 
kegg_data_Flattening <- read.table("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/all_contig_kegg_raw_flattening.txt",header = T,row.names = 1,check.names = F)

color1 <- brewer.pal(12,"Set3")
tdata1 <- t(kegg_data_Flattening)
distance <- vegdist(tdata1,method = "bray")
pcoa <- cmdscale(distance,k=(nrow(tdata1)-1),eig = T)
#提取前两个分类进行解释
plot_data <- data.frame({pcoa$point})[1:2]
#前两个分类解释命名
names(plot_data)[1:2] <- c('PCoA1','PCoA2')
eig=pcoa$eig
plot_data$id <- rownames(plot_data)
data <- merge(plot_data,location_kegg,by = "row.names")
data$Location <- factor(data$Location,levels= c("Sanya","Beihai","Zhuhai","Ningbo","Yancheng",
                                                   "Lianyungang","Shantou","Xiamen","Wenzhou",
                                                   "Qingdao","Dongying","Dandong"))
#基于bray-curtis距离进行计算
asv1.div <- adonis2(tdata1 ~ Location,data=location_kegg,permutations = 999,method = "bray")
#把统计结果放在pcoa上
asv1_adonis <- paste0("adonis R2: ",round(asv1.div$R2,2)," P-value:",asv1.div$`Pr(>F)`)

p3<- ggplot(data,aes(x=data$PCoA1,y=data$PCoA2,color=Location))+
  geom_point(alpha=1,size=3)+
  stat_ellipse(level = 0.95,size=1)+
  scale_color_manual(values = color1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""),title = asv1_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+geom_hline(aes(yintercept=0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white',colour = 'black'),axis.title.x=element_text(colour = 'black',size=20),axis.title.y = element_text(colour = 'black',size=20),legend.text = element_text(size = 15))
p3
ggsave("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/coverm_count/kegg/pcoa-kegg-25.5.9.pdf",p3,width = 8,height = 6)
