library(ggplot2)
library(RColorBrewer)
library(vegan)
library(ggpubr)
library(ggpmisc)
#微观多样性
setwd("~/Desktop/article1/宏基因组分析/基于contig微观多样性/10.Microdiversity/")
data1 <- read.table('metapop-pi-average.txt',header = T,sep = "\t")
data1$group <-factor(data1$group,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
color1 <- brewer.pal(12,"Set3")

#箱线图
ggplot(data1,aes(group,value,fill=group))+
  geom_boxplot()+
  scale_fill_manual(values = color1)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+ylab("value of π")+
  xlab("Location")


#ldg

ggplot(data1,aes(x=Latitude,y=value))+
  geom_point(aes(colour=data1$group))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2),color="black")+
  scale_color_manual(values = color1)+ggprism::theme_prism()+
  ylab("value of π")+
  stat_cor()
  

#  stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = '~~~~')),label.y = "bottom", label.x = "left")


#读入richness 相关数据
data2<- read.table("read_richness",header = T,row.names = 1)
data2$ID <-row.names(data2)

data3 <- merge(data1,data2,by="ID")
ggplot(data3,aes(x=Richness11,y=value))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ poly(x),color="black")+
  ggprism::theme_prism()+
  ylab("value of pi")+
  xlab("value of richness")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )
