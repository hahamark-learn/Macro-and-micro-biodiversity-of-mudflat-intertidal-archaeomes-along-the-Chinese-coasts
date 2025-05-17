'''
The correlation between microdiversity and macrodiversity.
'''
#Load the corresponding package
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(ggpubr)
library(ggpmisc)
#Load the corresponding folder 
setwd("~Microdiversity/")
data1 <- read.table('metapop-pi-average-1.txt',header = T,sep = "\t")
data1$group <-factor(data1$group,levels=c("Sanya", "Beihai" ,"Zhuhai","Shantou","Xiamen","Wenzhou","Ningbo","Yancheng","Lianyungang",
                                          "Qingdao","Dongying","Dandong"))
color1 <- brewer.pal(12,"Set3")

#Plot a boxplot
p2 <- ggplot(data1,aes(group,value,fill=group))+
  geom_boxplot()+
  scale_fill_manual(values = color1)+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))+ylab("Nucleotide diversity")+
  xlab("Location")

ggsave("micro_2025.pdf",p2,width = 8,height = 6)

#Load the richness-related data
data2<- read.table("read_richness",header = T,row.names = 1)
data2$ID <-row.names(data2)
data3 <- merge(data1,data2,by="ID")
p1<-ggplot(data3,aes(x=Richness11,y=value))+
  geom_point()+
  geom_smooth(method = "lm",formula =y ~ poly(x),color="black")+
  ggprism::theme_prism()+
  ylab("Nucleotide diversity")+
  xlab("Community richness")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = T
  )
p1
ggsave("relationship_macro_micro_2025.5.9.pdf",p1,width = 8,height = 6)
