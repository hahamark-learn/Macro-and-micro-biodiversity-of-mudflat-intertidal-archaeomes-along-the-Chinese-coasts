library(reshape2)
library(dplyr)
library(tidyr)
library(circlize)
library(RColorBrewer)
mat = matrix(sample(18,18),3,6)
rownames(mat) <- paste0("R",1:3)
colnames(mat) <- paste0("C",1:6)


chordDiagram(mat)
# 长数据变短数据
setwd('~/Desktop/mag-hgt/')
data1 <- read.table('donor_recipient',header = T)
data2<- dcast(
  data = data1,
  donor~recipient
)
rownames(data2)<-data2$donor
data2<- data2[-1]
chordDiagram(data2)


#COG  数量
library(ggplot2)
library(scales)
library(cowplot)
data1 <- read.table('hgt_cog_number',header = T)
p6<-ggplot(data1,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                      name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")
p6
ggsave(p6,file="hgt-cog-1-25.2.23.pdf",width=12,height=8)



#hgtector 可视化分析
#Thaumarchaeota
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/hgtector")
data2 <- read.table('thaumarchaeota_cog_number_10',header = T)
p1<-ggplot(data2,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")
p1
ggsave(p1,file="hgt-thaumarchaeota-cog-25.2.23.pdf",width=12,height=8)






#Euryarchaeota
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/hgtector")
data3 <- read.table('euryarchaeota_cog_number_10',header = T)
p2<-ggplot(data3,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")
p2
ggsave(p,file="hgt-euryarchaeota-cog-25.2.23.pdf",width=12,height=8)


#Nitrososphaerota
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/hgtector")
data4 <- read.table('Nitrososphaerota_raw_number_10',header = T)
p3<-ggplot(data4,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")
p3
ggsave(p3,file="hgt-Nitrososphaerota-cog-25.2.23.pdf",width=12,height=8)





#metachip_archaea_bacteria
data5 <- read.table("/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/metachip_tidal_archaea_bacteria/hgt_count_clean",sep = " ",header = T)
p4<-ggplot(data5,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")
ggsave(p4,file="/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/metachip_tidal_archaea_bacteria/metachip-archaea-bacteria-cog-25.2.23.pdf",width=12,height=8)


#metachip_archaea_to_archaea
data6 <- read.table("/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/metachip_tidal_archaea/cog_clean",sep = " ",header = T)
p5<-ggplot(data6,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")
ggsave(p5,file="/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/metachip_tidal_archaea/metachip-archaea-archaea-cog-25.2.23.pdf",width=12,height=8)


#绘制分面图
data2$phylum <- "Thaumarchaeota"
data3$phylum <- "Euryarchaeota"
data4$phylum <- "Nitrososphaerota"
hgtector_all<- rbind(data2,data3,data4)

hgtector_all_1<- hgtector_all %>%
  group_by(phylum) %>%
  mutate(group = factor(group, levels = unique(group[order(-value)])))


hgtector_all_1$phylum<- factor(hgtector_all_1$phylum,levels = c("thaumarchaeota",
                                                                "euryarchaeota",
                                                                "Nitrososphaerota"))
p7<-ggplot(hgtector_all_1,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")+
  facet_grid(~phylum,scales = "free_x")

ggsave(p7,file="/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/hgtector_all-cog-25.2.23.pdf",width=12,height=8)


#绘制metachip的图片

data5$phylum <- "AB"
data6$phylum <- "AA"
metachip_all<- rbind(data5,data6)

metachip_all_1<- metachip_all %>%
  group_by(phylum) %>%
  mutate(group = factor(group, levels = unique(group[order(-value)])))

metachip_all_1$phylum<- factor(metachip_all_1$phylum,levels = c("AB",
                                                                "AA"))

p8<-ggplot(metachip_all_1,aes(x=reorder(group,value,decreasing =T),y=value))+
  geom_col(aes(fill=value))+
  scale_fill_gradient2(low=muted("#3E89C2"),mid="#FEC554",high = muted("#E47779"),
                       name=expression(italic("number")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x="",y="Number of COG functions")+
  facet_grid(~phylum,scales = "free_x")

ggsave(p8,file="/Users/hahamark/Desktop/article1/宏基因组分析/基于MAG-HGT分析/metachip_all-cog-25.2.23.pdf",width=12,height=8)




#修改做图代码，修改为散点图

hgtector_all$type <- "hgtector"
metachip_all$type <- "metachip"
hgt_all<- rbind(hgtector_all,metachip_all)

#设置渐变性颜色：
color1 <- brewer.pal(10,"Paired")
color2 <- colorRampPalette(color1)(19)
# 绘制散点图
p1 <- ggplot(hgt_all, aes(x = phylum, y = group, size = value, color = group)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.15)) +  # 添加抖动避免点重叠
  scale_size_continuous(range = c(3, 18))+
  labs(x = "Category", y = "Value", size = "Size", color = "Color") +  # 设置轴标签和图例标题
  facet_grid(~type,scales = "free_x")+
  scale_color_manual(values=color2)

ggsave("hgt_plot_25.3.11.pdf", plot = p1, width = 210 / 25.4, height = 297 / 25.4)

