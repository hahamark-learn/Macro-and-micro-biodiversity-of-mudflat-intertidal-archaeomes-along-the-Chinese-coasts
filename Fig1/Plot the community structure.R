'''
Plot the community structure of the entire coastal wetland region of China based on amplicon data
'''
#Load the corresponding R packages
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
#Load the normalized data (exclude the row where the sum is 0)
asv_clean <- read.table("ASVs_counts_clean_Flattening-rdp-new.txt",row.names = 1,header = T)
asv <- asv_clean
dat <- merge(x=asv,y=taxa,by='row.names')
dat=dplyr::rename(dat,OTUID=Row.names)
aa<-aggregate(dat[,2:145],by=list(dat$Phylum),FUN=sum)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
#Delete the row where the value is 0.
aa <- aa[-5,]
head(aa, n = 3)
#Sort the data based on the row sum results.
order<-sort(rowSums(aa[,1:ncol(aa)]),index.return=TRUE,decreasing=T)   
cc<-aa[order$ix,]
head(cc, n = 3)
bb<-merge(t(cc),dplyr::select(metadata,id,Location),
          by.x = "row.names",by.y ="id")
kk<-tidyr::gather(bb,Phylum,Abundance,-c(Location,Row.names))
kk$Phylum<-ifelse(kk$Phylum=='Unassigned','Others',kk$Phylum)
#Plot the data for 144 locations
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

