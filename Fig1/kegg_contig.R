#Based on information obtained from contigs, the functional profile of overall tidal flat archaea at the continental scale was plotted
#Load the required R packages
library(tidyverse)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
#合并所有的kegg_pathway注释
setwd("/Users/hahamark/Desktop/article1/ghostkoala/")
files <- list.files(pattern = "kegg_clean_result_standardized_pathway")
combined_data <- data.frame()
for (i in files){
  data <- read.table(i,header = T,row.names = 1,check.names = F)
  if (nrow(combined_data)==0){
    combined_data <- data
  }else{
    combined_data <- cbind(combined_data,data)
  }
  
}

#按照数量大小排序：
order<-sort(rowSums(combined_data[,1:ncol(combined_data)]),index.return=TRUE,decreasing=T) 
combined_data1<-combined_data[order$ix,]
combined_data2 <- combined_data1[1:20,]
combined_data3 <- rownames_to_column(combined_data2,var = "kegg_id")
combined_data4 <- melt(combined_data3,id.vars= "kegg_id",value.name = "value",variable.name = "place")


combined_data5 <- combined_data4 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)


my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")

combined_data5$place <- factor(combined_data5$place,levels=c("SY11_FDME210455407-1r_1.fq.gz","SY12_FDME210455408-1r_1.fq.gz","SY13_FDME210455409-2a_1.fq.gz",
                                                              "SY15_FDME210455411-2a_1.fq.gz","SY1_FDME210455397-2a_1.fq.gz","SY3_FDME210455399-2a_1.fq.gz",
                                                              "SY7_FDME210455403-2a_1.fq.gz","SY9_FDME210455405-2a_1.fq.gz","BH10_FDSW220017130-1r_1.fq.gz",
                                                              "BH12_FDSW220017132-1r_1.fq.gz","BH14_FDSW220017134-1r_1.fq.gz","BH15_FDSW220017135-1r_1.fq.gz",
                                                              "BH2_FDSW220017122-1r_1.fq.gz","BH4_FDSW220017124-1r_1.fq.gz","BH6_FDSW220017126-1r_1.fq.gz",
                                                              "BH8_FDSW220017128-1r_1.fq.gz","ZH11_FDSW220008794-2r_1.fq.gz","ZH13_FDSW220008796-2r_1.fq.gz",
                                                              "ZH15_FDSW220008798-2r_1.fq.gz","ZH1_FDSW220008784-1r_1.fq.gz","ZH3_FDSW220008786-1r_1.fq.gz",
                                                              "ZH5_FDSW220008788-1r_1.fq.gz","ZH7_FDSW220008790-1r_1.fq.gz","ZH9_FDSW220008792-2r_1.fq.gz",
                                                              "ST10_L1_1.fq.gz","ST12_L1_1.fq.gz","ST13_L1_1.fq.gz","ST15_L1_1.fq.gz","ST1_FDME220007465-1r_1.fq.gz",
                                                              "ST3_L1_1.fq.gz","ST7_FDME220007471-1r_1.fq.gz","ST9_L1_1.fq.gz","XM10_FDME220007444-1r_1.fq.gz",
                                                              "XM12_FDME220007446-1r_1.fq.gz","XM14_FDME220007448-1r_1.fq.gz","XM15_FDME220007449-1r_1.fq.gz",
                                                              "XM2_FDME220007436-1r_1.fq.gz","XM4_FDME220007438-1r_1.fq.gz","XM6_FDME220007440-1r_1.fq.gz",
                                                              "XM8_FDME220007442-1r_1.fq.gz","WZ10_L1_1.fq.gz","WZ11_L1_1.fq.gz","WZ14_L1_1.fq.gz",
                                                              "WZ15_FDME210455471-1r_1.fq.gz","WZ1_L1_1.fq.gz","WZ3_L1_1.fq.gz","WZ5_L1_1.fq.gz",
                                                              "WZ7_FDME210455463-1r_1.fq.gz","NB11_FDME210455362-1b_1.fq.gz","NB12_FDME210455363-1r_1.fq.gz",
                                                              "NB15_FDME210455366-1r_1.fq.gz","NB1_FDME210455352-1r_1.fq.gz","NB2_FDME210455353-1r_1.fq.gz",
                                                              "NB4_FDME210455355-1r_1.fq.gz","NB6_FDME210455357-1r_1.fq.gz","NB7_FDME210455358-1r_1.fq.gz",
                                                              "YC11_FDME220007415-1r_1.fq.gz","YC13_FDME220007417-1r_1.fq.gz","YC15_FDME220007419-1r_1.fq.gz",
                                                              "YC1_FDME220007405-1r_1.fq.gz","YC3_FDME220007407-1r_1.fq.gz","YC5_FDME220007409-1r_1.fq.gz",
                                                              "YC7_FDME220007411-1a_1.fq.gz","YC9_FDME220007413-1r_1.fq.gz","LYG11_FDME210455392-2r_1.fq.gz",
                                                              "LYG13_FDME210455394-1a_1.fq.gz","LYG15_FDME210455396-1a_1.fq.gz","LYG1_FDME210455382-1r_1.fq.gz",
                                                              "LYG3_FDME210455384-1r_1.fq.gz","LYG5_FDME210455386-1r_1.fq.gz","LYG6_FDME210455387-1r_1.fq.gz",
                                                              "LYG7_FDME210455388-1r_1.fq.gz","CXD_1412_FDSW210423229-1r_1.fq.gz","CXD_1414_FDSW210423231-1r_1.fq.gz",
                                                              "CXD_1604_FDSW210423236-1r_1.fq.gz","CXD_1609_FDSW210423241-1r_1.fq.gz","CXD_1610_FDSW210423242-1r_1.fq.gz",
                                                              "CXD_1611_FDSW210423243-1r_1.fq.gz","CXD_1612_FDSW210423244-1r_1.fq.gz","CXD_1613_FDSW210423245-1r_1.fq.gz",
                                                              "DY10_FDME210438608-1r_1.fq.gz","DY11_FDME210438609-1r_1.fq.gz","DY13_FDME210438611-1b_1.fq.gz","DY15_FDME210438613-1r_1.fq.gz",
                                                              "DY2_FDME210438600-1b_1.fq.gz","DY3_FDME210438601-1r_1.fq.gz","DY6_FDME210438604-1r_1.fq.gz","DY7_FDME210438605-1r_1.fq.gz",
                                                              "DD11_FDSW220008398-1r_1.fq.gz","DD13_FDSW220008400-2r_1.fq.gz","DD15_FDSW220008402-2r_1.fq.gz","DD1_FDSW220008388-1r_1.fq.gz",
                                                              "DD3_FDSW220008390-2r_1.fq.gz","DD5_FDSW220008392-2r_1.fq.gz","DD7_FDSW220008394-1r_1.fq.gz","DD9_FDSW220008396-1r_1.fq.gz"
  
))

ggplot(combined_data5,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  scale_fill_manual(values = my_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave("kegg_contig_count.pdf", width = 8.27, height = 11.69)





#ko统计，绘制前50 ko

#合并所有的kegg注释
setwd("/Users/hahamark/Desktop/article1/ghostkoala/")
files <- list.files(pattern = "kegg_clean_result_standardized$")
combined_data <- data.frame()
for (i in files){
  data <- read.table(i,header = T,row.names = 1,check.names = F)
  if (nrow(combined_data)==0){
    combined_data <- data
  }else{
    combined_data <- cbind(combined_data,data)
  }
  
}

#写入所有的kegg,到文件中，方便后续生物地理学模式的计算
write.table(combined_data,"all_kegg_contig.txt",quote = F,sep = "\t")
#按照数量大小排序：
order<-sort(rowSums(combined_data[,1:ncol(combined_data)]),index.return=TRUE,decreasing=T) 
combined_data1<-combined_data[order$ix,]
combined_data2 <- combined_data1[1:30,]
combined_data3 <- rownames_to_column(combined_data2,var = "kegg_id")
combined_data4 <- melt(combined_data3,id.vars= "kegg_id",value.name = "value",variable.name = "place")

combined_data5 <- combined_data4 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)

my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#6b6ecf", "#9c9ede", "#bd9e39", "#ad494a"
)

combined_data5$place <- factor(combined_data5$place,levels=c("SY11_FDME210455407-1r_1.fq.gz","SY12_FDME210455408-1r_1.fq.gz","SY13_FDME210455409-2a_1.fq.gz",
                                                             "SY15_FDME210455411-2a_1.fq.gz","SY1_FDME210455397-2a_1.fq.gz","SY3_FDME210455399-2a_1.fq.gz",
                                                             "SY7_FDME210455403-2a_1.fq.gz","SY9_FDME210455405-2a_1.fq.gz","BH10_FDSW220017130-1r_1.fq.gz",
                                                             "BH12_FDSW220017132-1r_1.fq.gz","BH14_FDSW220017134-1r_1.fq.gz","BH15_FDSW220017135-1r_1.fq.gz",
                                                             "BH2_FDSW220017122-1r_1.fq.gz","BH4_FDSW220017124-1r_1.fq.gz","BH6_FDSW220017126-1r_1.fq.gz",
                                                             "BH8_FDSW220017128-1r_1.fq.gz","ZH11_FDSW220008794-2r_1.fq.gz","ZH13_FDSW220008796-2r_1.fq.gz",
                                                             "ZH15_FDSW220008798-2r_1.fq.gz","ZH1_FDSW220008784-1r_1.fq.gz","ZH3_FDSW220008786-1r_1.fq.gz",
                                                             "ZH5_FDSW220008788-1r_1.fq.gz","ZH7_FDSW220008790-1r_1.fq.gz","ZH9_FDSW220008792-2r_1.fq.gz",
                                                             "ST10_L1_1.fq.gz","ST12_L1_1.fq.gz","ST13_L1_1.fq.gz","ST15_L1_1.fq.gz","ST1_FDME220007465-1r_1.fq.gz",
                                                             "ST3_L1_1.fq.gz","ST7_FDME220007471-1r_1.fq.gz","ST9_L1_1.fq.gz","XM10_FDME220007444-1r_1.fq.gz",
                                                             "XM12_FDME220007446-1r_1.fq.gz","XM14_FDME220007448-1r_1.fq.gz","XM15_FDME220007449-1r_1.fq.gz",
                                                             "XM2_FDME220007436-1r_1.fq.gz","XM4_FDME220007438-1r_1.fq.gz","XM6_FDME220007440-1r_1.fq.gz",
                                                             "XM8_FDME220007442-1r_1.fq.gz","WZ10_L1_1.fq.gz","WZ11_L1_1.fq.gz","WZ14_L1_1.fq.gz",
                                                             "WZ15_FDME210455471-1r_1.fq.gz","WZ1_L1_1.fq.gz","WZ3_L1_1.fq.gz","WZ5_L1_1.fq.gz",
                                                             "WZ7_FDME210455463-1r_1.fq.gz","NB11_FDME210455362-1b_1.fq.gz","NB12_FDME210455363-1r_1.fq.gz",
                                                             "NB15_FDME210455366-1r_1.fq.gz","NB1_FDME210455352-1r_1.fq.gz","NB2_FDME210455353-1r_1.fq.gz",
                                                             "NB4_FDME210455355-1r_1.fq.gz","NB6_FDME210455357-1r_1.fq.gz","NB7_FDME210455358-1r_1.fq.gz",
                                                             "YC11_FDME220007415-1r_1.fq.gz","YC13_FDME220007417-1r_1.fq.gz","YC15_FDME220007419-1r_1.fq.gz",
                                                             "YC1_FDME220007405-1r_1.fq.gz","YC3_FDME220007407-1r_1.fq.gz","YC5_FDME220007409-1r_1.fq.gz",
                                                             "YC7_FDME220007411-1a_1.fq.gz","YC9_FDME220007413-1r_1.fq.gz","LYG11_FDME210455392-2r_1.fq.gz",
                                                             "LYG13_FDME210455394-1a_1.fq.gz","LYG15_FDME210455396-1a_1.fq.gz","LYG1_FDME210455382-1r_1.fq.gz",
                                                             "LYG3_FDME210455384-1r_1.fq.gz","LYG5_FDME210455386-1r_1.fq.gz","LYG6_FDME210455387-1r_1.fq.gz",
                                                             "LYG7_FDME210455388-1r_1.fq.gz","CXD_1412_FDSW210423229-1r_1.fq.gz","CXD_1414_FDSW210423231-1r_1.fq.gz",
                                                             "CXD_1604_FDSW210423236-1r_1.fq.gz","CXD_1609_FDSW210423241-1r_1.fq.gz","CXD_1610_FDSW210423242-1r_1.fq.gz",
                                                             "CXD_1611_FDSW210423243-1r_1.fq.gz","CXD_1612_FDSW210423244-1r_1.fq.gz","CXD_1613_FDSW210423245-1r_1.fq.gz",
                                                             "DY10_FDME210438608-1r_1.fq.gz","DY11_FDME210438609-1r_1.fq.gz","DY13_FDME210438611-1b_1.fq.gz","DY15_FDME210438613-1r_1.fq.gz",
                                                             "DY2_FDME210438600-1b_1.fq.gz","DY3_FDME210438601-1r_1.fq.gz","DY6_FDME210438604-1r_1.fq.gz","DY7_FDME210438605-1r_1.fq.gz",
                                                             "DD11_FDSW220008398-1r_1.fq.gz","DD13_FDSW220008400-2r_1.fq.gz","DD15_FDSW220008402-2r_1.fq.gz","DD1_FDSW220008388-1r_1.fq.gz",
                                                             "DD3_FDSW220008390-2r_1.fq.gz","DD5_FDSW220008392-2r_1.fq.gz","DD7_FDSW220008394-1r_1.fq.gz","DD9_FDSW220008396-1r_1.fq.gz"
                                                             
))

ggplot(combined_data5,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  scale_fill_manual(values = my_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave("kegg_contig.pdf", width = 8.27, height = 11.69)









#使用统计完的数据进行，并且平均后的数据

library(tidyverse)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
#合并所有的kegg_pathway注释
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/cd_hit_result/kegg/")
data1 <- read.table("all_pathway_kegg",header = T,row.names = 1,check.names = F)
order<-sort(rowSums(data1[,1:ncol(data1)]),index.return=TRUE,decreasing=T) 
combined_data1<-data1[order$ix,]
combined_data2 <- combined_data1[1:20,]
data2 <- rownames_to_column(combined_data2,var = "kegg_id")
data3 <- melt(data2,id.vars= "kegg_id",value.name = "value",variable.name = "place")


data4 <- data3 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)


my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")

data4$place <- factor(data4$place,levels=c("SY11_FDME210455407-1r_1.fq.gz","SY12_FDME210455408-1r_1.fq.gz","SY13_FDME210455409-2a_1.fq.gz",
                                                             "SY15_FDME210455411-2a_1.fq.gz","SY1_FDME210455397-2a_1.fq.gz","SY3_FDME210455399-2a_1.fq.gz",
                                                             "SY7_FDME210455403-2a_1.fq.gz","SY9_FDME210455405-2a_1.fq.gz","BH10_FDSW220017130-1r_1.fq.gz",
                                                             "BH12_FDSW220017132-1r_1.fq.gz","BH14_FDSW220017134-1r_1.fq.gz","BH15_FDSW220017135-1r_1.fq.gz",
                                                             "BH2_FDSW220017122-1r_1.fq.gz","BH4_FDSW220017124-1r_1.fq.gz","BH6_FDSW220017126-1r_1.fq.gz",
                                                             "BH8_FDSW220017128-1r_1.fq.gz","ZH11_FDSW220008794-2r_1.fq.gz","ZH13_FDSW220008796-2r_1.fq.gz",
                                                             "ZH15_FDSW220008798-2r_1.fq.gz","ZH1_FDSW220008784-1r_1.fq.gz","ZH3_FDSW220008786-1r_1.fq.gz",
                                                             "ZH5_FDSW220008788-1r_1.fq.gz","ZH7_FDSW220008790-1r_1.fq.gz","ZH9_FDSW220008792-2r_1.fq.gz",
                                                             "ST10_L1_1.fq.gz","ST12_L1_1.fq.gz","ST13_L1_1.fq.gz","ST15_L1_1.fq.gz","ST1_FDME220007465-1r_1.fq.gz",
                                                             "ST3_L1_1.fq.gz","ST7_FDME220007471-1r_1.fq.gz","ST9_L1_1.fq.gz","XM10_FDME220007444-1r_1.fq.gz",
                                                             "XM12_FDME220007446-1r_1.fq.gz","XM14_FDME220007448-1r_1.fq.gz","XM15_FDME220007449-1r_1.fq.gz",
                                                             "XM2_FDME220007436-1r_1.fq.gz","XM4_FDME220007438-1r_1.fq.gz","XM6_FDME220007440-1r_1.fq.gz",
                                                             "XM8_FDME220007442-1r_1.fq.gz","WZ10_L1_1.fq.gz","WZ11_L1_1.fq.gz","WZ14_L1_1.fq.gz",
                                                             "WZ15_FDME210455471-1r_1.fq.gz","WZ1_L1_1.fq.gz","WZ3_L1_1.fq.gz","WZ5_L1_1.fq.gz",
                                                             "WZ7_FDME210455463-1r_1.fq.gz","NB11_FDME210455362-1b_1.fq.gz","NB12_FDME210455363-1r_1.fq.gz",
                                                             "NB15_FDME210455366-1r_1.fq.gz","NB1_FDME210455352-1r_1.fq.gz","NB2_FDME210455353-1r_1.fq.gz",
                                                             "NB4_FDME210455355-1r_1.fq.gz","NB6_FDME210455357-1r_1.fq.gz","NB7_FDME210455358-1r_1.fq.gz",
                                                             "YC11_FDME220007415-1r_1.fq.gz","YC13_FDME220007417-1r_1.fq.gz","YC15_FDME220007419-1r_1.fq.gz",
                                                             "YC1_FDME220007405-1r_1.fq.gz","YC3_FDME220007407-1r_1.fq.gz","YC5_FDME220007409-1r_1.fq.gz",
                                                             "YC7_FDME220007411-1a_1.fq.gz","YC9_FDME220007413-1r_1.fq.gz","LYG11_FDME210455392-2r_1.fq.gz",
                                                             "LYG13_FDME210455394-1a_1.fq.gz","LYG15_FDME210455396-1a_1.fq.gz","LYG1_FDME210455382-1r_1.fq.gz",
                                                             "LYG3_FDME210455384-1r_1.fq.gz","LYG5_FDME210455386-1r_1.fq.gz","LYG6_FDME210455387-1r_1.fq.gz",
                                                             "LYG7_FDME210455388-1r_1.fq.gz","CXD_1412_FDSW210423229-1r_1.fq.gz","CXD_1414_FDSW210423231-1r_1.fq.gz",
                                                             "CXD_1604_FDSW210423236-1r_1.fq.gz","CXD_1609_FDSW210423241-1r_1.fq.gz","CXD_1610_FDSW210423242-1r_1.fq.gz",
                                                             "CXD_1611_FDSW210423243-1r_1.fq.gz","CXD_1612_FDSW210423244-1r_1.fq.gz","CXD_1613_FDSW210423245-1r_1.fq.gz",
                                                             "DY10_FDME210438608-1r_1.fq.gz","DY11_FDME210438609-1r_1.fq.gz","DY13_FDME210438611-1b_1.fq.gz","DY15_FDME210438613-1r_1.fq.gz",
                                                             "DY2_FDME210438600-1b_1.fq.gz","DY3_FDME210438601-1r_1.fq.gz","DY6_FDME210438604-1r_1.fq.gz","DY7_FDME210438605-1r_1.fq.gz",
                                                             "DD11_FDSW220008398-1r_1.fq.gz","DD13_FDSW220008400-2r_1.fq.gz","DD15_FDSW220008402-2r_1.fq.gz","DD1_FDSW220008388-1r_1.fq.gz",
                                                             "DD3_FDSW220008390-2r_1.fq.gz","DD5_FDSW220008392-2r_1.fq.gz","DD7_FDSW220008394-1r_1.fq.gz","DD9_FDSW220008396-1r_1.fq.gz"
                                                             
))

ggplot(data4,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  scale_fill_manual(values = my_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave("kegg_contig_pathway_count_25.2.27.pdf", width = 8.27, height = 11.69)



library(tidyverse)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
#合并所有的kegg_pathway注释
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/cd_hit_result/kegg/")
data1 <- read.table("all_contig_kegg.txt",header = T,row.names = 1,check.names = F)
order<-sort(rowSums(data1[,1:ncol(data1)]),index.return=TRUE,decreasing=T) 
combined_data1<-data1[order$ix,]
combined_data2 <- combined_data1[1:20,]
data2 <- rownames_to_column(combined_data2,var = "kegg_id")
data3 <- melt(data2,id.vars= "kegg_id",value.name = "value",variable.name = "place")


data4 <- data3 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)

my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#6b6ecf", "#9c9ede", "#bd9e39", "#ad494a"
)


data4$place <- factor(data4$place,levels=c("SY11_FDME210455407-1r_1.fq.gz","SY12_FDME210455408-1r_1.fq.gz","SY13_FDME210455409-2a_1.fq.gz",
                                           "SY15_FDME210455411-2a_1.fq.gz","SY1_FDME210455397-2a_1.fq.gz","SY3_FDME210455399-2a_1.fq.gz",
                                           "SY7_FDME210455403-2a_1.fq.gz","SY9_FDME210455405-2a_1.fq.gz","BH10_FDSW220017130-1r_1.fq.gz",
                                           "BH12_FDSW220017132-1r_1.fq.gz","BH14_FDSW220017134-1r_1.fq.gz","BH15_FDSW220017135-1r_1.fq.gz",
                                           "BH2_FDSW220017122-1r_1.fq.gz","BH4_FDSW220017124-1r_1.fq.gz","BH6_FDSW220017126-1r_1.fq.gz",
                                           "BH8_FDSW220017128-1r_1.fq.gz","ZH11_FDSW220008794-2r_1.fq.gz","ZH13_FDSW220008796-2r_1.fq.gz",
                                           "ZH15_FDSW220008798-2r_1.fq.gz","ZH1_FDSW220008784-1r_1.fq.gz","ZH3_FDSW220008786-1r_1.fq.gz",
                                           "ZH5_FDSW220008788-1r_1.fq.gz","ZH7_FDSW220008790-1r_1.fq.gz","ZH9_FDSW220008792-2r_1.fq.gz",
                                           "ST10_L1_1.fq.gz","ST12_L1_1.fq.gz","ST13_L1_1.fq.gz","ST15_L1_1.fq.gz","ST1_FDME220007465-1r_1.fq.gz",
                                           "ST3_L1_1.fq.gz","ST7_FDME220007471-1r_1.fq.gz","ST9_L1_1.fq.gz","XM10_FDME220007444-1r_1.fq.gz",
                                           "XM12_FDME220007446-1r_1.fq.gz","XM14_FDME220007448-1r_1.fq.gz","XM15_FDME220007449-1r_1.fq.gz",
                                           "XM2_FDME220007436-1r_1.fq.gz","XM4_FDME220007438-1r_1.fq.gz","XM6_FDME220007440-1r_1.fq.gz",
                                           "XM8_FDME220007442-1r_1.fq.gz","WZ10_L1_1.fq.gz","WZ11_L1_1.fq.gz","WZ14_L1_1.fq.gz",
                                           "WZ15_FDME210455471-1r_1.fq.gz","WZ1_L1_1.fq.gz","WZ3_L1_1.fq.gz","WZ5_L1_1.fq.gz",
                                           "WZ7_FDME210455463-1r_1.fq.gz","NB11_FDME210455362-1b_1.fq.gz","NB12_FDME210455363-1r_1.fq.gz",
                                           "NB15_FDME210455366-1r_1.fq.gz","NB1_FDME210455352-1r_1.fq.gz","NB2_FDME210455353-1r_1.fq.gz",
                                           "NB4_FDME210455355-1r_1.fq.gz","NB6_FDME210455357-1r_1.fq.gz","NB7_FDME210455358-1r_1.fq.gz",
                                           "YC11_FDME220007415-1r_1.fq.gz","YC13_FDME220007417-1r_1.fq.gz","YC15_FDME220007419-1r_1.fq.gz",
                                           "YC1_FDME220007405-1r_1.fq.gz","YC3_FDME220007407-1r_1.fq.gz","YC5_FDME220007409-1r_1.fq.gz",
                                           "YC7_FDME220007411-1a_1.fq.gz","YC9_FDME220007413-1r_1.fq.gz","LYG11_FDME210455392-2r_1.fq.gz",
                                           "LYG13_FDME210455394-1a_1.fq.gz","LYG15_FDME210455396-1a_1.fq.gz","LYG1_FDME210455382-1r_1.fq.gz",
                                           "LYG3_FDME210455384-1r_1.fq.gz","LYG5_FDME210455386-1r_1.fq.gz","LYG6_FDME210455387-1r_1.fq.gz",
                                           "LYG7_FDME210455388-1r_1.fq.gz","CXD_1412_FDSW210423229-1r_1.fq.gz","CXD_1414_FDSW210423231-1r_1.fq.gz",
                                           "CXD_1604_FDSW210423236-1r_1.fq.gz","CXD_1609_FDSW210423241-1r_1.fq.gz","CXD_1610_FDSW210423242-1r_1.fq.gz",
                                           "CXD_1611_FDSW210423243-1r_1.fq.gz","CXD_1612_FDSW210423244-1r_1.fq.gz","CXD_1613_FDSW210423245-1r_1.fq.gz",
                                           "DY10_FDME210438608-1r_1.fq.gz","DY11_FDME210438609-1r_1.fq.gz","DY13_FDME210438611-1b_1.fq.gz","DY15_FDME210438613-1r_1.fq.gz",
                                           "DY2_FDME210438600-1b_1.fq.gz","DY3_FDME210438601-1r_1.fq.gz","DY6_FDME210438604-1r_1.fq.gz","DY7_FDME210438605-1r_1.fq.gz",
                                           "DD11_FDSW220008398-1r_1.fq.gz","DD13_FDSW220008400-2r_1.fq.gz","DD15_FDSW220008402-2r_1.fq.gz","DD1_FDSW220008388-1r_1.fq.gz",
                                           "DD3_FDSW220008390-2r_1.fq.gz","DD5_FDSW220008392-2r_1.fq.gz","DD7_FDSW220008394-1r_1.fq.gz","DD9_FDSW220008396-1r_1.fq.gz"
                                           
))

ggplot(data4,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  scale_fill_manual(values = my_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave("kegg_contig_ko_25.2.27.pdf", width = 8.27, height = 11.69)


#绘制一个整体的数值

#合并所有的kegg_pathway注释
setwd("/Users/hahamark/Desktop/article1/宏基因组分析/contig_taxnonmy/cd_hit_result/kegg/")
data1 <- read.table("all_contig_kegg.txt",header = T,row.names = 1,check.names = F)
order<-sort(rowSums(data1[,1:ncol(data1)]),index.return=TRUE,decreasing=T) 
combined_data1<-data1[order$ix,]
data2 <- rownames_to_column(combined_data1,var = "kegg_id")
data3 <- melt(data2,id.vars= "kegg_id",value.name = "value",variable.name = "place")

colors <- brewer.pal(12,"Set3")
colors1 <- colorRampPalette(colors)(7510)
data4 <- data3 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)

data4$place <- factor(data4$place,levels=c("SY11_FDME210455407-1r_1.fq.gz","SY12_FDME210455408-1r_1.fq.gz","SY13_FDME210455409-2a_1.fq.gz",
                                           "SY15_FDME210455411-2a_1.fq.gz","SY1_FDME210455397-2a_1.fq.gz","SY3_FDME210455399-2a_1.fq.gz",
                                           "SY7_FDME210455403-2a_1.fq.gz","SY9_FDME210455405-2a_1.fq.gz","BH10_FDSW220017130-1r_1.fq.gz",
                                           "BH12_FDSW220017132-1r_1.fq.gz","BH14_FDSW220017134-1r_1.fq.gz","BH15_FDSW220017135-1r_1.fq.gz",
                                           "BH2_FDSW220017122-1r_1.fq.gz","BH4_FDSW220017124-1r_1.fq.gz","BH6_FDSW220017126-1r_1.fq.gz",
                                           "BH8_FDSW220017128-1r_1.fq.gz","ZH11_FDSW220008794-2r_1.fq.gz","ZH13_FDSW220008796-2r_1.fq.gz",
                                           "ZH15_FDSW220008798-2r_1.fq.gz","ZH1_FDSW220008784-1r_1.fq.gz","ZH3_FDSW220008786-1r_1.fq.gz",
                                           "ZH5_FDSW220008788-1r_1.fq.gz","ZH7_FDSW220008790-1r_1.fq.gz","ZH9_FDSW220008792-2r_1.fq.gz",
                                           "ST10_L1_1.fq.gz","ST12_L1_1.fq.gz","ST13_L1_1.fq.gz","ST15_L1_1.fq.gz","ST1_FDME220007465-1r_1.fq.gz",
                                           "ST3_L1_1.fq.gz","ST7_FDME220007471-1r_1.fq.gz","ST9_L1_1.fq.gz","XM10_FDME220007444-1r_1.fq.gz",
                                           "XM12_FDME220007446-1r_1.fq.gz","XM14_FDME220007448-1r_1.fq.gz","XM15_FDME220007449-1r_1.fq.gz",
                                           "XM2_FDME220007436-1r_1.fq.gz","XM4_FDME220007438-1r_1.fq.gz","XM6_FDME220007440-1r_1.fq.gz",
                                           "XM8_FDME220007442-1r_1.fq.gz","WZ10_L1_1.fq.gz","WZ11_L1_1.fq.gz","WZ14_L1_1.fq.gz",
                                           "WZ15_FDME210455471-1r_1.fq.gz","WZ1_L1_1.fq.gz","WZ3_L1_1.fq.gz","WZ5_L1_1.fq.gz",
                                           "WZ7_FDME210455463-1r_1.fq.gz","NB11_FDME210455362-1b_1.fq.gz","NB12_FDME210455363-1r_1.fq.gz",
                                           "NB15_FDME210455366-1r_1.fq.gz","NB1_FDME210455352-1r_1.fq.gz","NB2_FDME210455353-1r_1.fq.gz",
                                           "NB4_FDME210455355-1r_1.fq.gz","NB6_FDME210455357-1r_1.fq.gz","NB7_FDME210455358-1r_1.fq.gz",
                                           "YC11_FDME220007415-1r_1.fq.gz","YC13_FDME220007417-1r_1.fq.gz","YC15_FDME220007419-1r_1.fq.gz",
                                           "YC1_FDME220007405-1r_1.fq.gz","YC3_FDME220007407-1r_1.fq.gz","YC5_FDME220007409-1r_1.fq.gz",
                                           "YC7_FDME220007411-1a_1.fq.gz","YC9_FDME220007413-1r_1.fq.gz","LYG11_FDME210455392-2r_1.fq.gz",
                                           "LYG13_FDME210455394-1a_1.fq.gz","LYG15_FDME210455396-1a_1.fq.gz","LYG1_FDME210455382-1r_1.fq.gz",
                                           "LYG3_FDME210455384-1r_1.fq.gz","LYG5_FDME210455386-1r_1.fq.gz","LYG6_FDME210455387-1r_1.fq.gz",
                                           "LYG7_FDME210455388-1r_1.fq.gz","CXD_1412_FDSW210423229-1r_1.fq.gz","CXD_1414_FDSW210423231-1r_1.fq.gz",
                                           "CXD_1604_FDSW210423236-1r_1.fq.gz","CXD_1609_FDSW210423241-1r_1.fq.gz","CXD_1610_FDSW210423242-1r_1.fq.gz",
                                           "CXD_1611_FDSW210423243-1r_1.fq.gz","CXD_1612_FDSW210423244-1r_1.fq.gz","CXD_1613_FDSW210423245-1r_1.fq.gz",
                                           "DY10_FDME210438608-1r_1.fq.gz","DY11_FDME210438609-1r_1.fq.gz","DY13_FDME210438611-1b_1.fq.gz","DY15_FDME210438613-1r_1.fq.gz",
                                           "DY2_FDME210438600-1b_1.fq.gz","DY3_FDME210438601-1r_1.fq.gz","DY6_FDME210438604-1r_1.fq.gz","DY7_FDME210438605-1r_1.fq.gz",
                                           "DD11_FDSW220008398-1r_1.fq.gz","DD13_FDSW220008400-2r_1.fq.gz","DD15_FDSW220008402-2r_1.fq.gz","DD1_FDSW220008388-1r_1.fq.gz",
                                           "DD3_FDSW220008390-2r_1.fq.gz","DD5_FDSW220008392-2r_1.fq.gz","DD7_FDSW220008394-1r_1.fq.gz","DD9_FDSW220008396-1r_1.fq.gz"
                                           
))

p2 <- ggplot(data4,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),legend.position = "none")+
  scale_fill_manual(values = colors1)

ggsave("kegg_contig_ko_25.4.11.pdf",p2,width = 8.27, height = 11.69)







