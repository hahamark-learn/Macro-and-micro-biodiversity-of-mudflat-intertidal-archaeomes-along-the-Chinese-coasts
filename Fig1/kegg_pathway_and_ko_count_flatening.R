‘’‘
Plotting information for fig1C and fig1D patterns.
’‘’
#kegg_contig pathway plotting
#Load the required package
library(tidyverse)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
#Read in the relevant pathway data
combined_data <- read.table("all_pathway_kegg_flatening",header = T,row.names = 1,check.names = F)
#Process the relevant pathway data
order<-sort(rowSums(combined_data[,1:ncol(combined_data)]),index.return=TRUE,decreasing=T) 
combined_data1<-combined_data[order$ix,]
combined_data2 <- combined_data1[1:20,]
combined_data3 <- rownames_to_column(combined_data2,var = "kegg_id")
combined_data4 <- melt(combined_data3,id.vars= "kegg_id",value.name = "value",variable.name = "place")
combined_data5 <- combined_data4 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)
#Select 20 color codes
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
#Sort the data based on location information
combined_data5$place <- factor(combined_data5$place,levels=c("SY11_FDME210455407-1r","SY12_FDME210455408-1r","SY13_FDME210455409-2a",
                                                             "SY15_FDME210455411-2a","SY1_FDME210455397-2a","SY3_FDME210455399-2a",
                                                             "SY7_FDME210455403-2a","SY9_FDME210455405-2a","BH10_FDSW220017130-1r",
                                                             "BH12_FDSW220017132-1r","BH14_FDSW220017134-1r","BH15_FDSW220017135-1r",
                                                             "BH2_FDSW220017122-1r","BH4_FDSW220017124-1r","BH6_FDSW220017126-1r",
                                                             "BH8_FDSW220017128-1r","ZH11_FDSW220008794-2r","ZH13_FDSW220008796-2r",
                                                             "ZH15_FDSW220008798-2r","ZH1_FDSW220008784-1r","ZH3_FDSW220008786-1r",
                                                             "ZH5_FDSW220008788-1r","ZH7_FDSW220008790-1r","ZH9_FDSW220008792-2r",
                                                             "ST10_L1","ST12_L1","ST13_L1","ST15_L1","ST1_FDME220007465-1r",
                                                             "ST3_L1","ST7_FDME220007471-1r","ST9_L1","XM10_FDME220007444-1r",
                                                             "XM12_FDME220007446-1r","XM14_FDME220007448-1r","XM15_FDME220007449-1r",
                                                             "XM2_FDME220007436-1r","XM4_FDME220007438-1r","XM6_FDME220007440-1r",
                                                             "XM8_FDME220007442-1r","WZ10_L1","WZ11_L1","WZ14_L1",
                                                             "WZ15_FDME210455471-1r","WZ1_L1","WZ3_L1","WZ5_L1",
                                                             "WZ7_FDME210455463-1r","NB11_FDME210455362-1b","NB12_FDME210455363-1r",
                                                             "NB15_FDME210455366-1r","NB1_FDME210455352-1r","NB2_FDME210455353-1r",
                                                             "NB4_FDME210455355-1r","NB6_FDME210455357-1r","NB7_FDME210455358-1r",
                                                             "YC11_FDME220007415-1r","YC13_FDME220007417-1r","YC15_FDME220007419-1r",
                                                             "YC1_FDME220007405-1r","YC3_FDME220007407-1r","YC5_FDME220007409-1r",
                                                             "YC7_FDME220007411-1a","YC9_FDME220007413-1r","LYG11_FDME210455392-2r",
                                                             "LYG13_FDME210455394-1a","LYG15_FDME210455396-1a","LYG1_FDME210455382-1r",
                                                             "LYG3_FDME210455384-1r","LYG5_FDME210455386-1r","LYG6_FDME210455387-1r",
                                                             "LYG7_FDME210455388-1r","CXD_1412_FDSW210423229-1r","CXD_1414_FDSW210423231-1r",
                                                             "CXD_1604_FDSW210423236-1r","CXD_1609_FDSW210423241-1r","CXD_1610_FDSW210423242-1r",
                                                             "CXD_1611_FDSW210423243-1r","CXD_1612_FDSW210423244-1r","CXD_1613_FDSW210423245-1r",
                                                             "DY10_FDME210438608-1r","DY11_FDME210438609-1r","DY13_FDME210438611-1b","DY15_FDME210438613-1r",
                                                             "DY2_FDME210438600-1b","DY3_FDME210438601-1r","DY6_FDME210438604-1r","DY7_FDME210438605-1r",
                                                             "DD11_FDSW220008398-1r","DD13_FDSW220008400-2r","DD15_FDSW220008402-2r","DD1_FDSW220008388-1r",
                                                             "DD3_FDSW220008390-2r","DD5_FDSW220008392-2r","DD7_FDSW220008394-1r","DD9_FDSW220008396-1r"
))
#Draw the pattern
p1 <-ggplot(combined_data5,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  scale_fill_manual(values = my_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#Display the pattern
p1
#Save the PDF file
ggsave("kegg_contig_count_flatening.pdf", p1 ,width = 8.27, height = 11.69)


#Plot the KEGG functions of all pathways in the tidal flat area (8021 KO)
#Process the relevant ko data
data1 <- read.table("all_contig_kegg_raw_flattening.txt",header = T,row.names = 1,check.names = F)
order<-sort(rowSums(data1[,1:ncol(data1)]),index.return=TRUE,decreasing=T) 
combined_data1<-data1[order$ix,]
data2 <- rownames_to_column(combined_data1,var = "kegg_id")
data3 <- melt(data2,id.vars= "kegg_id",value.name = "value",variable.name = "place")
data4 <- data3 %>% 
  group_by(place) %>%
  mutate(precent = value/sum(value)*100)
#Select 8012 color codes
colors <- brewer.pal(12,"Set3")
colors1 <- colorRampPalette(colors)(8012)
#Sort the data based on location information
data4$place <- factor(data4$place,levels=c("SY11_FDME210455407-1r","SY12_FDME210455408-1r","SY13_FDME210455409-2a",
                                           "SY15_FDME210455411-2a","SY1_FDME210455397-2a","SY3_FDME210455399-2a",
                                           "SY7_FDME210455403-2a","SY9_FDME210455405-2a","BH10_FDSW220017130-1r",
                                           "BH12_FDSW220017132-1r","BH14_FDSW220017134-1r","BH15_FDSW220017135-1r",
                                           "BH2_FDSW220017122-1r","BH4_FDSW220017124-1r","BH6_FDSW220017126-1r",
                                           "BH8_FDSW220017128-1r","ZH11_FDSW220008794-2r","ZH13_FDSW220008796-2r",
                                           "ZH15_FDSW220008798-2r","ZH1_FDSW220008784-1r","ZH3_FDSW220008786-1r",
                                           "ZH5_FDSW220008788-1r","ZH7_FDSW220008790-1r","ZH9_FDSW220008792-2r",
                                           "ST10_L1","ST12_L1","ST13_L1","ST15_L1","ST1_FDME220007465-1r",
                                           "ST3_L1","ST7_FDME220007471-1r","ST9_L1","XM10_FDME220007444-1r",
                                           "XM12_FDME220007446-1r","XM14_FDME220007448-1r","XM15_FDME220007449-1r",
                                           "XM2_FDME220007436-1r","XM4_FDME220007438-1r","XM6_FDME220007440-1r",
                                           "XM8_FDME220007442-1r","WZ10_L1","WZ11_L1","WZ14_L1",
                                           "WZ15_FDME210455471-1r","WZ1_L1","WZ3_L1","WZ5_L1",
                                           "WZ7_FDME210455463-1r","NB11_FDME210455362-1b","NB12_FDME210455363-1r",
                                           "NB15_FDME210455366-1r","NB1_FDME210455352-1r","NB2_FDME210455353-1r",
                                           "NB4_FDME210455355-1r","NB6_FDME210455357-1r","NB7_FDME210455358-1r",
                                           "YC11_FDME220007415-1r","YC13_FDME220007417-1r","YC15_FDME220007419-1r",
                                           "YC1_FDME220007405-1r","YC3_FDME220007407-1r","YC5_FDME220007409-1r",
                                           "YC7_FDME220007411-1a","YC9_FDME220007413-1r","LYG11_FDME210455392-2r",
                                           "LYG13_FDME210455394-1a","LYG15_FDME210455396-1a","LYG1_FDME210455382-1r",
                                           "LYG3_FDME210455384-1r","LYG5_FDME210455386-1r","LYG6_FDME210455387-1r",
                                           "LYG7_FDME210455388-1r","CXD_1412_FDSW210423229-1r","CXD_1414_FDSW210423231-1r",
                                           "CXD_1604_FDSW210423236-1r","CXD_1609_FDSW210423241-1r","CXD_1610_FDSW210423242-1r",
                                           "CXD_1611_FDSW210423243-1r","CXD_1612_FDSW210423244-1r","CXD_1613_FDSW210423245-1r",
                                           "DY10_FDME210438608-1r","DY11_FDME210438609-1r","DY13_FDME210438611-1b","DY15_FDME210438613-1r",
                                           "DY2_FDME210438600-1b","DY3_FDME210438601-1r","DY6_FDME210438604-1r","DY7_FDME210438605-1r",
                                           "DD11_FDSW220008398-1r","DD13_FDSW220008400-2r","DD15_FDSW220008402-2r","DD1_FDSW220008388-1r",
                                           "DD3_FDSW220008390-2r","DD5_FDSW220008392-2r","DD7_FDSW220008394-1r","DD9_FDSW220008396-1r"
))
#Draw the pattern
p2 <- ggplot(data4,aes(x=place,y=precent,fill=kegg_id))+
  geom_bar(stat = "identity",position = "stack",width = 1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),legend.position = "none")+
  scale_fill_manual(values = colors1)
#Display the pattern
p2
#Save the PDF file
ggsave("kegg_contig_ko_flatening_25.5.9.pdf",p2,width = 8.27, height = 11.69)
