‘’‘
Here are the plotting data and scripts related to Figure 4
’‘’
#Load the corresponding package
library(reshape2)
library(dplyr)
library(tidyr)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(cowplot)

#Read in the relevant data
#Load the relevant HGT data and carry out necessary processing
setwd('~/Desktop/mag-hgt/')
data2 <- read.table("metachip_hgt",header = T)
data3 <- gather_set_data(data2,3:5)
#Set color list.
color_list <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
  "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
  "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
  "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#6b6ecf","green","blue"
)

#Draw the pattern
ggplot(data3, aes(x = x, id = id, split = y, value = Value)) +
  geom_parallel_sets(aes(fill = Target), alpha = 0.6, axis.width = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.3, colour = "black",aes(fill = data3$y))+
  geom_parallel_sets_labels(colour="black",angle = 0)+
  scale_fill_manual(values = color_list)+
  theme_minimal()


#Read the relevant COG data
data1 <- read.table('hgt_cog_number',header = T)
#hgtector visualization analysis
#Thaumarchaeota
data2 <- read.table('thaumarchaeota_cog_number_10',header = T)
#Euryarchaeota
data3 <- read.table('euryarchaeota_cog_number_10',header = T)
#Nitrososphaerota
data4 <- read.table('Nitrososphaerota_raw_number_10',header = T)
#metachip_archaea_bacteria
data5 <- read.table("ATOB_hgt_count_clean",sep = " ",header = T)
#metachip_archaea_to_archaea
data6 <- read.table("ATOA_cog_clean",sep = " ",header = T)
#Merge all data
data2$phylum <- "Thaumarchaeota"
data3$phylum <- "Euryarchaeota"
data4$phylum <- "Nitrososphaerota"
hgtector_all<- rbind(data2,data3,data4)
data5$phylum <- "AB"
data6$phylum <- "AA"
metachip_all<- rbind(data5,data6)
hgtector_all$type <- "hgtector"
metachip_all$type <- "metachip"
hgt_all<- rbind(hgtector_all,metachip_all)
#Set gradient colors
color1 <- brewer.pal(10,"Paired")
color2 <- colorRampPalette(color1)(19)
#Create a faceted plot
p1 <- ggplot(hgt_all, aes(x = phylum, y = group, size = value, color = group)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.15)) +  #Add jitter to avoid overlapping points
  scale_size_continuous(range = c(3, 18))+
  labs(x = "Category", y = "Value", size = "Size", color = "Color") +  # Set axis labels and legend titles
  facet_grid(~type,scales = "free_x")+
  scale_color_manual(values=color2)
#Save the PDF file
ggsave("hgt_plot_25.3.11.pdf", plot = p1, width = 210 / 25.4, height = 297 / 25.4)

