‘’‘
Figure 5 plotting information, some of the data can be obtained from the data sources in the article
’‘’
#Load the corresponding package
library(tidyverse)
library(reshape2)
library(tibble)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggpmisc)
library(RColorBrewer)
#Load the corresponding folder
setwd("~/coverm_count/kegg/")
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
#Filter out outliers using the coverm mean
kegg_data_integer<- read.table("all_kegg_count_reads_clean.txt",header = T,row.names = 1,check.names = F)

#Perform normalization before conducting the analysis
colSums(kegg_data_integer)
kegg_data_Flattening = as.data.frame(t(rrarefy(t(kegg_data_integer), min(colSums(kegg_data_integer)))))
rownames(kegg_data_Flattening) <- row.names(kegg_data_integer)
tkegg_data_Flattening <- t(kegg_data_Flattening)
write.table(file = "all_contig_kegg_raw_flattening.txt",kegg_data_Flattening,quote = F,sep = "\t")
#Calculate richness and merge the data with the corresponding locations
richness_kegg <- specnumber(kegg_data_Flattening,MARGIN = 2)
richness_kegg<- data.frame(richness_kegg)
colnames(richness_kegg) <- "richness"
location_kegg <- read.table('location-ldg.txt',header = T,row.names = 1,check.names = F)
data_kegg_richness_location <- merge(richness_kegg,location_kegg,by = "row.names")
colnames(data_kegg_richness_location)[3] <-c("Latitude")
colnames(data_kegg_richness_location)[4] <-c("location")
data_kegg_richness_location$type <- "function"
#Load the 16S richness-related data and merge it with the corresponding location data
data_16s <- read.table("ASVs_counts_clean_Flattening-rdp-new.txt",header = T,row.names = 1,check.names = F)
richness_16s <- data.frame(specnumber(data_16s,MARGIN = 2))
tdata_16s_Flattening <- t(data_16s)
colnames(richness_16s) <- "richness"
location_16s <- read.table("location-ldg-16s.txt",header = T,row.names = 1)
data_16s_richness_location <- merge(richness_16s,location_16s,by="row.names")
colnames(data_16s_richness_location)[3] <-c("Latitude")
colnames(data_16s_richness_location)[4] <-c("location")
data_16s_richness_location$type <- "taxonomy"
#Integrate the data by merging the 16S and read data.
data_all_richness_location <- rbind(data_16s_richness_location,data_kegg_richness_location)
#Assign colors
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
ggsave("16s_kegg_ldg.pdf",p1,width = 8,height = 6)

#Calculate the DDR
bray_distance_tasv <- vegdist(tkegg_data_Flattening,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
write.table(data2,"kegg_flattening_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#Load the data after processing with the for_ddr Python script.
data_kegg_ddr <- read.table("kegg_flattening_bray_similirty_get_location.txt",header = T)

#Calculate the 16S DDR
bray_distance_tasv <- vegdist(tdata_16s_Flattening,method = "bray")
data1 <- as.matrix(bray_distance_tasv,row.names=1)
data1[upper.tri(data1)]=NA
data2 <- melt(data1,na.rm = T)
write.table(data2,"16s_flattening_bray_similirty.txt",quote = FALSE, row.names = FALSE,sep = "\t")
#Load the data after processing with the for_ddr Python script for 16S DDR calculation
data_16s_ddr <- read.table("16s_flattening_bray_similirty_get_location.txt",header = T)
#Merge the data
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
ggsave("16s_kegg_ddr.pdf",p2,width = 8,height = 6)
