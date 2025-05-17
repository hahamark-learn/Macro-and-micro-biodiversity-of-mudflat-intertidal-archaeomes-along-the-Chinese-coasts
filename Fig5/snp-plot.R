‘’‘
Plot data for a single SNP
’‘’
setwd("~/Microdiversity/")
data1 <- read.table("snp-ploting",header = T,sep = "\t")
library(trackViewer)
library(ggplot2)
SNP <- data1$pos
yaxis <- c(0,5,10,15,20,25,30,35,40,45,50,55,60)#The numbers here represent the differences between different locations
xaxis <- c(1,400,800,1200,1400,1573)

sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)))
sample.gr$color <- data1$pos_in_codon
sample.gr$score <- data1$location
features <- GRanges("chr1", IRanges(c(1), 
                                    width=c(1434),
                                    names=paste0("CDS")))
features$fill <- c("mistyrose")
features$height <- c(0.05)
lolliplot(sample.gr, features,yaxis = yaxis,xaxis = xaxis,jitter=c("node", "label"),legendPosition = 'right',lollipop_style_switch_limit=5)

ggsave(file="snp-ceshi6.pdf",p1,width = 20,height = 10)
