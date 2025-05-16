setwd("Desktop/article1/宏基因组分析/基于kraken_contig微观多样性/10.Microdiversity/")
data1 <- read.table("snp-ploting",header = T,sep = "\t")
library(trackViewer)
library(ggplot2)
SNP <- data1$pos
yaxis <- c(0,5,10,15,20,25,30,35,40,45,50,55,60)
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


lolliplot(sample.gr, features,jitter="label")

data1<-sample.int(15, length(sample.gr), replace=TRUE)


SNP <- c(10,11,100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)w
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)))
lolliplot(sample.gr, features)

c("#FF8833", "#F9712A", "#DFA32D", 
  "#51C6E6", "#009DDA", "#4B9CDF")