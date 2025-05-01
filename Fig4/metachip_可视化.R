# 加载所需的包
library(ggforce)
library(tidyr)
library(ggplot2)

# 创建示例数据
data <- data.frame(
  Source = c("古菌1", "古菌1", "古菌2"),
  Target = c("古菌2", "细菌1", "细菌1"),
  Type = c("得到", "给予", "得到"),
  Value = c(10, 5, 8)
)

# 整理数据为适合平行集合图的格式
data1 <- gather_set_data(data, 1:2)  # 1:2 表示 Source 和 Target 列

# 绘制平行集合图
ggplot(data1, aes(x = x, id = id, split = y, value = Value)) +
  geom_parallel_sets(aes(fill = Type), alpha = 0.6, axis.width = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.3, colour = "black", fill = "white") +
  scale_x_discrete(labels = c("Source", "Target")) +  # 设置 x 轴标签
  labs(
    title = "水平基因转移可视化",
    x = "分类",
    y = "基因转移数量",
    fill = "类型"
  ) +
  theme_minimal()


setwd('~/Desktop/mag-hgt/')
head(as.data.frame(UCBAdmissions))
rm(list = ls())

data2 <- read.table("hgt_ceshi_1",header = T)
data3 <- gather_set_data(data2,3:5)


color_list <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
  "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
  "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
  "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#6b6ecf","green","blue"
)


ggplot(data3, aes(x = x, id = id, split = y, value = Value)) +
  geom_parallel_sets(aes(fill = Target), alpha = 0.6, axis.width = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.3, colour = "black",aes(fill = data3$y))+
  geom_parallel_sets_labels(colour="black",angle = 0)+
  scale_fill_manual(values = color_list)+
  theme_minimal()

