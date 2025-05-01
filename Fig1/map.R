setwd("~/Desktop/article1/third_version/")
# 安装和加载必要的包
install.packages("ggplot2")
install.packages("maps")
install.packages("cowplot")
install.packages("rnaturalearth")
install.packages("ggspatial")
#加载相关的包
library(ggplot2)
library(maps)
library(cowplot)
library(rnaturalearth)
library(sf)
library(ggspatial)
#在线读入数据，读入全世界数据
tiny_countries <- ne_countries(type = "sovereignty", scale = 10)
# 获取中国省份数据
china_provinces <- ne_states(country = "China", returnclass = "sf")
# 获取全球海洋数据
oceans <- ne_download(scale = 10, type = "ocean", category = "physical", returnclass = "sf")
#全世界河流数据
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
#十断线数据读取
ceshi2<- "十段线.shp"
ten_line<- st_read(ceshi2)



map_all <- ggplot() +
  geom_sf(data = tiny_countries,color= "#AA9F93",fill="#FBF9F5") +
  geom_sf(data = china_provinces, color = "#44AFDA",fill="#E6F2EE",linewidth = 0.03)+
  geom_sf_text(data = china_provinces, aes(label = name),size = 2, color = "darkred", check_overlap = TRUE) +
  geom_sf_text(data = tiny_countries, aes(label = name), size = 2, color = "darkred", check_overlap = TRUE) +
  geom_sf(data = oceans, fill = "#BCD5EF", color = NA)+
  geom_sf(data = rivers, color = "#BCD5EF", size = 0.001) +
  geom_sf(data=ten_line,color= "black")+
  coord_sf(xlim = c(70, 136), ylim = c(14, 55)) +
  annotation_north_arrow(
    location = "tr",  # 指北针位置：右上角（top-right）
    style = north_arrow_fancy_orienteering(),  # 指北针样式
    pad_x = unit(0.1, "cm"),  # 水平方向的内边距
    pad_y = unit(0.1, "cm"),
    width = unit(1, "cm"),height = unit(1, "cm")
  )+
  theme_minimal() 

map1 <-ggplot() +
  geom_sf(data = tiny_countries,color= "#AA9F93",fill="#FBF9F5") +
  geom_sf(data = china_provinces, color = "#44AFDA",fill="#E6F2EE",linewidth = 0.03)+
  geom_sf_text(data = china_provinces, aes(label = name),size = 2, color = "darkred", check_overlap = TRUE) +
  geom_sf_text(data = tiny_countries, aes(label = name), size = 2, color = "darkred", check_overlap = TRUE) +
  geom_sf(data = oceans, fill = "#BCD5EF", color = NA)+
  geom_sf(data = rivers, color = "#BCD5EF", size = 0.001) +
  geom_sf(data=ten_line,color= "black")+
  coord_sf(xlim = c(105, 120), ylim = c(0,25 )) +
  theme_minimal() 

final_map <- ggdraw()+
  draw_plot(map_all)+
  draw_plot(map1,x=0.7,y=0.15,width = 0.3,height = 0.3)



print(final_map)


ggsave("china_provinces_map_13.pdf", width = 8.27, height = 11.69)

