#Load the directory where the file is located
setwd("~/Desktop/article1/third_version/")
#Install necessary packages (if not already installed)
install.packages("ggplot2")
install.packages("maps")
install.packages("cowplot")
install.packages("rnaturalearth")
install.packages("ggspatial")
#Load the packages
library(ggplot2)
library(maps)
library(cowplot)
library(rnaturalearth)
library(sf)
library(ggspatial)
#Read data online, load global data
tiny_countries <- ne_countries(type = "sovereignty", scale = 10)
#Retrieve data for Chinese provinces
china_provinces <- ne_states(country = "China", returnclass = "sf")
#Retrieve global ocean data
oceans <- ne_download(scale = 10, type = "ocean", category = "physical", returnclass = "sf")
#Global river data
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
#Ten-dash line data retrieval
ceshi2<- "the-ten-dasg-line.shp"
ten_line<- st_read(ceshi2)
#Generate the overall plot
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
    location = "tr",  # Compass position: top-right
    style = north_arrow_fancy_orienteering(),  # Compass style
    pad_x = unit(0.1, "cm"),  # Horizontal padding
    pad_y = unit(0.1, "cm"),
    width = unit(1, "cm"),height = unit(1, "cm")
  )+
  theme_minimal() 
#Plot the ten-dash line data
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
#Merge the map of the ten-dash line with the overall map
final_map <- ggdraw()+
  draw_plot(map_all)+
  draw_plot(map1,x=0.7,y=0.15,width = 0.3,height = 0.3)
#Show the pattern
print(final_map)
#Save the PDF file"
ggsave("china_provinces_map.pdf", width = 8.27, height = 11.69)

