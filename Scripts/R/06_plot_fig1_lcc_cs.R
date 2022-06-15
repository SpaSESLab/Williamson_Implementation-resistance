library(raster)
library(gdistance)
library(rgdal)
library(sf)
library(magrittr)
library(ggplot2)
library(purrr)
library(dplyr)
library(rasterVis)
library(patchwork)
library(classInt)
library(cowplot)
library(ggsci)
library(stringr)
library(egg)
# load Circuitscape results -----------------------------------------------

biophys.cs <- raster(here::here('data/ProcessedData/biophys/biophys_cum_curmap.asc'))


# Load centroids ----------------------------------------------------------
origins <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Weminuche Wilderness") %>% 
  as(. , "Spatial")
goals <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Yellowstone National Park") %>% 
  as(. , "Spatial")
origin.proj <- spTransform(origins, crs(biophys.cs))
goals.proj <- spTransform(goals, crs(biophys.cs))


# Load k cost paths -------------------------------------------------------

biophys <- readRDS(here::here('Data/ProcessedData/TransitionLayers/biophystop5.rds'))
biophys.lst <- biophys[[1]]


# Custom map theme --------------------------------------------------------

theme_map <- function(...) {
  theme_minimal() +
    theme(
      #text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      axis.line = element_line(color = "black", size=0.02),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "white", size = 0.002),
      panel.grid.minor = element_line(color = "white", size = 0.002),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.background = element_rect(fill = "white", color = NA), 
      legend.background = element_rect(fill = "white", color = NA),
      #panel.border = element_rect(fill=NA, color = "black"),
      ...
    )
}
# Load Elev data and calc hillshade ---------------------------------------
elev <- getData('alt', country = 'USA')
elev.proj <- projectRaster(elev[[1]], biophys.cs)


elev.crop <- crop(elev.proj, biophys.cs)
elev.mod <- elev.crop *10
slope <- terrain(elev.mod, opt='slope')
aspect <- terrain(elev.mod, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)
hill2 <- aggregate(hill , fact = 5 , method = "bilinear" )
hills3 <- focal(hill2, w=matrix(1/9, nc=3, nr=3), mean)


# Get vectors for maps ----------------------------------------------------

PAs <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAs.shp")) %>% 
  dplyr::filter(. , Unit_Nm == "Yellowstone National Park" | Unit_Nm == "Weminuche Wilderness") %>% 
  st_transform(. , crs = crs(biophys.cs))

sts <- tigris::states() %>% 
  dplyr::filter(., STUSPS %in% c("ID", "MT", "UT", "WY", "CO", "AZ", "NM")) %>% 
  st_transform(., crs(biophys.cs)) %>% 
  as(., "Spatial")

conus <-  tigris::states() %>%
  dplyr::filter(,, !STUSPS %in% c("AK","AS", "HI", "GU", "VI", "DC", "PR", "MP")) %>% 
  st_transform(., crs(biophys.cs))


sts.crop <- crop(sts, biophys.cs)

pa.cents <- rbind(as(origin.proj,"sf"), as(goals.proj, "sf"))
pa.cents$lab <- c("Weminche \n Wilderness Area", "Yellowstone \n National Park")
# Plot Circuitscape results -----------------------------------------------


b_df <- biophys.cs %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "biophys"))

biophys.lcps <- lapply(biophys.lst, function (x) rasterToPolygons(x=x, n=8, dissolve = TRUE))

biophys.lcp.sf <- list(biophys.lcps, makeUniqueIDs = T) %>% 
  flatten() %>% 
  do.call(rbind, .) %>% 
  as(., "sf") %>% 
  dplyr::mutate(., rank = row_number())

my_blues = RColorBrewer::brewer.pal(n = 9, "Blues")[4:8]

p.cs <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = b_df,
    aes(
      x=x,
      y=y,
      fill=biophys
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_viridis_c(option="B", 
                       limits=c(-0.000000001,0.005), 
                       breaks = c(0,0.005),
                       labels = c("Low","High")) +
  geom_sf(data=biophys.lcp.sf, aes(color=as.character(rank)), fill=NA, lwd=1.15) +
  scale_colour_manual(values=my_blues) +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="white")+
  ggrepel::geom_text_repel(data = pa.cents, size = 1.75, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "white")+
  guides(fill = guide_colorbar(nbin = 8,  title.position = "left", 
                               label.position="right", title = "Current flow", barheight = 2.25, barwidth = 1), 
         color = guide_legend(title.position = "left", 
                              label.position="right", title = "Cost rank", keywidth = 1, keyheight = 1))+
  theme_map() + theme(legend.position = "left", 
                      legend.box = "vertical", 
                      legend.title = element_text(angle = 90, size = 3.75),
                      legend.text = element_text(size = 3.5),
                      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "pt"),
                      legend.box.margin = unit(c(0, -0.1, 0.1, 0.1), units = "pt"))

inset <- ggplot()+
  geom_sf(data=conus, fill="white") +
  geom_sf(data =st_as_sfc(st_bbox(biophys.cs)), fill=NA, color="red") +
  theme_map() +
  theme(panel.background = element_rect(fill = "gray", color = "black"),
        plot.background = element_rect(fill = NA, color = NA))  


p.combined <- ggdraw(p.cs) +
  draw_plot(inset, x = 0.33, y = 0,  
            width = 0.3, height = 0.15)



cowplot::save_plot(here::here("plots/fig1.png"), plot =p.combined, base_height = 4)


