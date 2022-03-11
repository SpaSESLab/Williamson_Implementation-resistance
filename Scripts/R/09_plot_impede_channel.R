library(raster)
library(rasterVis)
library(ggplot2)
library(tidyverse)
library(cmocean)
library(viridis)
library(sf)
library(patchwork)
# load cs results ---------------------------------------------------------
implement.cs <- raster(here::here('data/ProcessedData/implement/implement_cum_curmap.asc'))
jurisdiction.cs <- raster(here::here('data/ProcessedData/jurisdiction/jurisdiction_cum_curmap.asc'))
cattle.cs <- raster(here::here("data/ProcessedData/cows/cows_cum_curmap.asc"))
implement.min <- raster(here::here('data/ProcessedData/implementnorm/implementnorm_cum_curmap.asc'))



# normalize ---------------------------------------------------------------

implement.norm <- implement.cs/implement.min
juris.norm <- jurisdiction.cs/implement.min
cattle.norm <- cattle.cs/implement.min

# Convert to df -----------------------------------------------------------

implement.norm.df <- implement.norm %>%
  projectRaster(., res=300, crs = crs(implement.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "implement")) %>% 
  mutate(., channel = ifelse(implement > 1.2, 1, 
                              ifelse(2 > implement & implement >= 0.8, NA, 0)))

juris.norm.df <- juris.norm %>%
  projectRaster(., res=300, crs = crs(implement.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "jurisdiction")) %>% 
  mutate(., channel = ifelse(jurisdiction > 1.2, 1, 
                             ifelse(2 > jurisdiction & jurisdiction >= 0.8, NA, 0)))

cattle.norm.df <- cattle.norm %>%
  projectRaster(., res=300, crs = crs(implement.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "cattle")) %>% 
  mutate(., channel = ifelse(cattle > 1.2, 1, 
                             ifelse(2 > cattle & cattle >= 0.8, NA, 0)))

# Add addtional map eleements ---------------------------------------------

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
elev.proj <- projectRaster(elev[[1]], implement.cs)


elev.crop <- crop(elev.proj, implement.cs)
elev.mod <- elev.crop *10
slope <- terrain(elev.mod, opt='slope')
aspect <- terrain(elev.mod, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)
hill2 <- aggregate(hill , fact = 5 , method = "bilinear" )
hills3 <- focal(hill2, w=matrix(1/9, nc=3, nr=3), mean)

# Load centroids ----------------------------------------------------------
origins <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Weminuche Wilderness") %>% 
  as(. , "Spatial")
goals <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Yellowstone National Park") %>% 
  as(. , "Spatial")
origin.proj <- spTransform(origins, crs(implement.cs))
goals.proj <- spTransform(goals, crs(implement.cs))
# Get vectors for maps ----------------------------------------------------

PAs <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAs.shp")) %>% 
  dplyr::filter(. , Unit_Nm == "Yellowstone National Park" | Unit_Nm == "Weminuche Wilderness") %>% 
  st_transform(. , crs = crs(implement.cs))

sts <- tigris::states() %>% 
  dplyr::filter(., STUSPS %in% c("ID", "MT", "UT", "WY", "CO", "AZ", "NM")) %>% 
  st_transform(., crs(implement.cs)) %>% 
  as(., "Spatial")

conus <-  tigris::states() %>%
  dplyr::filter(,, !STUSPS %in% c("AK","AS", "HI", "GU", "VI", "DC", "PR", "MP")) %>% 
  st_transform(., crs(implement.cs))


sts.crop <- crop(sts, implement.cs)

pa.cents <- rbind(as(origin.proj,"sf"), as(goals.proj, "sf"))
pa.cents$lab <- c("Weminche \n Wilderness Area", "Yellowstone \n National Park")


# Make the plots ----------------------------------------------------------

p1 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = implement.norm.df[!is.na(implement.norm.df$channel),],
    aes(
      x=x,
      y=y,
      fill = factor(channel)
    ),
    interpolate = TRUE,
  ) +
  scale_fill_viridis(labels = c("Impeded", "Channelized"), discrete = TRUE, begin = 0.33, end = 0.8,option = "D") +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), 
                           nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "white")+
  guides(fill = guide_legend(title.position = "top", 
                               label.position="bottom", title = NULL)) +
  theme_map() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = "center"  )

p2 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = juris.norm.df[!is.na(juris.norm.df$channel),],
    aes(
      x=x,
      y=y,
      fill = factor(channel)
    ),
    interpolate = TRUE,
  ) +
  scale_fill_viridis(labels = c("Impeded", "Channelized"), discrete = TRUE, begin = 0.33, end = 0.8,option = "D") +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), 
                           nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "white")+
  guides(fill = guide_legend(title.position = "top", 
                             label.position="bottom", title = NULL)) +
  theme_map() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = "center"  )


p3 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = cattle.norm.df[!is.na(cattle.norm.df$channel),],
    aes(
      x=x,
      y=y,
      fill = factor(channel)
    ),
    interpolate = TRUE,
  ) +
  scale_fill_viridis(labels = c("Impeded", "Channelized"), discrete = TRUE, begin = 0.33, end = 0.8,option = "D") +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), 
                           nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "white")+
  guides(fill = guide_legend(title.position = "top", 
                             label.position="bottom", title = NULL)) +
  theme_map() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = "center"  )

combined.plot <- p1+p2+p3 + plot_layout(guides = "collect") +  plot_annotation(tag_levels = "A", tag_suffix = ")") & theme(legend.position = "bottom", legend.justification = "center",
                                                                                                         legend.margin=margin(0,0,0,0),
                                                                                                         legend.box.margin=margin(-15,-10,-10,-10)) 
