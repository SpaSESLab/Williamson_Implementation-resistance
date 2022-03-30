library(raster)
library(gdistance)
library(rgdal)
library(sf)
library(magrittr)
library(ggplot2)
library(rasterVis)
library(patchwork)
library(classInt)
#library(cowplot)
library(ggsci)
library(stringr)
library(egg)
# Load transition layers --------------------------------------------------
social.tr1 <- readRDS(here::here('data/ProcessedData/TransitionLayers/socialtrans1.rds'))
biophys.tr <- readRDS(here::here('data/ProcessedData/TransitionLayers/biophystrans.rds'))
jurisdiction.tr <- readRDS(here::here('data/ProcessedData/TransitionLayers/socialtrans_jurisdiction.rds'))
cattle.tr <- readRDS(here::here('data/ProcessedData/TransitionLayers/socialtrans_cattle.rds'))

# Load Resistance surfaces ------------------------------------------------
biophys.resist <- raster(here::here("data/ProcessedData/ResistanceSurfaces/biophys_resist.tif"))
implementation.resist1 <- raster(here::here("data/ProcessedData/ResistanceSurfaces/implement_resist.tif"))
implementation.resist2 <-  raster(here::here('data/ProcessedData/ResistanceSurfaces/implement_resist_jurisdiction.tif'))
implementation.resist3 <-  raster(here::here('data/ProcessedData/ResistanceSurfaces/implement_resist_cattle.tif'))
# Load centroids ----------------------------------------------------------
origins <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Weminuche Wilderness") %>% 
  as(. , "Spatial")
goals <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Yellowstone National Park") %>% 
  as(. , "Spatial")
origin.proj <- spTransform(origins, crs(biophys.resist))
goals.proj <- spTransform(goals, crs(biophys.resist))


# Load k cost paths -------------------------------------------------------

social1 <- readRDS(here::here('Data/ProcessedData/TransitionLayers/socialtop5.rds'))
biophys <- readRDS(here::here('Data/ProcessedData/TransitionLayers/biophystop5.rds'))

# Load biophys circuitscape run -------------------------------------------

biophys.cs <- raster(here::here("data/ProcessedData/biophys/biophys_cum_curmap.asc"))

# Create Euclidean distance-----------------------------------------------------

euclidean.resist <- biophys.resist
values(euclidean.resist) <- 1
euclidean.tr <- transition(1/euclidean.resist, transitionFunction = mean, 16)
euclidean.tr <- geoCorrection(euclidean.tr, "c")
eucdist <- accCost(euclidean.tr, origin.proj)

#extract values for last 50km
euc.dist <- st_distance(as(origin.proj, "sf"), as(goals.proj, "sf"))
units(euc.dist) <- NULL
threshold <- as.vector(euc.dist - 10000)
# Extract cost and current values -----------------------------------------


social1.lst <- social1[[1]]
#social2.lst <- social2[[1]]
biophys.lst <- biophys[[1]]
#all.lst <- c(social1.lst, social2.lst, biophys.lst)
all.lst <- c(social1.lst, biophys.lst)
soc.costdist1 <- accCost(social.tr1, origin.proj)
jurisdiction.costdist <- accCost(jurisdiction.tr, origin.proj)
cattle.costdist <- accCost(cattle.tr, origin.proj)
bio.costdist <- accCost(biophys.tr, origin.proj)
#dist.stack <- stack(eucdist,bio.costdist, soc.costdist1, soc.costdist2)
dist.stack <- stack(eucdist,bio.costdist, soc.costdist1, jurisdiction.costdist, cattle.costdist, biophys.cs)

distance.extract <- lapply(1:length(all.lst), function(x) raster::extract(dist.stack, rasterToPolygons(all.lst[[x]], dissolve = TRUE)))
#names(distance.extract) <- c("s1_1", "s1_2", "s1_3", "s1_4", "s1_5","s2_1", "s2_2", "s2_3", "s2_4", "s2_5","b1", "b2", "b3", "b4", "b5")
names(distance.extract) <- c("s1", "s2", "s3", "s4", "s5","b1", "b2", "b3", "b4", "b5")

dist.df.list <- lapply(1:length(distance.extract), function(x) as.data.frame(distance.extract[[x]]))
names(dist.df.list) <- names(distance.extract)

dist.df <- dplyr::bind_rows(dist.df.list, .id="column_label")
#colnames(dist.df) <- c("LCPID","eucdist","bio.costdist", "soc.costdist1", "soc.costdist2")
colnames(dist.df) <- c("LCPID","eucdist","bio.costdist", "soc.costdist1", "jurisdiction.costdist","cattle.costdist", "biophys.cur")

dist.cln <- do.call(data.frame,lapply(dist.df, function(x) replace(x, is.infinite(x),NA)))

# Estimate current density ------------------------------------------------

current.dens <- dist.cln %>% 
  dplyr::group_by(., LCPID) %>%
  dplyr::filter(., startsWith(LCPID, "b")) %>% 
  dplyr::summarise(., sum.current = round(sum(biophys.cur, na.rm = TRUE), 2),
                   max.dist = round(max(eucdist, na.rm = TRUE)/1000, 2))  %>%
  dplyr::mutate(., ID = str_remove_all( LCPID, fixed("b"))) %>% 
  dplyr::select(., c(ID, max.dist, sum.current)) %>% 
  flextable::flextable() %>% 
  flextable::set_header_labels(., ID = "Cost Rank", max.dist = "Distance (km)", sum.current = "Total Current") 

flextable::save_as_docx(current.dens, path = here::here("plots/table1.docx"))

# Estimate cost.ratio -----------------------------------------------------
library(tidyr)
library(dplyr)
dist.90 <- dist.cln %>% 
  dplyr::group_by(LCPID) %>% 
  dplyr::filter(., eucdist > threshold)
#dist.90[,c(2:5)] <- apply(dist.90[,c(2:5)], 2, log10)
#dist.90[,c(2:4)] <- apply(dist.90[,c(2:4)], 2, log10)

max.bphys.lcp <- dist.90 %>% 
  dplyr::group_by(., LCPID) %>%
  dplyr::summarise(., max_biophys = max(bio.costdist, na.rm=TRUE))
max.difs <- max.bphys.lcp
max.difs[-1] <- lapply(max.bphys.lcp[-1], function(x) x - min(x, na.rm=TRUE))
colnames(max.difs)[2] <- "delta_bio.costdist"

dist.group <- dist.90 %>% 
  dplyr::left_join(., max.difs) %>%
  dplyr::mutate(., abase = soc.costdist1/bio.costdist,
                ajuris = jurisdiction.costdist/bio.costdist,
                bcattle = cattle.costdist/bio.costdist) %>%
  dplyr::filter(.,  startsWith(LCPID, "b")) %>% 
  dplyr::select(., c(1,8:11)) %>% 
  pivot_longer(., !c(LCPID, delta_bio.costdist)) %>% 
  dplyr::mutate(., ID = str_remove_all( LCPID, fixed("b")))

dist.sum <- dist.group %>% 
  group_by(., ID, name) %>% 
  summarise(mn = mean(value),
            std = sd(value),
            costdist = log(delta_bio.costdist, 10))


p1 <- ggplot(data = filter(dist.sum, name == "abase"), mapping = aes(x = ID, 
                                                                     y = mn, 
                                                                     color = costdist)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_viridis_c(direction = -1, begin= 0, end =0.8) + 
  ylim(110, 150) +
  labs(x= "Cost Rank", y= "Cost ratio") + 
  guides(color = guide_colorbar(title = "\u0394 log(biophysical \ncost)")) +
  theme_bw()

p2 <- ggplot(data = filter(dist.sum, name == "ajuris"), mapping = aes(x = ID, 
                                                                      y = mn, 
                                                                      color = costdist)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_viridis_c(direction = -1, begin= 0, end =0.8) + 
  ylim(110, 150) +
  labs(x= "Cost Rank", y= "Cost ratio") + 
  guides(color = guide_colorbar(title = "\u0394 log(biophysical \ncost)")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p3 <- ggplot(data = filter(dist.sum, name == "bcattle"), mapping = aes(x = ID, 
                                                                       y = mn, 
                                                                       color = costdist)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_viridis_c(direction = -1, begin= 0, end =0.8) + 
  ylim(110, 150) +
  labs(x= "Cost Rank", y= "Cost ratio") + 
  guides(color = guide_colorbar(title = "\u0394 log(biophysical \ncost)")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )

combined <- p1 + p2 + p3 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A", tag_suffix = ")")

ggsave('plots/fig2.png', combined)
