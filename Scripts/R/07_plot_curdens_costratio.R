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

eucdist <- distanceFromPoints(biophys.resist, origin.proj)
#extract values for last 10km
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


# Estimate current density ------------------------------------------------

current.dens <- dist.df %>% 
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
dist.90 <- dist.df %>% 
  dplyr::group_by(LCPID) %>% 
  dplyr::filter(., eucdist > threshold)
#dist.90[,c(2:5)] <- apply(dist.90[,c(2:5)], 2, log10)
#dist.90[,c(2:4)] <- apply(dist.90[,c(2:4)], 2, log10)

max.biophys <- max(dist.90$bio.costdist)
max.social <- max(dist.90$soc.costdist1)
max.cattle <- max(dist.90$cattle.costdist)
max.juris <- max(dist.90$jurisdiction.costdist)

dist.tbl <- dist.90 %>% 
  mutate(., biocostkm = bio.costdist/eucdist,
               soccostkm = soc.costdist1/eucdist,
               juriscostkm = jurisdiction.costdist/eucdist,
               cattlecostkm = cattle.costdist/eucdist,
               tobase = (1+(max.biophys - bio.costdist)/(1+(max.social-soc.costdist1))),
               tocattle = (1+(max.biophys - bio.costdist)/(1+(max.cattle - cattle.costdist))),
               tojuris = (1+(max.biophys - bio.costdist)/(1+(max.juris - jurisdiction.costdist)))
            ) %>%
  select(., c(1, 8:14)) %>% 
  dplyr::filter(.,  startsWith(LCPID, "b")) %>% 
  pivot_longer(., !LCPID)%>% 
  dplyr::mutate(., ID = str_remove_all( LCPID, fixed("b")))
  
dist.sum <- dist.tbl %>% 
  group_by(., ID, name) %>% 
  summarise(mn = mean(value),
            std = sd(value))


my_blues = RColorBrewer::brewer.pal(n = 9, "Blues")[4:8]
p1 <- ggplot(data = filter(dist.sum, name == "biocostkm"), mapping = aes(x = ID, 
                                                                     y = mn, 
                                                                     color = ID)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_manual(values=my_blues) + 
  labs(x= "Cost Rank", y= "Cost per km") + 
  guides(color = guide_legend(title = "Cost Rank", direction="horizontal")) +
  ggtitle("Biophysical costs") +
  theme_bw(base_size = 12, base_family = "Times" ) + theme(legend.position = "bottom", axis.text = element_text(size=13), axis.title = element_text(size=15))

p2 <- ggplot(data = filter(dist.sum, name == "soccostkm"), mapping = aes(x = ID, 
                                                                         y = mn, 
                                                                         color = ID)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_manual(values=my_blues) + 
  ylim(c(470, 630)) +
  labs(x= "Cost Rank", y= NULL) + 
  guides(color = guide_legend(title = "Cost Rank", direction="horizontal")) +
  ggtitle("Implementation costs\n (Baseline)") +
  theme_bw(base_size = 12, base_family = "Times") + theme(legend.position = "bottom", axis.text = element_text(size=13), axis.title = element_text(size=15))


p3 <- ggplot(data = filter(dist.sum, name == "juriscostkm"), mapping = aes(x = ID, 
                                                                         y = mn, 
                                                                         color = ID)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_manual(values=my_blues) + 
  ylim(c(470, 630)) +
  labs(x= "Cost Rank", y= NULL) + 
  guides(color = guide_legend(title = "Cost Rank", direction="horizontal")) +
  ggtitle("Implementation costs \n(policy scenario)") +
  theme_bw(base_size = 12, base_family = "Times" )+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=13), axis.title.x = element_text(size=15), legend.position = "bottom")

p4 <-ggplot(data = filter(dist.sum, name == "cattlecostkm"), mapping = aes(x = ID, 
                                                                          y = mn, 
                                                                          color = ID)) +
  geom_pointrange(mapping = aes(ymin = mn - std,
                                ymax = mn + std)) +
  scale_color_manual(values=my_blues) + 
  ylim(c(470, 630)) +
  labs(x= "Cost Rank", y= NULL) + 
  guides(color = guide_legend(title = "Cost Rank", direction="horizontal")) +
  ggtitle("Implementation costs \n(land use scenario)") +
  theme_bw(base_size = 12, base_family = "Times" ) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=13), axis.title.x = element_text(size=15), legend.position = "bottom")



combined <- (p1 | p2 | p3  | p4) + plot_annotation(tag_levels = "A", tag_suffix = ")") + 
  plot_layout(guides='collect') &
  theme(legend.position='bottom')

ggsave('plots/fig2.png', combined, width = 9.5, height = 7.5, units = "in")

to <- dist.tbl %>% 
  group_by(., ID, name) %>% 
  summarise(mn = mean(value),
            std = sd(value)) %>% 
  filter(., str_detect(name, "^to")) %>% 
  pivot_wider(., id_cols = ID, names_from = name, values_from=c(mn,std)) %>% 
  select(., 1:4) %>% 
  flextable::flextable() %>% 
  flextable::set_header_labels(., ID = "Cost Rank", mn_tobase = "Trade-off Index (baseline)", mn_tojuris = "Trade-off Index (policy scenario)", mn_tocattle = "Trade-off Index (land use scenario)") 

flextable::save_as_docx(to, path = here::here("plots/table1.docx"))
