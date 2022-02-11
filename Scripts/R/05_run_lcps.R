library(raster)
library(gdistance)
library(yenpathy)
library(rgdal)
library(sf)
library(magrittr)

# Load top paths function -------------------------------------------------


source(here::here("scripts/R/04_fun_klcp.R"))


# Load Data ---------------------------------------------------------------

hmi <- raster(here::here("data/ProcessedData/hmi_crop_resamp.tif"))
slope <- raster(here::here("data/ProcessedData/slope_crop_resamp.tif"))
imp.fuz.sum1 <- raster(here::here('data/ProcessedData/imp_fuz_sum_nohouse.tif'))
#imp.fuz.sum2 <- raster(here::here('data/ProcessedData/imp_fuz_sum.tif'))
house.val <- raster(here::here("data/ProcessedData/landval_crop_unscaled.tif"))
# Create resistance surfaces  ---------------------------------------------

biophys.resist <- (hmi + 1)^10 + (slope/4) #Dickson et al 2017
implementation.resist1 <- (imp.fuz.sum1 + 1)^10 + (house.val/4) 
implementation.resist1[is.na(implementation.resist1[])] <- 5* cellStats(implementation.resist1, max)## drop NAs for costodistance
# Export resistance surfaces  for use in Circuitscape---------------------------------------------

biophys.resist <- (hmi + 1)^10 + (slope/4) #Dickson et al 2017
implementation.resist1 <- (imp.fuz.sum1 + 1)^10 + (house.val/4) 
writeRaster(biophys.resist, here::here("data/ProcessedData/ResistanceSurfaces/biophys_resist.tif"), overwrite=TRUE)
writeRaster(implementation.resist1, here::here("data/ProcessedData/ResistanceSurfaces/implement_resist.tif"), overwrite=TRUE)

origins <- st_read(here::here("data/ProcessedData/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Weminuche Wilderness" | Unit_Nm == "Yellowstone National Park") %>%
  st_transform(. , crs(biophys.resist))

origins$ID <- c(1,2)

pa.rast <- rasterize(as(origins, "Spatial"), biophys.resist, field="ID")
writeRaster(pa.rast, here::here("data/ProcessedData/pa_locs.tif"), overwrite=TRUE)


# Create Transition Matrix ------------------------------------------------
#need to convert resistance to conductance
biophys.tr <- transition(1/biophys.resist, transitionFunction = mean, 16)
biophys.tr <- geoCorrection(biophys.tr, "c")
saveRDS(biophys.tr, here::here('data/ProcessedData/TransitionLayers/biophystrans.rds'))

social.tr1 <- transition(1/implementation.resist1, transitionFunction = mean, 16)
social.tr1 <- geoCorrection(social.tr1, "c")
saveRDS(social.tr1, here::here('data/ProcessedData/TransitionLayers/socialtrans1.rds'))

#social.tr2 <- transition(1/implementation.resist2, transitionFunction = mean, 16)
#social.tr2 <- geoCorrection(social.tr2, "c")
#saveRDS(social.tr2, here::here('Data/ProcessedData/TransitionLayers/socialtrans2.rds'))
# Estimate k low cost paths -----------------------------------------------

origins <- st_read(here::here("data/ProcessedData/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Weminuche Wilderness") %>% 
  as(. , "Spatial")
goals <- st_read(here::here("data/ProcessedData/studyPAcentroids.shp")) %>% 
  dplyr::filter(., Unit_Nm == "Yellowstone National Park") %>% 
  as(. , "Spatial")
origin.proj <- spTransform(origins, crs(biophys.resist))
goals.proj <- spTransform(goals, crs(biophys.resist))

social1 <- gen_top_paths(tr = social.tr1, resist = implementation.resist1, numpath = 5, bufdist = 4000, orig = origin.proj, goal=goals.proj)

#social2 <- gen_top_paths(tr = social.tr2, resist = implementation.resist2, numpath = 5, bufdist = 4000, orig = origin.proj, goal=goals.proj)

biophys <- gen_top_paths(tr = biophys.tr, resist = biophys.resist, numpath = 5, bufdist = 4000, orig = origin.proj, goal=goals.proj)

saveRDS(social1, here::here('data/ProcessedData/TransitionLayers/socialtop5.rds'))
saveRDS(biophys, here::here('data/ProcessedData/TransitionLayers/biophystop5.rds'))
