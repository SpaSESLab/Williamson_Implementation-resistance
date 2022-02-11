library(sf)
library(rgdal)
library(tidyverse)
library(tigris)
library(gdalUtilities)
library(raster)
library(terra)

# Function to fix PADUS geometries ----------------------------------------

ensure_multipolygons <- function(X) {
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geom = st_geometry(Y))
}
ogrListLayers(here::here("Data/OriginalData/PAD_US2_1.gdb/"))

# generate PADUS datasets -------------------------------------------------

state.nm <- c("ID", "MT", "UT", "WY", "CO", "AZ", "NM")
public.land <- c("FED","TRIB","STAT","DIST","LOC")


public.land.region <- read_sf(dsn = here::here("Data/OriginalData/PAD_US2_1.gdb/"), layer = "PADUS2_1Combined_Proclamation_Marine_Fee_Designation_Easement") %>% 
  filter(., State_Nm %in% state.nm & Mang_Type %in% public.land)


public.land.corrected <-  ensure_multipolygons(public.land.region)
a <- which(is.na(st_is_valid(public.land.corrected)))

public.land.fail <- public.land.corrected[a,]

#figure out how many pts is in each feature that throws an NA in st_is_valid
wrong.pts <- lapply(1:nrow(public.land.fail), function(x) sapply(st_geometry(public.land.fail)[[x]], function(y) nrow(y[[1]])))
#make sure the number of pts is the reason
few.pts <- lapply(1:length(wrong.pts), function(x) any(wrong.pts[[x]] < 4))

g <- st_geometry(public.land.corrected)
g[[a]][1] = NULL #this is the feature with less than 4 pts
st_geometry(public.land.corrected) = g


public.land.valid <- st_make_valid(public.land.corrected)

st_write(public.land.valid[,c(1:13, 23)], here::here('data/ProcessedData/shapefiles/studypubliclands.shp'), append=FALSE)

fed.land <- public.land.valid %>% 
  filter(., Mang_Type=="FED")


st_write(fed.land[,c(1:13, 23)], here::here('data/ProcessedData/shapefiles/fedlands.shp'), append=FALSE)

unit.nm <- c("Yellowstone National Park", "Weminuche Wilderness")

study.areas <- public.land.corrected %>% 
  filter(., Unit_Nm %in% unit.nm) %>% 
  aggregate(., by = list(.$Unit_Nm), FUN = dplyr::first)

st_write(study.areas[,c(1:13, 23)], here::here('data/ProcessedData/shapefiles/studyPAs.shp'), append = FALSE)

study.area.point <- st_point_on_surface(study.areas) #generates a strange point for FCRNR
study.area.centroid <- st_centroid(study.areas)
st_write(study.area.centroid[,c(1:13, 23)], here::here('data/ProcessedData/shapefiles/studyPAcentroids.shp'), append = FALSE)


# Crop rasters to study area and resample to 1km ----------------------------------------------
house.val <- raster::raster(here::here("data/OriginalData/LandValue/places_fmv_all.tif"))
study.pas <- st_read( here::here('data/ProcessedData/shapefiles/studyPAs.shp'))

padus.buf <-  study.pas %>%  
  st_buffer(., dist = 100000) %>% 
  st_transform(., proj4string(house.val)) %>% 
  as(., "Spatial")

house.val.crop <- raster::crop(house.val, padus.buf)

s <- raster(ext=extent(house.val.crop), res=1000, crs=proj4string(house.val.crop))
s.rast <- terra::rast(s)
house.rast <- terra::rast(house.val.crop)
house.resample <- resample(house.rast, s.rast, method="bilinear")
house.resamp <- raster(house.resample)
house.res.test <- house.resamp 

#set public land cost to 0
public.land <- st_read(here::here('Data/ProcessedData/shapefiles/studypubliclands.shp')) %>% 
  st_transform(., proj4string(house.resamp)) %>% 
  as(., "Spatial")


house.res.test[public.land] <- 1 #set public land costs to 0 because value is log

housing.rescale <- (house.res.test - cellStats(house.res.test, min))/(cellStats(house.res.test, max) - cellStats(house.res.test, min))

raster::writeRaster(housing.rescale, here::here("data/ProcessedData/rasters/landval_crop_resamp.tif"), overwrite=TRUE)
raster::writeRaster(house.res.test, here::here("data/ProcessedData/rasters/landval_crop_unscaled.tif"), overwrite=TRUE)

# Crop HMI to study area and resample to 1km--------------------------------------------------

hmi <- raster::raster(here::here("Data/OriginalData/HMI/hm_fsum3_270/"))

padus.buf <- as(padus.buf, "sf") %>% 
  st_transform(., crs(hmi)) %>% 
  as(., "Spatial")

hmi.crop <- crop(hmi, padus.buf)
hmi.rast <- terra::rast(hmi.crop)
hmi.resample <- resample(hmi.rast, terra::rast(housing.rescale), method="bilinear")
hmi.resamp <- raster(hmi.resample)
writeRaster(hmi.resamp, here::here('data/ProcessedData/rasters/hmi_crop_resamp.tif'), overwrite=TRUE)

# Crop slope to study area and resample to 1km ----------------------------

alt<- getData('alt', country='USA', path=here::here("data/OriginalData/elev/"))[[1]]
slope <- terrain(alt, opt='slope', unit = 'degrees')
slope.rast <- terra::rast(slope)
slope.proj <- terra::project(slope.rast, crs(hmi.resample))

padus.buf <- as(padus.buf, "sf") %>% 
  st_transform(., crs(slope.proj)) %>% 
  as(., "Spatial")


slope.crop <- terra::crop(slope.proj, padus.buf)
slope.resample <- terra::resample(slope.crop, terra::rast(housing.rescale))
slope.resamp <- raster(slope.resample)
writeRaster(slope.resamp, here::here('data/ProcessedData/rasters/slope_crop_resamp.tif'), overwrite=TRUE)

# Crop wolf prefs to study area and resample to 1km ----------------------------
wolf.increase <- rast(here::here('data/OriginalData/wolves/wolf.increase.tif'))

padus.buf <-  as(padus.buf, "sf") %>% 
  st_transform(., crs(wolf.increase)) %>% 
  as(., "Spatial")

wolf.inc.crop <- crop(wolf.increase, padus.buf)
wolf.inc.proj <- terra::resample(wolf.inc.crop, terra::rast(housing.rescale), method="bilinear")

wolf.inc.focal <- focal(wolf.inc.proj, w = 9, fun = "mean", na.rm = TRUE, na.policy = "only") #smoothing things out a bit
wolf.resist <-1 - wolf.inc.focal

writeRaster(wolf.resist, here::here('data/ProcessedData/rasters/wolf_crop_resamp.tif'), overwrite=TRUE)


# Cattle numbers ----------------------------------------------------------

cattle.tab <- read.csv("data/OriginalData/cattle_inventory2017.csv")
cty <- tigris::counties(unique(cattle.tab$State), year=2017, refresh = TRUE)
#Need to padd the FIPS codes
cattle.pad <- cattle.tab %>% 
  mutate(., STATEFP =str_pad(string=State.ANSI, width=2, side='left', pad = '0'),
         COUNTYFP = str_pad(string=County.ANSI, width=3, side='left', pad = '0'),
         GEOID = paste0(STATEFP, COUNTYFP))
cattle.pad$Value <- parse_number(cattle.pad$Value)

cattle.sf <- cty %>% 
  left_join(., cattle.pad, by="GEOID") %>% 
  st_transform(., crs(housing.rescale)) 



#need to interp values
cattle.na <- cattle.sf %>% 
  filter(., is.na(Value)) %>% 
  dplyr::select(., GEOID, Value)

cattle.touches <- cattle.na %>% 
  st_touches(., cattle.sf, sparse=TRUE)

names(cattle.touches) <- cattle.na$GEOID

cattle.na.fill <- lapply(cattle.na$GEOID, function(x){
  slice(cattle.sf, cattle.touches[[x]]) %>% 
    summarise(., GEOID = x, 
              meanNB = round(mean(Value, na.rm=TRUE))) %>% st_drop_geometry()
}) %>% do.call(rbind, .)

colnames(cattle.na.fill)[2] <- "Value"

na.records <- which(is.na(cattle.sf$Value))
cattle.sf[na.records, "Value"] <- cattle.na.fill$Value

cattle.sp <- cattle.sf %>% 
  as(., "Spatial")

cattle.crop <- crop(cattle.sp, housing.rescale)
cattle.rst <- terra::rasterize(cattle.crop, housing.rescale, field="Value")
log.cattle <- log(cattle.rst, 10)
cattle.rescale <- (log.cattle - cellStats(log.cattle, min))/(cellStats(log.cattle, max) - cellStats(log.cattle, min))
writeRaster(cattle.rescale, here::here('data/ProcessedData/rasters/cattle_num.tif'), overwrite = TRUE)


# institutional capacity --------------------------------------------------
fedlands <- st_read("Data/ProcessedData/shapefiles/fedlands.shp")

unique(fedlands$Mang_Name) #Use to generate manager names for resistance reclass

# reclassify based on connectivity mandates

fedland.reclass <- data.frame(Mang_Name = unique(fedlands$Mang_Name),
                              resist = c(0.3, 0.2, 0.3, 0.7, 0.5, 1, 0.8, 1, 1,0.5,0.5))
fedlands.resist <- fedlands %>% 
  left_join(. , fedland.reclass) %>% st_cast(.,"MULTIPOLYGON") %>% 
  st_transform(., crs(housing.rescale)) %>% 
  as(., "Spatial")

fedlands.crop <- crop(fedlands.resist, housing.rescale)
fedlands.resist.rst <- terra::rasterize(fedlands.crop, housing.rescale, field="resist")
fedlands.resist.rst[is.na(fedlands.resist.rst)] <- 0.5
writeRaster(fedlands.resist.rst, here::here('data/ProcessedData/rasters/federal_instcap_rast.tif'), overwrite=TRUE)


# Institutional Capacity Scenario -----------------------------------------
fedland.reclass <- data.frame(Mang_Name = unique(fedlands$Mang_Name),
                              resist = c(0.3, 0.2, 0.3, 0.7, 0.5, 1, 0.2, 1, 1,0.5,0.5))
fedlands.resist <- fedlands %>% 
  left_join(. , fedland.reclass) %>% st_cast(.,"MULTIPOLYGON") %>% 
  st_transform(., crs(housing.rescale)) %>% 
  as(., "Spatial")

fedlands.crop <- crop(fedlands.resist, housing.rescale)
fedlands.resist.rst <- terra::rasterize(fedlands.crop, housing.rescale, field="resist")
fedlands.resist.rst[is.na(fedlands.resist.rst)] <- 0.5
writeRaster(fedlands.resist.rst, here::here('data/ProcessedData/rasters/federal_instcap_rast_scnenario.tif'), overwrite=TRUE)


