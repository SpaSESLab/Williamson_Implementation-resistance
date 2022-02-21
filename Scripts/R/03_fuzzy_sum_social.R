library(terra)

house.val <- rast(here::here("data/ProcessedData/rasters/landval_crop_resamp.tif"))  #already 0-1
wolves <- rast(here::here('data/ProcessedData/rasters/wolf_crop_resamp.tif'))   #0-1 but not scaled
cattle <- rast(here::here('data/ProcessedData/rasters/cattle_num.tif'))  #0-1 rescale of log(cattle num)
jurisdiction <- rast(here::here('data/ProcessedData/rasters/federal_instcap_rast.tif'))
jurisdiction2 <- rast(here::here('data/ProcessedData/rasters/federal_instcap_rast_scnenario.tif'))
wolves.rscl <- (raster::raster(wolves) - raster::cellStats(raster::raster(wolves), min))/(raster::cellStats(raster::raster(wolves), max)-raster::cellStats(raster::raster(wolves), min))

house.val.1m <- 1- house.val
wolves.1m <- 1 - rast(wolves.rscl)
cattle.1m <- 1 - cattle
jurisdiction.1m <- 1 - jurisdiction
jurisdiction2.1m <- 1 - jurisdiction2

imp.fuz.sum1 <- 1 - ((wolves.1m) * (cattle.1m) * (jurisdiction.1m))
imp.fuz.sum2 <- 1 - ((wolves.1m) * (cattle.1m) * (jurisdiction2.1m))

#imp.fuz.sum[is.na(imp.fuz.sum)] <- 5 #give NAs (lots of water) high resistance vals
writeRaster(imp.fuz.sum1, here::here('data/ProcessedData/rasters/imp_fuz_sum_nohouse.tif'), overwrite=TRUE)
writeRaster(imp.fuz.sum2, here::here('data/ProcessedData/rasters/imp_fuz_sum_jurisdiction_scenario.tif'), overwrite=TRUE)

