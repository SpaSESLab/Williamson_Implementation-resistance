library(terra)

house.val <- rast(here::here("data/ProcessedData/landval_crop_resamp.tif"))  #already 0-1
wolves <- rast(here::here('data/ProcessedData/wolf_crop_resamp.tif'))   #0-1 but not scaled
cattle <- rast(here::here('data/ProcessedData/cattle_num.tif'))  #0-1 rescale of log(cattle num)
jurisdiction <- rast(here::here('data/ProcessedData/federal_instcap_rast.tif'))

wolves.rscl <- (raster::raster(wolves) - raster::cellStats(raster::raster(wolves), min))/(raster::cellStats(raster::raster(wolves), max)-raster::cellStats(raster::raster(wolves), min))

house.val.1m <- 1- house.val
wolves.1m <- 1 - rast(wolves.rscl)
cattle.1m <- 1 - cattle
jurisdiction.1m <- 1 - jurisdiction

imp.fuz.sum1 <- 1 - ((wolves.1m) * (cattle.1m) * (jurisdiction.1m))
#imp.fuz.sum2 <- 1 - ((wolves.1m) * (cattle.1m) * (jurisdiction.1m) * (house.val.1m))
#imp.fuz.sum[is.na(imp.fuz.sum)] <- 5 #give NAs (lots of water) high resistance vals
writeRaster(imp.fuz.sum1, here::here('data/ProcessedData/imp_fuz_sum_nohouse.tif'), overwrite=TRUE)
#writeRaster(imp.fuz.sum2, here::here('Data/ProcessedData/imp_fuz_sum.tif'), overwrite=TRUE)

