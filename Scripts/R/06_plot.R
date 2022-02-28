library(raster)
library(gdistance)
library(rgdal)
library(sf)
library(magrittr)
library(ggplot2)
library(rasterVis)
library(patchwork)
library(classInt)
library(cowplot)
library(ggsci)
library(stringr)

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
# load Circuitscape results -----------------------------------------------


biophys.cs <- raster(here::here('data/ProcessedData/biophys/biophys_cum_curmap.asc'))
implement.cs <- raster(here::here('data/ProcessedData/implement/implement_cum_curmap.asc'))
jurisdiction.cs <- raster(here::here('data/ProcessedData/jurisdiction/jurisdiction_cum_curmap.asc'))
cattle.cs <- raster(here::here("data/ProcessedData/cows/cows_cum_curmap.asc"))
biophys.norm <- (biophys.cs - cellStats(biophys.cs, min))/(cellStats(biophys.cs, max)-cellStats(biophys.cs, min))
implement.norm <- (implement.cs - cellStats(implement.cs, min))/(cellStats(implement.cs, max)-cellStats(implement.cs, min))
jurisdiction.norm <- (jurisdiction.cs - cellStats(jurisdiction.cs, min))/(cellStats(jurisdiction.cs, max)-cellStats(jurisdiction.cs, min))
cattle.norm <- (cattle.cs - cellStats(cattle.cs, min))/(cellStats(cattle.cs, max)-cellStats(cattle.cs, min))

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
  dplyr::summarise(., sum.current = round(sum(biophys.cur, na.rm = TRUE), 2)) %>%
  dplyr::mutate(., ID = str_remove_all( LCPID, fixed("b"))) %>% 
  dplyr::select(., c(ID, sum.current)) %>% 
  flextable::flextable() %>% 
  flextable::set_header_labels(., ID = "Cost Rank", sum.current = "Total Current") %>% 
  flextable::as_raster(zoom = 3)


#dens.plot <- ggplot(current.dens, aes(x=LCPID, y= sum.current, color=LCPID), size=8) +
#  scale_color_locuszoom() +
#  geom_point(stat = "identity") +
#  coord_flip() +
#  theme_cowplot()


# Estimate cost.ratio -----------------------------------------------------

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
  dplyr::mutate(., ratio = soc.costdist1/bio.costdist,
                ratio_jur = jurisdiction.costdist/bio.costdist,
                ratio_cattle = cattle.costdist/bio.costdist) %>% 
  dplyr::group_by(., LCPID) %>% 
  dplyr::filter(.,  startsWith(LCPID, "b")) %>% 
  dplyr::summarise(., mn_ratio = mean(ratio, na.rm = TRUE), 
                   sd_ratio = sd(ratio, na.rm=TRUE),
                   mean_delta = mean(delta_bio.costdist, na.rm=TRUE),
                   mn_ratio_jur = mean(ratio_jur, na.rm = TRUE), 
                   sd_ratio_jur = sd(ratio_jur, na.rm=TRUE),
                   mn_ratio_cattle = mean(ratio_cattle, na.rm = TRUE), 
                   sd_ratio_cattle = sd(ratio_cattle, na.rm=TRUE)) %>% 
  dplyr::mutate(., ID = str_remove_all( LCPID, fixed("b")))
  


p1 <- ggplot(data = dist.group, mapping = aes(x = ID, 
                                             y = mn_ratio, color = log(mean_delta, 10))) +
  scale_color_viridis_c(direction = -1, begin= 0, end =0.8)

ratio.plot1 <- p1 + geom_pointrange(mapping = aes(ymin = mn_ratio - sd_ratio,
                                  ymax = mn_ratio + sd_ratio)) +
  labs(x= "Cost Rank", y= "Cost ratio") + 
  guides(color = guide_colorbar(title = "\u0394 log(biophysical \ncost)")) +
  theme_cowplot()


p2 <- ggplot(data = dist.group, mapping = aes(x = ID, 
                                              y = mn_ratio_jur, color = log(mean_delta, 10))) +
  scale_color_viridis_c(direction = -1, begin= 0, end =0.8)

ratio.plot2 <- p2 + geom_pointrange(mapping = aes(ymin = mn_ratio_jur - sd_ratio_jur,
                                                  ymax = mn_ratio_jur + sd_ratio_jur)) +
  labs(x= "Cost Rank", y= "Cost ratio") + 
  guides(color = guide_colorbar(title = "\u0394 log(biophysical \ncost)")) +
  theme_cowplot()

p3 <- ggplot(data = dist.group, mapping = aes(x = ID, 
                                              y = mn_ratio_cattle, color = log(mean_delta, 10))) +
  scale_color_viridis_c(direction = -1, begin= 0, end =0.8)

ratio.plot3 <- p3 + geom_pointrange(mapping = aes(ymin = mn_ratio_cattle - sd_ratio_cattle,
                                                  ymax = mn_ratio_cattle + sd_ratio_cattle)) +
  labs(x= "Cost Rank", y= "Cost ratio") + 
  guides(color = guide_colorbar(title = "\u0394 log(biophysical \ncost)")) +
  theme_cowplot()

#possibly facet wrap this
ratio.plot1 + ratio.plot2 + ratio.plot3



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
elev.proj <- projectRaster(elev[[1]], dist.stack)


elev.crop <- crop(elev.proj, dist.stack)
elev.mod <- elev.crop *10
slope <- terrain(elev.mod, opt='slope')
aspect <- terrain(elev.mod, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)
hill2 <- aggregate(hill , fact = 5 , method = "bilinear" )
hills3 <- focal(hill2, w=matrix(1/9, nc=3, nr=3), mean)


# Get vectors for maps ----------------------------------------------------

PAs <- st_read(here::here("Data/ProcessedData/shapefiles/studyPAs.shp")) %>% 
  dplyr::filter(. , Unit_Nm == "Yellowstone National Park" | Unit_Nm == "Weminuche Wilderness") %>% 
  st_transform(. , crs = crs(dist.stack))

sts <- tigris::states() %>% 
  dplyr::filter(., STUSPS %in% c("ID", "MT", "UT", "WY", "CO", "AZ", "NM")) %>% 
  st_transform(., crs(dist.stack)) %>% 
  as(., "Spatial")

conus <-  tigris::states() %>%
  dplyr::filter(,, !STUSPS %in% c("AK","AS", "HI", "GU", "VI", "DC", "PR", "MP")) %>% 
  st_transform(., crs(dist.stack))

cs.proj <- biophys.cs %>%
  projectRaster(., res=300, crs = crs(dist.stack)) 

sts.crop <- crop(sts, cs.proj)

pa.cents <- rbind(as(origin.proj,"sf"), as(goals.proj, "sf"))
pa.cents$lab <- c("Weminche \n Wilderness Area", "Yellowstone \n National Park")
# Plot Circuitscape results -----------------------------------------------
library(purrr)
library(dplyr)


b_df <- cs.proj %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "biophys"))



biophys.lcps <- lapply(biophys.lst, function (x) rasterToPolygons(x=x, n=8, dissolve = TRUE))

biophys.lcp.sf <- list(biophys.lcps, makeUniqueIDs = T) %>% 
  flatten() %>% 
  do.call(rbind, .) %>% 
  as(., "sf") %>% 
  dplyr::mutate(., rank = row_number())



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
  scale_color_locuszoom() +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="white")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "white")+
  guides(fill = guide_colorbar(nbin = 8,  title.position = "left", 
                               label.position="right", title = "Current flow"), 
         color = guide_legend(title.position = "left", 
                              label.position="right", title = "Cost rank"))+
  theme_map() + theme(legend.position = "left", 
                      legend.box = "vertical", 
                      legend.title = element_text(angle = 90),
                      plot.margin = unit(c(5.5, 1.5, 0, 1.5), units = "pt"),
                      legend.box.margin = unit(c(0, 0.1, 0.1, 0.1), units = "pt"))

inset <- ggplot()+
  geom_sf(data=conus, fill="white") +
  geom_sf(data =st_as_sfc(st_bbox(dist.stack)), fill=NA, color="red") +
  theme_map() +
  theme(panel.background = element_rect(fill = "gray", color = "black"),
        plot.background = element_rect(fill = NA, color = NA))  


p.combined <- ggdraw(p.cs) +
  draw_plot(inset, x = 0.65, y = 0,  
            width = 0.2, height = 0.1)

panel_a <- plot_grid(p.combined,
                     nrow = 1,
                     labels = c("A)"),
                     label_size = 16,
                     label_fontfamily = "sans")

p2 <- ggplot() + 
  theme_void() +
  theme(plot.margin = unit(c(1.5, 1.5, 0, 1.5), units = "pt")) +
  annotation_custom(grid::rasterGrob(current.dens), xmin=-0.07, xmax=0.5, ymin=-Inf, ymax=Inf)


panel_b <- plot_grid(ratio.plot, p2,
                               # A negative rel_height shrinks space between elements
                               rel_heights = c( 1,1),
                               align = "v",
                               axis = "l",
                               ncol = 1,
                               labels = c("B)","C)"),
                               label_size = 16,
                               label_fontfamily = "sans")


panel_ABC <- plot_grid(panel_a, panel_b,
                     rel_widths = c(3,2),
                     nrow = 1)

cowplot::save_plot(here::here("plots/fig1.png"), plot = panel_ABC, base_height = 7.5, base_width = 10)

combined <- p.combined | ratio.plot + patchwork::plot_layout(widths = c(2, 1))
ggsave(here::here("plots/fig1.png"),  plot = combined, height = 7, width = 10, units = "in", dpi = "screen")

# Plot bivariate raster ---------------------------------------------------
library(tidyr) #don't load this initially as it masks all versions of extract
# convert gridded raster dato dataframe
b_df <- biophys.cs %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "biophys")) %>% 
  mutate(., biophys.scl = ((biophys - min(biophys, na.rm=TRUE))/(max(biophys, na.rm = TRUE) - min(biophys, na.rm = TRUE))) )

s_df <- implement.cs %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "implement")) %>% 
  mutate(., implement.scl = ((implement - min(implement, na.rm=TRUE))/(max(implement, na.rm = TRUE) - min(implement, na.rm = TRUE))) )


j_df <- jurisdiction.cs %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "jurisdiction")) %>% 
  mutate(., juris.scl = ((jurisdiction - min(jurisdiction, na.rm=TRUE))/(max(jurisdiction, na.rm = TRUE) - min(jurisdiction, na.rm = TRUE))) )


c_df <- cattle.norm %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "cattle"))%>% 
  mutate(., cattle.scl = ((cattle - min(cattle, na.rm=TRUE))/(max(cattle, na.rm = TRUE) - min(cattle, na.rm = TRUE))) )


b_reclass <- b_df %>% 
  mutate(PCT = ntile(biophys.scl, 5))  %>%  data.frame()

s_reclass <- s_df %>% 
  mutate(PCT = ntile(implement.scl, 5))  %>%  data.frame()

rc <- matrix(c(quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[1:2],1, 
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[2:3], 2,
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[3:4], 3,
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[4:5],4,
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[5], Inf, 5), 
             ncol = 3, byrow = TRUE)
implement_quant <- quantile(s_reclass$implement.scl, probs = seq(0, 1, 0.2))[1:5]
j_reclass <- j_df %>% 
  mutate(PCT = as.integer(cut(juris.scl, breaks = c(implement_quant, Inf), right = FALSE)))  %>%  data.frame()

c_reclass <- c_df %>% 
  mutate(PCT = as.integer(cut(cattle.scl, breaks = c(implement_quant, Inf), right=FALSE)))  %>%  data.frame()


df_base <- left_join(b_reclass, s_reclass, by = c("x","y")) %>%
  mutate(
    group = paste(PCT.y, PCT.x, sep = " - ")
  ) %>%
  dplyr::select(-c(PCT.x, PCT.y))

df_juris <- left_join(b_reclass, j_reclass, by = c("x","y")) %>%
  mutate(
    group = paste(PCT.y, PCT.x, sep = " - ")
  ) %>%
  dplyr::select(-c(PCT.x, PCT.y))

df_cattle <-left_join(b_reclass, c_reclass, by = c("x","y")) %>%
  mutate(
    group = paste(PCT.y, PCT.x, sep = " - ")
  ) %>%
  dplyr::select(-c(PCT.x, PCT.y))

# Build a palette ---------------------------------------------------------
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}


col.matrix<-colmat(nquantiles=9)




legend_5 <- tibble(
  "5 - 5" = "#820050",
  "4 - 5" = "#FF05BF",
  "3 - 5" = "#FF0580", 
  "2 - 5" = "#FF0540",
  "1 - 5" = "#FF0500",
  "5 - 4" = "#BF05FF",
  "4 - 4" = "#BF05BF",
  "3 - 4" = "#BF0580", 
  "2 - 4" = "#BF0540",
  "1 - 4" = "#BF0500",
  "5 - 3" = "#8005FF",
  "4 - 3" = "#8005BF",
  "3 - 3" = "#800580", 
  "2 - 3" = "#800540",
  "1 - 3" = "#800500",
  "5 - 2" = "#4005FF",
  "4 - 2" = "#4005BF",
  "3 - 2" = "#400580", 
  "2 - 2" = "#400540",
  "1 - 2" = "#400500",
  "5 - 1" = "#0005FF",
  "4 - 1" = "#0005BF",
  "3 - 1" = "#000580", 
  "2 - 1" = "#000540",
  "1 - 1" = "#000500"
) %>%
  tidyr::gather("group", "fill")

legend_5$fill[1:5] <- col.matrix[10,seq(10, 0, -2)]
legend_5$fill[6:10] <- col.matrix[8,seq(10, 0, -2)]
legend_5$fill[11:15] <- col.matrix[6,seq(10, 0, -2)]
legend_5$fill[16:20] <- col.matrix[4,seq(10, 0, -2)]
legend_5$fill[21:25] <- col.matrix[2,seq(10, 0, -2)]
df_base <- left_join(df_base, legend_5)
df_juris <- left_join(df_juris, legend_5)
df_cattle <- left_join(df_cattle, legend_5)

p1 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = df_base,
    aes(
      x=x,
      y=y,
      fill=fill
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_identity() +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "white")+
  theme_map() +
  theme(legend.position = 'none',
        plot.margin = unit(c(5.5, 1.5, 0, 1.5), units = "pt"))


p2 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = df_juris,
    aes(
      x=x,
      y=y,
      fill=fill
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_identity() +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  theme_map() +
  theme(legend.position = 'none',
        plot.margin = unit(c(5.5, 1.5, 0, 1.5), units = "pt"))

p3 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = df_cattle,
    aes(
      x=x,
      y=y,
      fill=fill
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_identity() +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  theme_map() +
  theme(legend.position = 'none',
        plot.margin = unit(c(5.5, 1.5, 0, 1.5), units = "pt"))


p1+p2+p3

l <- legend_5 %>%
  separate(group,
           into = c("biophys", "implement"),
           sep = " - ") %>%
  mutate(biophys = as.integer(biophys),
         implement = as.integer(implement)) %>%
  ggplot() +
  geom_tile(mapping = aes(
    x = biophys,
    y = implement,
    fill = fill)) +
  scale_fill_identity() +
  labs(x = "Biophysical →\n probability",
       y = "Implementation →\n probability") +
  theme_void() +
  theme(
    axis.title = element_text(
      size = 8,
    ),
    axis.title.y = element_text(angle = 90),
    plot.background = element_rect(fill="white", color = "black")) +
  coord_fixed()

p.combined <- ggdraw(p1) +
  draw_plot(inset, x = 0.48, y = 0,  
            width = 0.2, height = 0.1)

# create final layout
p <- ggdraw(p.combined)  +
  draw_plot(p2, x = 0.53, y = 0.74, 
            width = 0.26, height = 0.26) 


cowplot::save_plot(here::here("plots/fig2.png"),  plot = p, base_height = 7.5, units = "in", dpi = "screen")

# need to think about a delta, maybe?
rc <- matrix(c(quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[1:2],1, 
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[2:3], 2,
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[3:4], 3,
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[4:5],4,
               quantile(combined.baseline, probs = seq(0, 1, 0.2), names = FALSE)[5], Inf, 5), 
             ncol = 3, byrow = TRUE)

combined.baseline <- biophys.norm * implement.norm
combined.reclass <- reclassify(combined.baseline, rc)
combined.jurisdiction <- reclassify(biophys.norm * jurisdiction.norm, rc)
combined.cattle <- reclassify(biophys.norm * cattle.norm, rc)

delta.juris <- combined.jurisdiction - combined.reclass
delta.cattle <- combined.cattle - combined.reclass

juris.delta.df <- delta.juris %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "delta")) %>% data.frame()

cattle.delta.df <- delta.cattle %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "delta")) %>% data.frame()


tst.1 <- ggplot() +
  geom_raster(
    data = juris.delta.df,
    aes(
      x=x,
      y=y,
      fill=delta
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_scico(palette = "cork", midpoint = 0) +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  theme_map() +
  theme(legend.position = 'bottom',
        plot.margin = unit(c(5.5, 1.5, 0, 1.5), units = "pt"))


tst.2 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = cattle.delta.df,
    aes(
      x=x,
      y=y,
      fill=delta
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_scico(palette = "vik") +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  theme_map() +
  theme(legend.position = 'bottom',
        plot.margin = unit(c(5.5, 1.5, 0, 1.5), units = "pt"))


juris.df <- combined.jurisdiction %>%
  projectRaster(., res=300, crs = crs(dist.stack)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "combined")) %>% 
  mutate(PCT = ntile(combined, 5))  %>%  data.frame()
