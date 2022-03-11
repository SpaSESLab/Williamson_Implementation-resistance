library(raster)
library(rasterVis)
library(ggplot2)
library(tidyverse)
library(cmocean)

# load inital cs layers ---------------------------------------------------
implement.cs <- raster(here::here('data/ProcessedData/implement/implement_cum_curmap.asc'))
biophys.cs <- raster(here::here('data/ProcessedData/biophys/biophys_cum_curmap.asc'))
jurisdiction.cs <- raster(here::here('data/ProcessedData/jurisdiction/jurisdiction_cum_curmap.asc'))
cattle.cs <- raster(here::here("data/ProcessedData/cows/cows_cum_curmap.asc"))


# convert to dataframes ---------------------------------------------------

biophys.df <- biophys.cs %>%
  projectRaster(., res=300, crs = crs(biophys.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "biophys"))

implement.df <- implement.cs %>%
  projectRaster(., res=300, crs = crs(biophys.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "implement")) %>% 
  mutate(imp.norm = (implement - min(implement, na.rm=TRUE))/(max(implement, na.rm=TRUE)-min(implement, na.rm=TRUE)))

jurisdiction.df <- jurisdiction.cs %>% 
  projectRaster(., res=300, crs = crs(biophys.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "jurisdiction")) %>% 
  mutate(juris.norm = (jurisdiction - min(jurisdiction, na.rm=TRUE))/(max(jurisdiction, na.rm=TRUE)-min(jurisdiction, na.rm=TRUE)))
  
cattle.df <- cattle.cs %>% 
  projectRaster(., res=300, crs = crs(biophys.cs)) %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "cattle"))  %>% 
  mutate(cattle.norm = (cattle - min(cattle, na.rm=TRUE))/(max(cattle, na.rm=TRUE)-min(cattle, na.rm=TRUE)))


# Generate quantiles for each layer ---------------------------------------
biophys.quant <- c(quantile(biophys.df$biophys, probs = seq(0, 1, 0.2), type = 5), Inf)
implement.quant <- c(quantile(implement.df$implement, probs = seq(0, 1, 0.2), type = 5), Inf)
implement.norm.quant <- c(quantile(implement.df$imp.norm, probs = seq(0, 1, 0.2), type = 5), Inf)

# reclass -----------------------------------------------------------------

biophys.reclass <- biophys.df %>% 
  dplyr::mutate(., PCT = as.integer(cut(biophys, breaks = biophys.quant, include.lowest = TRUE, ordered_result=TRUE)))

implement.reclass <- implement.df %>%
  dplyr::mutate(., PCT = as.integer(cut(implement, breaks = implement.quant, include.lowest = TRUE, ordered_result = TRUE)),
                PCTnorm = as.integer(cut(imp.norm, breaks = implement.norm.quant, include.lowest = TRUE, ordered_result = TRUE)))

juris.reclass <- jurisdiction.df %>% 
  dplyr::mutate(., PCT = as.integer(cut(jurisdiction, breaks = implement.quant, include.lowest = TRUE, ordered_result = TRUE)),
                PCTnorm = as.integer(cut(juris.norm, breaks = implement.norm.quant, include.lowest = TRUE, ordered_result = TRUE)),
                diff = PCT - implement.reclass$PCT,
                difnorm = PCTnorm - implement.reclass$PCTnorm)


cattle.reclass <- cattle.df %>% 
  dplyr::mutate(., PCT = as.integer(cut(cattle, breaks = implement.quant, include.lowest = TRUE, ordered_result = TRUE)),
                PCTnorm = as.integer(cut(cattle.norm, breaks = implement.norm.quant, include.lowest = TRUE, ordered_result = TRUE)),
                diff = PCT - implement.reclass$PCT,
                difnorm = PCTnorm - implement.reclass$PCTnorm)


df.join.base <- left_join(biophys.reclass, implement.reclass, by = c("x","y")) %>%
  mutate(
    group = paste(PCT.x, PCT.y, sep = " - ")
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


col.matrix<-colmat(nquantiles=5)


# set up legend -----------------------------------------------------------
groups <- unique(df.join.base$group)


legend.5<- data.frame(matrix(ncol=2, nrow=25))
colnames(legend.5) <- c("group", "fill")  
legend.5$group <- groups[order(groups)]  
legend.5$fill <- as.vector(col.matrix[2:6,2:6])


df.join.base <- left_join(df.join.base, legend.5)

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
elev.proj <- projectRaster(elev[[1]], biophys.cs)


elev.crop <- crop(elev.proj, biophys.cs)
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
origin.proj <- spTransform(origins, crs(biophys.cs))
goals.proj <- spTransform(goals, crs(biophys.cs))
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


# Make the plots ----------------------------------------------------------


p1 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = df.join.base,
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
  theme(legend.position = 'none')


l <- legend.5 %>%
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
    plot.background = element_rect(fill="white", color = "white")) +
  coord_fixed()



p2 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = juris.reclass,
    aes(
      x=x,
      y=y,
      fill=difnorm
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_cmocean(name="curl") +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "black")+
  guides(fill = guide_colorbar(nbin = 8,  title.position = "top", 
                               label.position="bottom", title = "\u0394 in rank")) +
  theme_map() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = "center"  )

p3 <- RStoolbox::ggR(hills3) +
  geom_raster(
    data = cattle.reclass,
    aes(
      x=x,
      y=y,
      fill=difnorm
    ),
    alpha=0.8,
    interpolate = TRUE
  ) +
  scale_fill_cmocean(name="curl", limits = c(-3,3)) +
  geom_sf(data = PAs, fill = "forestgreen") +
  geom_sf(data = as(sts.crop, "sf"), fill = NA, color="black")+
  ggrepel::geom_text_repel(data = pa.cents, aes(x = st_coordinates(pa.cents)[,1], y = st_coordinates(pa.cents)[,2], label=lab), nudge_x = -40000 , nudge_y = c(80000,90000), fontface="bold", color = "black")+
  guides(fill = guide_colorbar(nbin = 8,  title.position = "top", 
                               label.position="bottom", title = "\u0394 in rank")) +
  theme_map() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = "center"  )

deltas <- p2 + p3 + plot_layout(guides = 'collect') & theme(legend.position = "bottom", legend.justification = "center",
                                                                                                                legend.margin=margin(0,0,0,0),
                                                                                                                legend.box.margin=margin(-15,-10,-10,-10))

triplot <- p1 + p2 + p3 + plot_annotation(tag_levels = "A", tag_suffix = ")")& theme(legend.position = "bottom", legend.justification = "center",
                                                                    legend.margin=margin(0,0,0,0),
                                                                    legend.box.margin=margin(-15,-10,-10,-10))   


p <- ggdraw(triplot)  +
  draw_plot(l, x = 0, y = 0, 
            width = 0.24, height = 0.24) 

ggsave(here::here("plots/fig3.png"), plot =p)
