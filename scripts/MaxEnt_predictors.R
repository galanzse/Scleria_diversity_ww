

library(tidyverse)
library(terra)
# library(rnaturalearthdata)
# library(ecospat)


# import occurrence data (already curated)
load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")


# import predictors (1km at equator)
worldclim <- rast(list.files('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio', full.names=T))
names(worldclim) <- substr(names(worldclim), 11, 16)

# explore correlation for a random subset of points
rdm_points <- data_final[['occurrences']][sample(1:nrow(data_final[['occurrences']]), 250),] %>% vect(geom=c('x','y'), 'epsg:4326')
rdm_points <- worldclim %>% terra::extract(rdm_points, ID=F) %>% as.data.frame()

# retain bio_1,5,7,12,13,15,18
cor(na.omit(rdm_points[,1:19]), method='pearson') # %>% View()

worldclim <- worldclim[[c('bio_1','bio_5','bio_7','bio_12','bio_13','bio_15','bio_18')]]


HPD <- rast('C:/Users/user/Desktop/Global_2020_PopulationDensity30sec_GPWv4.tiff') %>%
  resample(worldclim$bio_1, method='bilinear')
names(HPD) <- 'HPD'


elevation <- rast('C:/Users/user/Desktop/worldclim/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif') %>%
  resample(worldclim$bio_1, method='bilinear')
names(elevation) <- 'elevation'


landcover <- rast('SDM/global_vegetation/gm_lc_v3.tif') %>%
  resample(worldclim$bio_1, method='near')
names(landcover) <- 'landcover'
landcover[landcover%in%c(19,20)] <- NA # ice, snow and water bodies


treecover <- rast('SDM/global_vegetation/gm_ve_v2.tif')
names(treecover) <- 'treecover'
treecover[treecover>100] <- NA # water bodies and no data
treecover <- treecover %>% resample(worldclim$bio_1, method='bilinear')


ecoregions <- vect('C:/Users/user/Desktop/Ecoregions2017/Ecoregions2017.shp')['ECO_ID']
names(ecoregions) <- 'ecoregions'
ecoregions <- terra::rasterize(x=ecoregions, y=worldclim$bio_1, field='ecoregions', fun='min')


biomes <- vect('C:/Users/user/Desktop/Ecoregions2017/Ecoregions2017.shp')['BIOME_NUM']
names(biomes) <- 'biomes'
biomes <- terra::rasterize(x=biomes, y=worldclim$bio_1, field='biomes', fun='min')


# combine
predictors <- c(worldclim, HPD, elevation, ecoregions, biomes, landcover, treecover)
names(predictors)


# filter and save
writeRaster(predictors, 'SDM/rst_predictors_MaxEnt.tiff', overwrite=T)



# aggregate to similar resolution to reduce computation time

# mean
agg_bio1 <- predictors$bio_1 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_bio_5 <- predictors$bio_5 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_bio_7 <- predictors$bio_7 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_bio_12 <- predictors$bio_12 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_bio_13 <- predictors$bio_13 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_bio_15 <- predictors$bio_15 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_bio_18 <- predictors$bio_18 %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_HPD <- predictors$HPD %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_elevation <- predictors$elevation %>% aggregate(fact=5, fun='mean', na.rm=T)
agg_treecover <- predictors$treecover %>% aggregate(fact=5, fun='mean', na.rm=T)

# modal
agg_ecoregions <- predictors$ecoregions %>% aggregate(fact=5, fun='modal', na.rm=T)
agg_biomes <- predictors$biomes %>% aggregate(fact=5, fun='modal', na.rm=T)
agg_landcover <- predictors$landcover %>% aggregate(fact=5, fun='modal', na.rm=T)


# combine
agg_predictors <- c(agg_bio1, agg_bio_5, agg_bio_7, agg_bio_12, agg_bio_13, agg_bio_15, agg_bio_18, agg_HPD, agg_elevation, agg_ecoregions, agg_biomes, agg_landcover, agg_treecover)  

rm(agg_bio1, agg_bio_5, agg_bio_7, agg_bio_12, agg_bio_13, agg_bio_15, agg_bio_18, agg_HPD, agg_elevation, agg_ecoregions, agg_biomes, agg_landcover, agg_treecover)


# filter and save
writeRaster(agg_predictors, 'SDM/agg_predictors_MaxEnt.tiff', overwrite=TRUE)


