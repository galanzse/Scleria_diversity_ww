

library(readxl)
library(tidyverse)
library(terra)
library(DescTools)



# import data
load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")


# vect
ecoregions <- vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp')

# create unique ID for each combination REALM x BIOME
ecoregions$REALMxBIOME <- paste(ecoregions$REALM, ecoregions$BIOME, sep='')

# rast
rast_ecoregions <- rasterize(x=project(ecoregions,'+proj=eqearth'), y=world_grid, field='REALMxBIOME', touches=TRUE)


# extract most frequent biome nested within realm
df_ecoregions <- data_final[['occurrences']][,c('scientific_name','x','y')]

scl_points2 <- vect(df_ecoregions, geom=c('x','y'),crs="epsg:4326") %>% project('+proj=eqearth')

df_ecoregions$REALMxBIOME <- rast_ecoregions %>% terra::extract(y=scl_points2) %>%
  select(REALMxBIOME) %>% deframe()

df_ecoregions <- unique(df_ecoregions[,c('scientific_name','REALMxBIOME')]) %>% na.omit()

rm(scl_points2)


# # correct levels
# df_ecoregions$REALM <- as.factor(df_ecoregions$REALM)
# levels(df_ecoregions$REALM) <- c("Australasia", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic")
# 
# df_ecoregions$BIOME <- as.factor(df_ecoregions$BIOME)
# levels(df_ecoregions$BIOME) <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Coniferous Forests","Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Tropical & Subtropical Grasslands, Savannas & Shrublands","Temperate Grasslands, Savannas & Shrublands","Flooded Grasslands & Savannas","Montane Grasslands & Shrublands","Deserts & Xeric Shrublands","Mangroves",NA)
# 
# df_ecoregions$GBL_STAT <- as.factor(df_ecoregions$GBL_STAT)
# levels(df_ecoregions$GBL_STAT) <- c(NA,"CRITICAL OR ENDANGERED", "VULNERABLE", "RELATIVELY STABLE OR INTACT")



# save
write.table(df_ecoregions, 'results/df_ecoregions.txt')
writeRaster(rast_ecoregions, 'results/rast_ecoregions.tiff', overwrite=TRUE)


