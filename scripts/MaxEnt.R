

# SDM to estimate potential localities of species
# 1/ Determine minimum sample size (at the moment N=15)
# 2/ Choose predictors considering we are working with herbarium data (i.e., biased) (at the moment elevation, worldclim and percentage tree cover)
# 3/ Set area to be predicted (at the moment ecoregions within the EOO)
# 4/ Selection of background points - sampling bias (target-group approach)
# 4/ Fit model and retain best parametization (ENMevaluate)
# 5/ Convert ROR into presence absence


library(tidyverse)
library(terra)
library(dismo)
library(ENMeval)
library(ks)
library(spatialEco)
library(rnaturalearth)


# import occurrence data (already curated)
load("data/data_final.RData")

# occurrences
scl_occurrences <- data_final[['occurrences']]

# remove isolated point in S. polycarpa
out <- which(scl_occurrences$scientific_name=='Scleria polycarpa Boeckeler' & scl_occurrences$x < c(-150))
scl_occurrences[out,] <- NA
scl_occurrences <- na.omit(scl_occurrences)


# import ecoregions maps (version Dinerstein et al 2017)
ecoregions <- terra::vect('data/Ecoregions2017/Ecoregions2017.shp')

# define AOI: ecoregions with observations
AOI <- ecoregions %>% terra::intersect(vect(data_final[['occurrences']], geom=c('x','y'), 'epsg:4326'))
AOI <- AOI$ECO_NAME %>% unique()
AOI <- terra::subset(ecoregions, ecoregions$ECO_NAME %in% AOI)

# coastlines and countries for plots
coastlines <- ne_coastline(scale = 110, returnclass = "sf") %>% vect()
countries <- ne_countries(scale = 110, type = "countries", returnclass = "sf") %>% vect()

# import predictors and crop
agg_predictors <- terra::rast('SDM/agg_predictors_MaxEnt.tiff')[[c("bio_1","bio_5","bio_7","bio_12","bio_13", "bio_15", "bio_18", "elevation", "treecover")]] %>% terra::crop(AOI, mask=T)
names(agg_predictors)
# cellSize(agg_predictors)


# # vector of species with at least 15 observations
# n_obs <- table(data_final[['occurrences']]$scientific_name) %>% as.data.frame() %>%
#   subset(Freq>=15)
# sdm_species <- n_obs[order(n_obs$Freq, decreasing=F),]
# colnames(sdm_species)[1] <- 'scientific_name'
# 
# 
# # loop to model distributions and obtain sites with high habitat suitability for posterior analyses
# suitable_locations <- list()
# allmodels <- list()
# summary_bestmodels <- sdm_species
# summary_bestmodels$parameters <- NA
# summary_bestmodels$auc.diff.avg <- NA
# summary_bestmodels$auc.val.avg <- NA
# summary_bestmodels$or.mtp.avg <- NA

load("SDM/suitable_locations.Rdata")
load("SDM/allmodels.Rdata")
summary_bestmodels <- read.csv("SDM/summary_bestmodels.txt", sep="")

for (s1 in 1:nrow(summary_bestmodels)) { # 1:nrow(summary_bestmodels)
  
  spp1 = summary_bestmodels$scientific_name[s1]

  # occurrences
  occ_temp <- scl_occurrences %>% filter(scientific_name==spp1)

  # points
  pts_temp <- vect(occ_temp, geom=c('x','y'), 'epsg:4326')

  # model IN and NEAR ecoregions within its EOO
  ext_temp <- terra::convHull(pts_temp)
  eco_temp <- AOI %>% terra::mask(ext_temp) # select geometries of x that intersect with the geometries of y

  # crop predictors
  temp_predictors <- agg_predictors %>% terra::crop(eco_temp, mask=T)
  

  # plot(temp_predictors$bio_1)
  # lines(ext_temp, col='red')
  # lines(countries)
  # points(pts_temp, col='blue')
  
  
  # Bias file: Target-group background (https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13442)
  target_points <- data_final[['occurrences']][,c('x','y')] %>% unique() %>% 
    terra::vect(geom=c('x','y'), 'epsg:4326') %>%
    terra::intersect(eco_temp) # all Scleria observations for our study area
  
  # 2d kernel density estimation 
  target_density <- ks::kde(geom(target_points)[,c('x','y')]) 
  
  # resample + crop to ground
  target_raster <- raster::raster(target_density) %>% terra::rast() %>%
    terra::resample(temp_predictors$bio_1, method='bilinear') %>%
    terra::crop(eco_temp, mask=T)
  
  # normalize
  target_raster <- target_raster - min(target_raster[], na.rm=T)
  target_raster <- raster.transformation(target_raster, trans="norm") 

  # select a maximum of 10000 background points
  bg_temp <- as.data.frame(target_raster, xy=T) %>% subset(layer>0.1)
  if (nrow(bg_temp) > 10000) { bg_temp <- bg_temp[sample(1:nrow(bg_temp), prob=bg_temp$layer, size=10000, replace=F),] }

  
  # plot(target_raster)
  # points(target_points, col='red')
  # points(bg_temp, col='orange')
  # lines(countries)
  # points(pts_temp, col='blue')
  
  
  # convert to Raster* file to run ENMevaluate
  stk_predictors <- raster::stack(temp_predictors) # Raster* file
  
  # Get best model: ENMeval
  modeval <- ENMevaluate(occs = occ_temp[,c('x','y')], 
                         envs = stk_predictors,
                         # categoricals = c('ecoregions','biomes','landcover'),
                         bg= bg_temp[,c('x','y')],
                         algorithm = 'maxnet',
                         RMvalues = 1:5,
                         tune.args = list(fc = c('L', 'Q', 'LQ', 'LQH')), 
                         partitions = "randomkfold", partition.settings = list(kfolds = 10), 
                         clamp = TRUE, 
                         parallel = TRUE, numCores = NULL) # 8
  
  
  # save results and select best model
  results1 <- modeval@results
  results1$nfeat <- nchar(as.character(results1$fc)) # number of features
  allmodels[[spp1]] <- results1

  # best model (lowest AICc)
  bestmod1 <- results1 %>% subset(delta.AICc<2)
  if (nrow(bestmod1)>1) { bestmod1 <- bestmod1 %>% subset(ncoef==min(bestmod1$ncoef, na.rm=T)) }
  if (nrow(bestmod1)>1) { bestmod1 <- bestmod1 %>% subset(nfeat==min(bestmod1$nfeat, na.rm=T)) }
  if (nrow(bestmod1)>1) { bestmod1 <- bestmod1[1,] }
  
  summary_bestmodels$parameters[s1] <- as.character(bestmod1$tune.args)
  summary_bestmodels$auc.diff.avg[s1] <- bestmod1$auc.diff.avg
  summary_bestmodels$auc.val.avg[s1] <- bestmod1$auc.val.avg
  summary_bestmodels$or.mtp.avg[s1] <- bestmod1$or.mtp.avg


  # obtain best prediction
  map1 <- terra::rast(modeval@predictions[[as.numeric(rownames(bestmod1))]])

  
  # plot(map1, main="Model with lowest AICc")
  # points(pts_temp, col='blue')
  
  
  # use 10th percentile as threshold
  q10 <- map1 %>% terra::extract(pts_temp, ID=F) %>% deframe() %>% quantile(probs=0.1, na.rm=T)
  map1[map1>=q10] <- 1; map1[map1<q10] <- 0
  
  
  # save map
  mypath <- file.path(paste(getwd(),'/SDM/maps/',spp1,'.tiff',sep=''))
  tiff(file=mypath, width=600, height=600, units='px')
  plot(map1, main=spp1, legend=F)
  points(pts_temp, col='black')
  lines(countries)
  dev.off()
  
  
  # retain suitable locations
  map1[map1!=1] <- NA
  sl1 <- as.data.frame(map1, xy=T)[,c('x','y')]
  sl1$scientific_name <- spp1
  suitable_locations[[s1]] <- sl1

  # save results
  save(suitable_locations, file='SDM/suitable_locations.Rdata')
  save(allmodels, file='SDM/allmodels.Rdata')
  write.table(summary_bestmodels, 'SDM/summary_bestmodels.txt')
  
  View(summary_bestmodels)
  print(spp1)

}



# transform habitat suitability into presence based on reference grid ####
# let's use 150*150km cells: contains 450 5*5km cells (approx.), thus we retains cells with min 450 obs + observations
world_grid <- rast('results/world_grid.tiff')
sqrt(cellSize(world_grid))

filtered_suitable <- list()

for (s1 in 106:length(suitable_locations)) {

  sp1 = unique(suitable_locations[[s1]]$scientific_name)

  # observations: reproject and save as df
  df_obs <- scl_occurrences[scl_occurrences$scientific_name==sp1,] %>%
    vect(geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth') %>%
    geom() %>% as.data.frame() %>% dplyr::select(x,y)
  df_obs$scientific_name <- sp1
  
  # expectation > reproject to equal area
  df_exp <- suitable_locations[[s1]] %>%
    vect(geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth')

  # filter expected values based on 50% habitat suitability
  df_exp <- world_grid %>% terra::extract(df_exp, xy=T, cells=T, ID=F)
  cll1 <- table(df_exp$cell)[table(df_exp$cell)>450] %>% names() %>% as.numeric()

  # filter suitable cells and merge with observations
  if (length(cll1)>0) {
    df_exp <- df_exp %>% subset(cell%in%cll1) %>% unique()
    df_exp$scientific_name <- sp1
    df_exp$reference_grid <- NULL; df_exp$cell <- NULL
    df_obs <- rbind(df_obs, df_exp[,colnames(df_obs)])
  }
  
  # remove observations from ecoregions that do not have Scleria
  pts_obs <- df_obs %>% vect(geom=c('x','y'), '+proj=eqearth') %>% project('epsg:4326')
  obs_in <- terra::extract(x=AOI, y=pts_obs)
  df_obs <- df_obs[which(!is.na(obs_in$ECO_NAME)),]
  
  # save
  filtered_suitable[[s1]] <- df_obs
  
  print(s1)
}

# collapse
df_filtered_suitable <- do.call("rbind", filtered_suitable)

# add original observations from other species
temp <- scl_occurrences %>% subset(!(scientific_name%in%df_filtered_suitable$scientific_name))
temp2 <- temp %>% vect(geom=c('x','y'), 'epsg:4326') %>% project('+proj=eqearth') %>%
  geom(df=T) %>% dplyr::select(x,y)
temp2$scientific_name <- temp$scientific_name

# bind
df_filtered_suitable <- rbind(df_filtered_suitable, temp2[,colnames(df_filtered_suitable)])

# save
save(filtered_suitable, file="results/filtered_suitable.RData")
write.table(df_filtered_suitable, 'results/df_filtered_suitable.txt')


