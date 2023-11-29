

# SDM to estimate potential localities of species
# STEPS:
# 1/ Determine minimum sample size (at the moment N=10)
# 2/ Choose predictors considering we are working with herbarium data (i.e., biased) (at the moment elevation and worldclim)
# 3/ Set area to be predicted (at the moment ecoregions within the EOO)
# 4/ Selection of background points - sampling bias (at the moment 100km around each observation)
# 4/ Fit model and retain best parametization (ENMevaluate)
# 5/ Convert ROR into presence absence (at the moment mean predicted value at a presence location)


library(tidyverse)
library(terra)
library(dismo)
library(ENMeval)
library(ks)
library(spatialEco)
library(rnaturalearth)


# import occurrence data (already curated)
load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")

# import ecoregions maps (version Dinerstein et al 2017)
ecoregions <- vect('C:/Users/user/Desktop/Ecoregions2017/Ecoregions2017.shp')

# coastlines
coastlines <- ne_coastline(scale = 110, returnclass = "sf") %>% vect()

# import predictors
agg_predictors <- rast('SDM/agg_predictors_MaxEnt.tiff')[[c("bio_1","bio_5","bio_7","bio_12","bio_13", "bio_15", "bio_18", "elevation", "treecover")]]
names(agg_predictors)
# cellSize(agg_predictors)


# vector of species with at least 15 observations
n_obs <- table(data_final[['occurrences']]$scientific_name) %>% as.data.frame() %>%
  subset(Freq>=15)
sdm_species <- n_obs[order(n_obs$Freq, decreasing=F),] %>% dplyr::select(Var1) %>% deframe() %>% droplevels()


# loop to model distributions and obtain sites with high habitat suitability for posterior analyses
suitable_locations <- list()
allmodels <- list()
summary_bestmodels <- data.frame(scientific_name=sdm_species,
                         parameters=NA,
                         auc.diff.avg=NA,
                         auc.val.avg=NA,
                         or.mtp.avg=NA)

for (s1 in 1:nrow(summary_bestmodels)) {
  
  spp1 <- summary_bestmodels$scientific_name[s1]
  
  # occurrences
  occ_temp <- data_final[['occurrences']] %>% filter(scientific_name==spp1)
  
  # points
  pts_temp <- vect(occ_temp, geom=c('x','y'), 'epsg:4326')
  
  # model IN and NEAR ecoregions within its EOO
  ext_temp <- convHull(pts_temp)
  eco_temp <- ecoregions %>% terra::mask(ext_temp) %>% as.data.frame() %>% dplyr::select(ECO_NAME) %>%
    unique() %>% deframe()
  eco_temp <- subset(ecoregions, ecoregions$ECO_NAME %in% eco_temp)
  
  # crop predictors
  temp_predictors <- agg_predictors %>% crop(eco_temp, mask=T)
  

  # plot(temp_predictors$elevation)
  # lines(ext_temp, col='red')
  # lines(eco_temp)
  # points(pts_temp, col='blue')
  
  
  # Bias file: Target-group background (https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13442)
  target_points <- data_final[['occurrences']][,c('x','y')] %>% unique() %>% 
    vect(geom=c('x','y'), 'epsg:4326') %>%
    intersect(eco_temp) # all Scleria observations for our study area
  
  # 2d kernel density estimation 
  target_density <- kde(geom(target_points)[,c('x','y')]) 
  
  # resample + crop to ground
  target_raster <- raster(target_density) %>% rast() %>%
    resample(temp_predictors$bio_1, method='bilinear') %>%
    crop(eco_temp, mask=T) 
  
  # normalize
  target_raster <- target_raster - min(target_raster[], na.rm=T)
  target_raster <- raster.transformation(target_raster, trans="norm") 

  # plot(target_raster)
  # points(target_points, col='red')
  # lines(eco_temp)
  # points(pts_temp, col='blue')
  
  
  # select a maximum of 10000 background points
  bg_temp <- as.data.frame(target_raster, xy=T) %>% subset(layer>0.1)
  if (nrow(bg_temp) > 10000) {
    bg_temp <- bg_temp[sample(1:nrow(bg_temp), prob=bg_temp$layer, size=10000, replace=F),]
  } else {
    bg_temp
  }

  
  # convert to Raster* file to run ENMevaluate
  stk_predictors <- stack(temp_predictors) # Raster* file
  
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
                         parallel = TRUE, numCores = NULL)
  
  
  # save results and select best model
  results1 <- modeval@results
  results1$nfeat <- nchar(as.character(results1$fc)) # number of features
  allmodels[[sdm_species[s1]]] <- results1

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
  map1 <- rast(modeval@predictions[[as.numeric(rownames(bestmod1))]])

  
  # plot(map1, main="Model with lowest AICc")
  # points(pts_temp, col='blue')
  
  
  # use 10th percentile as threshold
  q10 <- map1 %>% terra::extract(pts_temp, ID=F) %>% deframe() %>% quantile(probs=0.1)
  map1[map1>=q10] <- 1; map1[map1<q10] <- 0
  
  
  # plot map
  mypath <- file.path(paste(getwd(),'/SDM/maps/',spp1,'.tiff',sep=''))
  tiff(file=mypath, width=600, height=600, units='px')
  plot(map1, main=spp1, legend=F)
  points(pts_temp, col='black')
  lines(coastlines)
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
  
  print(spp1)

}


