

# Function optimizeRR: explore relationships between raster cell area, prevalence, richness and Moran's I
# optimizeRR takes a raster and species occurrences to compute several metrics after applying different aggregation factors


library(tidyverse)
library(ggpubr)
library(terra)
library(ape) # moran.i
library(geosphere) # distm
library(GeoThinneR) # thin_points


# function autocor to calculate Moran's I
url <- "https://raw.githubusercontent.com/galanzse/Scleria_diversity_ww/refs/heads/master/scripts/fun_autocor.R"
source(url)


# Arguments:
# raster - SpatRaster in epsg:4326 with the maximum resolution that wants to be considered for analysis
# occurrences - dataframe of species occurrences (y/x) 
# max_agg - Maximum aggregation factor. All integer values between 1 and max_agg will be used as aggregation factors in terra::aggregate
# aoi.I - For large datasets. If TRUE, Moran.I will be computed for a given region which has its centroid in cen.I c(lat,lon) and a width of buff.i (in kms)


optimizeRR <- function(raster=NULL, occurrences=NULL, max_agg=NULL, aoi.I=NULL, cen.I=NULL, buff.i=NULL) {

  # check data format
  occurrences <- as.data.frame(occurrences)
  if (ncol(occurrences)!=3) stop("Occurrence dataframe should have 3 columns only: 'species', 'x' and 'y'")
  if (all(colnames(occurrences) %in% c('x','y','species'))==FALSE) stop("Colnames are not 'species', 'x', 'y'")
  
  # project raster to equal area
  raster_ll <- raster
  raster_ee <- raster %>% terra::project('+proj=eqearth')
  
  # project occurrences same as raster
  occurrences_ee <- occurrences %>% terra::vect(geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth')
  
  print('-- Raster and observations projected to eqearth --')
  
  # loop
  optimizeRR_out <- matrix(nrow=max_agg, ncol=8) %>% as.data.frame()
  colnames(optimizeRR_out) <- c('fact', 'mean_area', 'sqrt_area',
                                'n_occ_cells', 'prop_occ_cells',
                                'mean_richness', 'sd_richness',
                                'I_Moran')
  optimizeRR_out$fact <- 1:max_agg
  
  for (r in 1:nrow(optimizeRR_out)) {

    # aggregate
    agg_raster_ll <- raster_ll %>% terra::aggregate(fact=optimizeRR_out$fact[r], na.rm=TRUE)
    agg_raster_ee <- raster_ee %>% terra::aggregate(fact=optimizeRR_out$fact[r], na.rm=TRUE)
    
    # calculate average area
    optimizeRR_out$mean_area[r] <- agg_raster_ee %>% terra::cellSize(transform=TRUE, unit='km') %>%
      as.data.frame() %>% colMeans()

    print(paste('-- Raster aggregated! Cell area is', round(optimizeRR_out$mean_area[r], 2), 'km2 --'))
    
    # calculate raster cell side length
    optimizeRR_out$sqrt_area[r] <- sqrt(optimizeRR_out$mean_area[r])
    
    # extract cell ID occupied per species
    temp <- agg_raster_ee %>% terra::extract(occurrences_ee, cells=T, xy=T) %>% dplyr::select(cell, x, y)
    temp$species <- occurrences[,-which(colnames(occurrences)%in%c('x','y'))] ### CHECK THISSSS

    # calculate number of occupied cells
    optimizeRR_out$n_occ_cells[r] <- temp %>% dplyr::select(cell) %>% unique() %>% nrow()
    
    # calculate proportion of occupied cells
    optimizeRR_out$prop_occ_cells[r] <- optimizeRR_out$n_occ_cells[r]/length(unique(cells(agg_raster_ee)))

    # calculate mean species richness
    optimizeRR_out$mean_richness[r] <- temp %>% group_by(cell) %>%
      summarise(temp=n_distinct(species)) %>% dplyr::select(temp) %>% deframe() %>% mean()
    optimizeRR_out$sd_richness[r] <- temp %>% group_by(cell) %>%
      summarise(temp=n_distinct(species)) %>% dplyr::select(temp) %>% deframe() %>% sd()
    
    # Estimate Moran.I for the entire dataset (retain one observation per species and cell)
    if (aoi.I==TRUE) {
      temp <- data.frame(y=cen.I[1], x=cen.I[2]) %>%
        terra::vect(geom=c('x','y'), crs='epsg:4326') %>% terra::buffer(width=buff.i*1000)
      temp <- occurrences[,c('x','y')] %>% vect(geom=c('x','y'), crs='epsg:4326') %>% terra::crop(temp) %>%
        terra::geom() %>% as.data.frame() %>% dplyr::select(x,y)
    } else {
      temp <- occurrences[,c('x','y')]
    }

    colnames(temp)[colnames(temp)=='y'] <- 'lat'; colnames(temp)[colnames(temp)=='x'] <- 'lon'
    temp <- temp[thin_points(temp, method="distance", thin_dist=optimizeRR_out$sqrt_area[r])$retained[[1]],]
    colnames(temp)[colnames(temp)=='lon'] <- 'x'; colnames(temp)[colnames(temp)=='lat'] <- 'y'
    optimizeRR_out$I_Moran[r] <- autocor(agg_raster_ll, temp)$observed
    
    # progress
    print(paste('-- ', round(r/nrow(optimizeRR_out),2)*100, '% --' , sep=''))
    
  }
  
  write.table(optimizeRR_out, 'optimizeRR_out.txt')
  print(paste('-- optimizeRR_out saved in directory', getwd(), '--'))

  
  
  # results
  g1 <- ggplot(aes(x=sqrt_area, y=n_occ_cells), data=optimizeRR_out) +
    geom_point() + theme_classic() +
    xlab('Cell side length (km)') + ylab('Number of occupied cells')
  
  g2 <- ggplot(aes(x=sqrt_area, y=prop_occ_cells), data=optimizeRR_out) +
    geom_point() + theme_classic() +
    xlab('Cell side length (km)') + ylab('Proportion of occupied cells')
  
  g3 <- ggplot(aes(x=sqrt_area, y=mean_richness), data=optimizeRR_out) +
    geom_point() +
    geom_pointrange(aes(ymin = mean_richness - sd_richness, ymax = mean_richness + sd_richness)) +
    theme_classic() +
    xlab('Cell side length (km)') + ylab('Richness (mean + sd)')
  
  g4 <- ggplot(aes(x=sqrt_area, y=I_Moran), data=optimizeRR_out) +
    geom_point() + theme_classic() +
    xlab('Cell side length (km)') + ylab('Moran I')
  

  print(ggarrange(g1, g2, g3, g4))
  
  png(file='optimizeRR_out.png')
  plot(ggarrange(g1, g2, g3, g4), main="optimizeRR_out", quality=100,res=200)
  dev.off()
  
  print(paste('-- optimizeRR_out.png saved in directory', getwd(), '--'))


}


