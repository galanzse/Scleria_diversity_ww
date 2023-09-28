

library(tidyverse)
library(ggpubr)
library(terra)


# import data
load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")

# occurrences
occ_scleria <- data_final[['occurrences']]
pts_scleria <- occ_scleria %>% vect(geom=c('x','y'), 'epsg:4326') %>% project('+proj=eqearth')


# reference raster
mygrid <- rast('C:/Users/user/Desktop/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif') %>% # 7 km2
  project('+proj=eqearth')

cellSize(mygrid, unit='km')


# loop
optim_res <- matrix(nrow=100, ncol=4) %>% as.data.frame()
colnames(optim_res) <- c('fact', 'mean_area', 'n_occupied_cells', 'mean_spp_cell')
optim_res$fact <- 2:101

# file to retain species ID
occ_temp <- occ_scleria

for (r in 1:nrow(optim_res)) {
  
  # empty cell IDs
  occ_temp$cell <- NA
  
  # aggregate
  grd1 <- mygrid %>% aggregate(fact=optim_res$fact[r], na.rm=TRUE) 
  
  # extract cell IDs
  occ_temp$cell <- grd1 %>% extract(pts_scleria, cells=TRUE) %>% dplyr::select(cell) %>% deframe()
  
  # calculate average area
  optim_res$mean_area[r] <- grd1 %>% cellSize(transform=TRUE, unit='km') %>% as.data.frame() %>% colMeans() %>% sqrt()
  
  # calculate number of occupied cells
  optim_res$n_occupied_cells[r] <- occ_temp$cell %>% unique() %>% length()
  
  # calculate average species richness
  optim_res$mean_spp_cell[r] <- occ_temp %>% group_by(cell) %>% summarise(spp_cell=n_distinct(scientific_name)) %>% dplyr::select(spp_cell) %>% colMeans()
  
  # progress
  print(r)
  
}

# write.table(optim_res, 'results/optim_res.txt')


# results
g1 <- ggplot(aes(x=mean_area, y=n_occupied_cells), data=optim_res) +
  geom_point() + theme_classic() +
  xlab('Cell size (km)') + ylab('Number of occupied cells') +
  theme( panel.grid.major = element_line(colour = "grey", size=.5),
         panel.grid.minor = element_line(colour = "grey", size=.5) )
  
g2 <- ggplot(aes(x=mean_area, y=mean_spp_cell), data=optim_res) +
  geom_point() + theme_bw() +
  xlab('Cell size (km)') + ylab('Mean richness') +
  theme( panel.grid.major = element_line(colour = "grey", size=.5),
         panel.grid.minor = element_line(colour = "grey", size=.5) )

ggarrange(g1, g2)


