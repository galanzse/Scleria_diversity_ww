


# CALCULATE ALPHA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY, AND SES


library(tidyverse)
library(rnaturalearth)
library(terra)
library(FD) # Fdisp
library(picante) # pd (Faith), mpd


# load data for analyses
source('scripts/import data alpha and beta diversity calc.R')


# null models and BIOMExREALM id
regional_pools <- read.csv("results/regional_pools.txt", sep="")
rast_ecoregions <- rast('results/rast_ecoregions.tiff')
load("results/global_pool_list.Rdata")
load("results/regional_pool_list.Rdata")



# compute alpha diversity and trait means ####

# unique cells
scl_occurrences$cell <- extract(world_grid, scl_points, cells=TRUE) %>% dplyr::select(cell) %>% deframe()
scl_occurrences$BIOMExREALM <- extract(rast_ecoregions, scl_points, cells=TRUE) %>% # info regional pools
  dplyr::select(REALMxBIOME) %>% deframe()
scl_indices <- scl_occurrences %>% dplyr::select(cell, BIOMExREALM) %>% unique()

table(is.na(scl_indices$BIOMExREALM)) # NAs in regional pool
scl_indices$BIOMExREALM[which(!(scl_indices$BIOMExREALM %in% ecoregions$REALMxBIOME))] %>% unique()
scl_indices$BIOMExREALM[scl_indices$BIOMExREALM %in% c('AT98','AA14','NA2','PA13','NA98')] <- NA

# site x species presence matrix to subset communities from
scl_presences <- matrix(ncol=nrow(scl_traits), nrow=1)
colnames(scl_presences) <- rownames(scl_traits)
scl_presences[1,] <- 1


# indices
scl_indices$richness <- NA # taxonomic richness

scl_indices$ann_per <- NA #
scl_indices$cwm_height <- NA # CWM
scl_indices$cwm_blade <- NA #
scl_indices$cwm_nutlet <- NA #

scl_indices$Frich_ <- NA # Convex Hull    OBSERVED
scl_indices$Fdisp_ <- NA # Fdisp
scl_indices$Faith_ <- NA # Faith
scl_indices$mpd_ <- NA # mpd

scl_indices$SES_Frich_gl <- NA #    SES GLOBAL
scl_indices$SES_Fdisp_gl <- NA #
scl_indices$SES_Faith_gl <- NA #
scl_indices$SES_mpd_gl <- NA #

scl_indices$SES_Frich_rg <- NA #    SES REGIONAL
scl_indices$SES_Fdisp_rg <- NA #
scl_indices$SES_Faith_rg <- NA #
scl_indices$SES_mpd_rg <- NA #


# loop
for (cl in 1:nrow(scl_indices)) {
  
  cell_spp <- scl_occurrences %>% subset(cell==scl_indices$cell[cl]) %>%
    dplyr::select(scientific_name) %>% unique() %>% deframe()
  
  # presence matrix
  cell_comm <- t(as.matrix(scl_presences[1,cell_spp]))

  
  # taxonomic richness
  scl_indices$richness[cl] <- length(cell_spp)


  # life form: from -1 (100% annual), 0 (50% each life form), +1 (100% perennial)
  temp_lifeform <- scl_traits[cell_spp,'life_form_simp']
  if (length(unique(temp_lifeform))==1) {
    if (unique(temp_lifeform)==1) {
      scl_indices$ann_per[cl] <- c(-1)
    } else {
      scl_indices$ann_per[cl] <- 1
    }
  } else {
    lf_temp <- table(temp_lifeform==1)
    if (lf_temp['TRUE']==lf_temp['FALSE']) {
      scl_indices$ann_per[cl] <- 0
    } else {
        if (lf_temp['TRUE']>lf_temp['FALSE']) { scl_indices$ann_per[cl] <- lf_temp['TRUE']/sum(lf_temp)*c(-1) }
        if (lf_temp['TRUE']<lf_temp['FALSE']) { scl_indices$ann_per[cl] <- lf_temp['FALSE']/sum(lf_temp) }
      }
  }
  
  
  # CWM
  cell_CWM <- scl_traits[cell_spp,c('height', 'blade_area', 'nutlet_volume')]
  scl_indices$cwm_height[cl] <- mean(cell_CWM$height)
  scl_indices$cwm_blade[cl] <- mean(cell_CWM$blade_area)
  scl_indices$cwm_nutlet[cl] <- mean(cell_CWM$nutlet_volume)


  # Convex Hull & Fdisp
  scl_dist <- scl_traits[cell_spp,]
  if (length(cell_spp)>1) {
      scl_dist <- daisy(as.matrix(scl_dist), metric='gower',
                    type=list('factor'='life_form_simp','numeric'=c('height','blade_area','nutlet_volume')))
      tempFD <- dbFD(x=scl_dist, stand.FRic=F, messages=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
      scl_indices$Frich_[cl] <- tempFD$FRic
      scl_indices$Fdisp_[cl] <- tempFD$FDis
  }
  
  # # calculate average Faith and mpd from imputed trees
  # # Xi=length(imputed_trees)
  Xi=100
  v_faith <- vector(length=Xi); v_faith[] <- NA
  v_mpd <- vector(length=Xi)
  for (p1 in 1:Xi) {
    v_faith[p1] <- pd(cell_comm, imputed_trees[[p1]], include.root=TRUE)[1,'PD']
    v_mpd[p1] <- mpd(cell_comm, cophenetic(imputed_trees[[p1]]), abundance.weighted=FALSE)
  }
  scl_indices$Faith_[cl] <- median(v_faith)
  scl_indices$mpd_[cl] <- median(v_mpd, na.rm=T)
  
  
  
  # standardized effect size: GLOBAL POOL
  if (length(cell_spp)>1) {
    pool1 <- global_pool_list[[as.character(scl_indices$richness[cl])]]
    percentile <- ecdf(pool1$Frich_); scl_indices$SES_Frich_gl[cl] <- percentile(scl_indices$Frich_[cl])
    percentile <- ecdf(pool1$Fdisp_); scl_indices$SES_Fdisp_gl[cl] <- percentile(scl_indices$Fdisp_[cl])
    percentile <- ecdf(pool1$Faith_); scl_indices$SES_Faith_gl[cl] <- percentile(scl_indices$Faith_[cl])
    percentile <- ecdf(pool1$mpd_); scl_indices$SES_mpd_gl[cl] <- percentile(scl_indices$mpd_[cl])
    
    
    # standardized effect size: REGIONAL POOL
    if (!(is.na(scl_indices$BIOMExREALM[cl]))) {
      pool1 <- regional_pool_list[[as.character(scl_indices$BIOMExREALM[cl])]]
      pool1 <- pool1[[as.character(scl_indices$richness[cl])]]
      percentile <- ecdf(pool1$Frich_); scl_indices$SES_Frich_rg[cl] <- percentile(scl_indices$Frich_[cl])
      percentile <- ecdf(pool1$Fdisp_); scl_indices$SES_Fdisp_rg[cl] <- percentile(scl_indices$Fdisp_[cl])
      percentile <- ecdf(pool1$Faith_); scl_indices$SES_Faith_rg[cl] <- percentile(scl_indices$Faith_[cl])
      percentile <- ecdf(pool1$mpd_); scl_indices$SES_mpd_rg[cl] <- percentile(scl_indices$mpd_[cl])
    }
  }
  
  

  # progress
  print(paste('--- ', round(cl/nrow(scl_indices)*100, 2), '% ---', sep=''))
  
}


scl_indices[is.na(scl_indices)] <- 0
write.table(scl_indices, 'results/scl_indices.txt')



# maps ####


map_richness <- world_grid; names(map_richness) <- 'richness'
map_func_richness <- world_grid; names(map_func_richness) <- 'F_rich'
map_func_dispersion <- world_grid; names(map_func_dispersion) <- 'F_disp'
map_properenn <- world_grid; names(map_properenn) <- 'prop_perennials'
map_cwm_height <- world_grid; names(map_cwm_height) <- 'CWM_height'
map_cwm_blade <- world_grid; names(map_cwm_blade) <- 'CWM_bladearea'
map_cwm_nutlet <- world_grid; names(map_cwm_nutlet) <- 'CWM_nutletvolume'
map_phylo_richness <- world_grid; names(map_phylo_richness) <- 'Faith'
map_phylo_diversity <- world_grid; names(map_phylo_diversity) <- 'mpd'


map_richness[scl_indices$cell[cl]] <- scl_indices$richness[cl]
map_func_richness[scl_indices$cell[cl]] <- scl_indices$Frich_[cl]
map_func_dispersion[scl_indices$cell[cl]] <- scl_indices$Fdisp_[cl]
map_properenn[scl_indices$cell[cl]] <- scl_indices$ann_per[cl]
map_propannual[scl_indices$cell[cl]] <- scl_indices$annual_prop[cl]
map_cwm_height[scl_indices$cell[cl]] <- scl_indices$cwm_height[cl]
map_cwm_blade[scl_indices$cell[cl]] <- scl_indices$cwm_blade[cl]
map_cwm_nutlet[scl_indices$cell[cl]] <- scl_indices$cwm_nutlet[cl]
map_phylo_richness[scl_indices$cell[cl]] <- scl_indices$Faith_[cl]
map_phylo_diversity[scl_indices$cell[cl]] <- scl_indices$mpd_[cl]

alpha_div_rasters <- c(map_richness, map_func_richness, map_func_dispersion, map_properenn, map_propannual, map_cwm_height, map_cwm_blade, map_cwm_nutlet, map_phylo_richness, map_phylo_diversity)
writeRaster(alpha_div_rasters, 'results/alpha_div_rasters.tiff', overwrite=TRUE)


# par(mfrow=c(1,1))
plot(map_richness, main='Taxonomic richness'); lines(world_lines)
plot(map_func_richness, main='Functional richness'); lines(world_lines)
plot(map_func_dispersion, main='Functional dispersion'); lines(world_lines)
plot(map_phylo_richness, main='Faith`s phylogenetic diversity'); lines(world_lines)
plot(map_phylo_diversity, main='Phylogenetic mean pairwise distance'); lines(world_lines)
plot(map_cwm_height, main='Community weighted mean Height'); lines(world_lines)
plot(map_cwm_blade, main='Community weighted mean Blade area'); lines(world_lines)
plot(map_cwm_nutlet, main='Community weighted mean Nutlet volume'); lines(world_lines)
plot(map_propannual, main='Proportion of annual species'); lines(world_lines)
plot(map_properenn, main='Proportion of perennial species'); lines(world_lines)


