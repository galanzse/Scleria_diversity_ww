

# CALCULATE ALPHA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY 


library(tidyverse)

library(rnaturalearth)
library(terra)

library(FD) # Fdisp
library(randtip) # imputation
library(picante) # pd (Faith), mpd



# import data ####

load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")

# occurrences
scl_occurrences <- data_final[['occurrences']]
scl_points <- scl_occurrences %>% vect(geom=c('x','y'), 'epsg:4326') %>% project('+proj=eqearth')

scl_taxonomy <- data_final[['assessments']] %>% filter(!(is.na(section))) # infrageneric classification, treat sections as genera for imputation


# reference raster (cell size previously optimised)
world_lines <- ne_countries(scale=10, type="countries", continent=NULL, country=NULL, geounit=NULL, sovereignty=NULL, returnclass="sf") %>% terra::vect() %>% terra::project('+proj=eqearth')

world_grid <- rast('C:/Users/user/Desktop/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif') %>%
  project('+proj=eqearth') %>%
  aggregate(fact=37, fun='modal', na.rm=TRUE) %>% # approx. 100km 
  crop(world_lines, touches=T, mask=T)

world_grid[] <- NA

names(world_grid) <- 'reference_grid'

cellSize(world_grid, unit='km')

plot(world_grid, col='grey', legend=NULL)
lines(world_lines)
points(scl_points, pch=1)


# null models
regional_pools <- read.csv("results/regional_pools.txt", sep="")
load("results/global_pool_list.Rdata")
load("results/regional_pool_list.Rdata")



# traits ####

scl_traits <- data_final[['traits']]

for (i in 1:nrow(scl_traits)) {
  for (t in c('life_form_simp', 'height', 'blade_area', 'nutlet_volume')) {
    if (is.na(scl_traits[i,t])) {
      temp <- scl_traits %>% subset(section==scl_traits$section[i] & life_form_simp==scl_traits$life_form_simp[i]) %>%
        dplyr::select(all_of(t)) %>% colMeans(na.rm=T)
      scl_traits[i,t] <- temp
    }
  }
} # impute by section x life form means

scl_traits[scl_traits=='NaN'] <- NA # correct format

scl_traits[scl_traits$scientific_name=='Scleria khasiana Boeckeler',c('height', 'blade_area', 'nutlet_volume')] <- scl_traits %>% subset(section=='Elatae') %>% dplyr::select(height, blade_area, nutlet_volume) %>% colMeans(na.rm=T) # S. khasiana is imputed from section data only

# transform and scale
scl_traits$height <- scl_traits$height %>% log() %>% scale(center=F)
scl_traits$blade_area <- scl_traits$blade_area %>% log() %>% scale(center=F)
scl_traits$nutlet_volume <- scl_traits$nutlet_volume %>% log() %>% scale(center=F)

# hist(scl_traits$height)
# hist(scl_traits$blade_area)
# hist(scl_traits$nutlet_volume)

# prepare matrix for analyses
scl_traits <- scl_traits %>% dplyr::select(scientific_name, life_form_simp, height, blade_area, nutlet_volume)
rownames(scl_traits) <- scl_traits$scientific_name

# species to work with (present in all datasets and with infrageneric information)
v_spp <- intersect(intersect(scl_occurrences$scientific_name, scl_traits$scientific_name), scl_taxonomy$scientific_name)
scl_occurrences <- scl_occurrences %>% filter(scientific_name %in% v_spp) # filter occurrences without trait data
scl_points <- scl_occurrences %>% vect(geom=c('x','y'), 'epsg:4326') %>% project('+proj=eqearth') # fix
scl_traits <- scl_traits %>% filter(scientific_name %in% v_spp) # fix
scl_traits$scientific_name <- NULL
scl_traits$life_form_simp <- as.numeric(as.factor(scl_traits$life_form_simp)) # fix



# create X imputed trees using randtip and save in list ####
# evolutionary diversity needs to be calculated considering uncertainty (40% species imputed)

scl_traits$scientific_name[is.na(scl_traits$height)]

scl_phylogeny <- data_final[['phylogeny']] # phylogeny

scl_phylogeny <- drop.tip(scl_phylogeny, scl_phylogeny$tip.label[which(!(scl_phylogeny$tip.label %in% scl_occurrences$scientific_name))]) # drop species without occurrences

scl_taxonomy <- data_final[['assessments']] %>% filter(scientific_name %in% v_spp) # infrageneric classification, treat sections as genera for imputation
scl_taxonomy$sectxspp <- paste(scl_taxonomy$section, str_split(scl_taxonomy$scientific_name, pattern=" ", simplify = TRUE)[,2], sep=' ')

phylogeny_sect <- scl_phylogeny # change names in tree
intree <- scl_taxonomy[scl_taxonomy$scientific_name %in% phylogeny_sect$tip.label,]
phylogeny_sect$tip.label <- intree$sectxspp[order(match(intree$scientific_name, phylogeny_sect$tip.label))]

back.tree <- phylogeny_sect # backbone phylogeny
class(back.tree)
is.ultrametric(back.tree)

sp.list <- scl_taxonomy$sectxspp # species list

# build info df, check and define ranks
my.info.noranks <- randtip::build_info(species=sp.list, tree=back.tree, mode='backbone', find.ranks=FALSE)
my.check <- randtip::check_info(my.info.noranks, tree=back.tree)
my.input.noranks <- randtip::info2input(my.info.noranks, back.tree)

imputed_trees <- list()
nr=500
scl_taxonomy$sectxspp <- gsub(' ', '_', scl_taxonomy$sectxspp) # to change to scientific names after imputation
for (i in 1:nr) {
  imp_tree <- randtip::rand_tip(input=my.input.noranks, tree=back.tree,
                                rand.type='random', # default
                                respect.mono=T, # all our clades are monophyletic so not relevant
                                prob=F, # branch selection probability is equiprobable
                                use.stem=T, # the stem branch can be considered as candidate for binding
                                prune=F, # we need the entire tree to run EDGE2
                                verbose=F)
  
  # (!!!) Nnode(imp_tree) == length(imp_tree$tip.label)
  
  imp_tree$tip.label <- scl_taxonomy$scientific_name[order(match(scl_taxonomy$sectxspp, imp_tree$tip.label))] # retrieve original names
  
  imputed_trees[[i]] <- imp_tree # save in list
  ck_n <- deframe(table(imp_tree$tip.label%in%scl_taxonomy$scientific_name))
  print(paste('x',i,'   ', 'Ntip=', ck_n, sep='')) # progress
}

save(imputed_trees, file="results/imputed_trees.RData")



# compute alpha diversity and trait means ####

# unique cells
scl_occurrences$cell <- extract(world_grid, scl_points, cells=TRUE) %>% dplyr::select(cell) %>% deframe()

# indices dataframe
slc_indices <- scl_occurrences %>% dplyr::select(cell) %>% unique()
slc_indices$taxon_richness <- NA # taxonomic richness
slc_indices$func_richness <- NA # Convex Hull
slc_indices$func_dispersion <- NA # Fdisp
slc_indices$perenn_prop <- NA #
slc_indices$cwm_height <- NA #
slc_indices$cwm_blade <- NA #
slc_indices$cwm_nutlet <- NA #
slc_indices$phylo_richness <- NA # Faith
slc_indices$phylo_diversity <- NA # mpd

# maps
map_richness <- world_grid; names(map_richness) <- 'richness'
map_func_richness <- world_grid; names(map_func_richness) <- 'F_rich'
map_func_dispersion <- world_grid; names(map_func_dispersion) <- 'F_disp'
map_properenn <- world_grid; names(map_properenn) <- 'prop_perennials'
map_cwm_height <- world_grid; names(map_cwm_height) <- 'CWM_height'
map_cwm_blade <- world_grid; names(map_cwm_blade) <- 'CWM_bladearea'
map_cwm_nutlet <- world_grid; names(map_cwm_nutlet) <- 'CWM_nutletvolume'
map_phylo_richness <- world_grid; names(map_phylo_richness) <- 'Faith'
map_phylo_diversity <- world_grid; names(map_phylo_diversity) <- 'mpd'

# site x species presence matrix to subset communities from
scl_presences <- matrix(ncol=length(scl_traits$scientific_name), nrow=1)
colnames(scl_presences) <- scl_traits$scientific_name
scl_presences[1,] <- 1

# loop
for (c in 1:nrow(slc_indices)) {
  
  cell_spp <- scl_occurrences %>% subset(cell==slc_indices$cell[c]) %>%
    dplyr::select(scientific_name) %>% unique() %>% deframe()

  # taxonomic richness
  slc_indices$taxon_richness[c] <- length(cell_spp)

  # presence matrix
  cell_comm <- t(as.matrix(scl_presences[1,cell_spp]))

  # FRic
  scl_dist <- scl_traits %>% subset(scientific_name%in%cell_spp) %>% select(-scientific_name)
  scl_dist$life_form_simp <- as.numeric(as.factor(scl_dist$life_form_simp))
  scl_dist <- daisy(as.matrix(scl_dist), metric='gower',
                    type=list('factor'='life_form_simp','numeric'=c('height','blade_area','nutlet_volume')))
  tempFD <- dbFD(x=scl_dist, stand.FRic=F, messages=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
  Matnull$func_richness[n] <- tempFD$FRic
  Matnull$func_dispersion[n] <- tempFD$FDis
  

  # life form
  temp_lifeform <- dummy_traits[cell_spp,c('annual','perennial')]
  if (class(temp_lifeform)[1]=='numeric') {
    slc_indices$perenn_prop[c] <- temp_lifeform['perennial'] / c(temp_lifeform['perennial'] + temp_lifeform['annual'])
    slc_indices$annual_prop[c] <- temp_lifeform['annual'] / c(temp_lifeform['perennial'] + temp_lifeform['annual'])
  } else {
    temp_lifeform <- colSums(temp_lifeform)
    slc_indices$perenn_prop[c] <- temp_lifeform['perennial'] / c(temp_lifeform['perennial'] + temp_lifeform['annual'])
    slc_indices$annual_prop[c] <- temp_lifeform['annual'] / c(temp_lifeform['perennial'] + temp_lifeform['annual'])
  }
  
  # CWM
  cell_CWM <- dummy_traits[cell_spp,c('height', 'blade_area', 'nutlet_volume')]
  if(class(cell_CWM)[1]=='numeric') {
    slc_indices$cwm_height[c] <- mean(cell_CWM['height'])
    slc_indices$cwm_blade[c] <- mean(cell_CWM['blade_area'])
    slc_indices$cwm_nutlet[c] <- mean(cell_CWM['nutlet_volume'])
  } else {
    slc_indices$cwm_height[c] <- mean(cell_CWM[,'height'])
    slc_indices$cwm_blade[c] <- mean(cell_CWM[,'blade_area'])
    slc_indices$cwm_nutlet[c] <- mean(cell_CWM[,'nutlet_volume'])
    }

  # # calculate average Faith and mpd from imputed trees
  # # c2=length(imputed_trees)
  c2=150
  v_faith <- vector(length=c2); v_faith[] <- NA
  v_mpd <- vector(length=c2)
  for (c3 in 1:c2) {
    v_faith[c3] <- pd(cell_comm, imputed_trees[[c3]], include.root=TRUE)[1,'PD']
    v_mpd[c3] <- mpd(cell_comm, cophenetic(imputed_trees[[c3]]), abundance.weighted=FALSE)
  }
  slc_indices$phylo_richness[c] <- mean(v_faith)
  slc_indices$phylo_diversity[c] <- mean(v_mpd, na.rm=T)


  # maps
  map_richness[slc_indices$cell[c]] <- slc_indices$taxon_richness[c]
  map_func_richness[slc_indices$cell[c]] <- slc_indices$func_richness[c]
  map_func_dispersion[slc_indices$cell[c]] <- slc_indices$func_dispersion[c]
  map_properenn[slc_indices$cell[c]] <- slc_indices$perenn_prop[c]
  map_propannual[slc_indices$cell[c]] <- slc_indices$annual_prop[c]
  map_cwm_height[slc_indices$cell[c]] <- slc_indices$cwm_height[c]
  map_cwm_blade[slc_indices$cell[c]] <- slc_indices$cwm_blade[c]
  map_cwm_nutlet[slc_indices$cell[c]] <- slc_indices$cwm_nutlet[c]
  map_phylo_richness[slc_indices$cell[c]] <- slc_indices$phylo_richness[c]
  map_phylo_diversity[slc_indices$cell[c]] <- slc_indices$phylo_diversity[c]
  
  
  # progress
  print(paste('--- ', round(c/nrow(slc_indices)*100, 2), '% ---', sep=''))
  
}


# slc_indices[is.na(slc_indices)] <- 0
# write.table(slc_indices, 'results/slc_indices.txt')

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


