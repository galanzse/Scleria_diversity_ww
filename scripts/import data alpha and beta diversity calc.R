

# PREPARE DATA FOR THE CALCULATION OF ALPHA AND BETA DIVERSITY, INCLUDING NULL MODELS


library(tidyverse)
library(rnaturalearth)
library(cluster)
library(terra)
library(randtip) # imputation



load("results/data_final.RData")

# occurrences
scl_occurrences <- data_final[['occurrences']]
scl_points <- scl_occurrences %>% vect(geom=c('x','y'), 'epsg:4326') %>% project('+proj=eqearth')

scl_taxonomy <- data_final[['assessments']] %>% filter(!(is.na(section))) # infrageneric classification, treat sections as genera for imputation


# reference raster (cell size previously optimised)
world_lines <- ne_countries(scale=10, type="countries", continent=NULL, country=NULL, geounit=NULL, sovereignty=NULL, returnclass="sf") %>% terra::vect() %>% terra::project('+proj=eqearth')

world_grid <-
  # rast('C:/Users/user/Desktop/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif') %>%
  rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m_elev.tif') %>%
  project('+proj=eqearth') %>%
  aggregate(fact=55, fun='modal', na.rm=TRUE) %>% # approx. 150km 
  crop(world_lines, touches=T, mask=T)

world_grid[] <- NA

names(world_grid) <- 'reference_grid'

# cellSize(world_grid, unit='km')
# 
# plot(world_grid, col='grey', legend=NULL)
# lines(world_lines)
# points(scl_points, pch=1)



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


# We will use a dendrogram so there are values of Frich for all species
str(scl_traits)
scl_dist <- daisy(as.matrix(scl_traits), metric="gower", type=list('factor'=1,'numeric'=2:4))
scl_dendrogram <- hclust(scl_dist, method="average") # clustering (UPGMA)


# create X imputed trees using randtip and save in list ####
# evolutionary diversity needs to be calculated considering uncertainty (40% species imputed)

scl_traits$scientific_name[is.na(scl_traits$height)]

scl_phylogeny <- data_final[['phylogeny']] # phylogeny

scl_phylogeny <- drop.tip(scl_phylogeny, scl_phylogeny$tip.label[which(!(scl_phylogeny$tip.label %in% scl_occurrences$scientific_name))]) # drop species without occurrences

scl_taxonomy <- data_final[['assessments']] %>% filter(scientific_name %in% v_spp) # infrageneric classification, treat sections as genera for imputation
# scl_taxonomy$sectxspp <- paste(scl_taxonomy$section, str_split(scl_taxonomy$scientific_name, pattern=" ", simplify = TRUE)[,2], sep=' ')
# 
# phylogeny_sect <- scl_phylogeny # change names in tree
# intree <- scl_taxonomy[scl_taxonomy$scientific_name %in% phylogeny_sect$tip.label,]
# phylogeny_sect$tip.label <- intree$sectxspp[order(match(intree$scientific_name, phylogeny_sect$tip.label))]
# 
# back.tree <- phylogeny_sect # backbone phylogeny
# class(back.tree)
# is.ultrametric(back.tree)
# 
# sp.list <- scl_taxonomy$sectxspp # species list
# 
# # build info df, check and define ranks
# my.info.noranks <- randtip::build_info(species=sp.list, tree=back.tree, mode='backbone', find.ranks=FALSE)
# my.check <- randtip::check_info(my.info.noranks, tree=back.tree)
# my.input.noranks <- randtip::info2input(my.info.noranks, back.tree)
# 
# imputed_trees <- list()
# nr=1000
# scl_taxonomy$sectxspp <- gsub(' ', '_', scl_taxonomy$sectxspp) # to change to scientific names after imputation
# for (i in 1:nr) {
#   imp_tree <- randtip::rand_tip(input=my.input.noranks, tree=back.tree,
#                                 rand.type='random', # default
#                                 respect.mono=T, # all our clades are monophyletic so not relevant
#                                 prob=F, # branch selection probability is equiprobable
#                                 use.stem=T, # the stem branch can be considered as candidate for binding
#                                 prune=F, # we need the entire tree to run EDGE2
#                                 verbose=F)
#   
#   # (!!!) Nnode(imp_tree) == length(imp_tree$tip.label)
#   
#   imp_tree$tip.label <- scl_taxonomy$scientific_name[order(match(scl_taxonomy$sectxspp, imp_tree$tip.label))] # retrieve original names
#   
#   imputed_trees[[i]] <- imp_tree # save in list
#   ck_n <- deframe(table(imp_tree$tip.label%in%scl_taxonomy$scientific_name))
#   print(paste('x',i,'   ', 'Ntip=', ck_n, sep='')) # progress
# }
# 
# save(imputed_trees, file="results/imputed_trees.RData")


load("results/imputed_trees.RData")


rm(i, t, temp)
# rm(i, t, temp, sp.list, imp_tree, ck_n, my.info.noranks, my.check,
#    my.input.noranks, phylogeny_sect, intree) #


