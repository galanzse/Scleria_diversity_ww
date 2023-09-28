

# Null models for alpha diversity indices (except for richness)


library(tidyverse)
library(cluster)
library(BAT) # convex.hull
library(FD) # Fdisp
library(randtip) # imputation
library(picante) # pd (Faith), mpd



# Xnull randomization: shuffling the taxon labels on the phylogeny & func.distance matrix + species richness constant

# we are going to define two pools:
  # global pool: any species can occur in any grid (Swenson et al. 2012)
  # regional pool: species ranges are constant (Spalink et al. 2018)



# global pool: the most efficient way to do this is to generate a distribution for each species richness value
scl_vrich <- read.csv("results/slc_indices.txt", sep="") %>% dplyr::select(taxon_richness) %>% deframe() %>% unique()
scl_vrich <- scl_vrich[-which(scl_vrich==1)] # at least two species

global_pool_list <- list()

# site x species presence matrix to subset communities from
scl_presences <- matrix(ncol=length(scl_traits$scientific_name), nrow=1)
colnames(scl_presences) <- scl_traits$scientific_name
scl_presences[1,] <- 1

Xnull=100

Matnull <- matrix(ncol=4, nrow=Xnull) %>% data.frame() # results for every combination
colnames(Matnull) <- c('func_richness','func_dispersion','phylo_richness','phylo_diversity')

for (r in 1:length(scl_vrich)) { # for every observed richness value
  
  Matnull[] <- NA
  
  for (n in 1:nrow(Matnull)) {  # create X random values
    
    cell_spp <- sample(rownames(dummy_traits), scl_vrich[r], replace=FALSE)
    
    # presence matrix
    cell_comm <- t(as.matrix(scl_presences[1,cell_spp]))
    
    # Convex Hull & Fdisp
    scl_dist <- scl_traits[cell_spp,]
    scl_dist <- daisy(as.matrix(scl_dist), metric='gower',
                      type=list('factor'='life_form_simp','numeric'=c('height','blade_area','nutlet_volume')))
    tempFD <- dbFD(x=scl_dist, stand.FRic=F, messages=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
    Matnull$func_richness[n] <- tempFD$FRic
    Matnull$func_dispersion[n] <- tempFD$FDis
    
    # calculate average Faith and mpd from imputed trees
    temp_tree <- imputed_trees[[sample(1:length(imputed_trees),1)]]
    Matnull$phylo_richness[n] <- pd(cell_comm, temp_tree, include.root=TRUE)[1,'PD']
    Matnull$phylo_diversity[n] <- mpd(cell_comm, cophenetic(temp_tree), abundance.weighted=FALSE)
    
  }
  
  # save
  global_pool_list[[scl_vrich[r]]] <- Matnull
  
  # progress
  print(paste('--- ', round(r/length(scl_vrich)*100, 2), '% ---', sep=''))
  
}

save(global_pool_list, file='results/global_pool_list.Rdata')



# regional pools: identify species that occur within the same BIOME in the same REALM
ecoregions <- read.csv("results/df_ecoregions.txt", sep="") %>%
  subset(scientific_name%in%v_scl_spp) %>%
  dplyr::select(scientific_name, REALM, BIOME) %>% unique() %>% na.omit()

regional_pools <- unique(regional_pool[,c('REALM','BIOME')]) # 39 potential pools, 17 species per pool on average
for (p in 1:nrow(spp_x_pool)) { spp_x_pool$n_species[p] <- regional_pool %>% subset(REALM==spp_x_pool$REALM[p] & BIOME==spp_x_pool$BIOME[p]) %>% dplyr::select(scientific_name) %>% deframe() %>% length() }
regional_pools <- spp_x_pool[spp_x_pool$n_species>2,] # remove combinations with 2 or less species
regional_pools$code <- 1:nrow(regional_pools) # add code to add list

write.table(regional_pools, 'results/regional_pools.txt')

regional_pool_list <- list()

Matnull <- matrix(ncol=4, nrow=Xnull) %>% data.frame() # results for every combination
colnames(Matnull) <- c('func_richness','func_dispersion','phylo_richness','phylo_diversity')

for (e in 1:nrow(regional_pools)) {
  
  reg1 <- regional_pools[e,]
  pool1 <- ecoregions %>% subset(BIOME==reg1$BIOME & REALM==reg1$REALM)
  pool1 <- pool1$scientific_name[which(pool1$scientific_name %in% rownames(dummy_traits))]
  
  reg_pool <- list()
  
  for (r in 1:length(scl_vrich)) { # for every observed richness value
  
    if (length(pool1) >= scl_vrich[r]) { # create null models up to the richness of the richest pixel
      
      Matnull[] <- NA
      
      for (n in 1:nrow(Matnull)) {  # create X random values
        
        cell_spp <- sample(pool1, scl_vrich[r], replace=FALSE)
        
        # presence matrix
        cell_comm <- t(as.matrix(scl_presences[,cell_spp]))
        
        # Convex Hull & Fdisp
        scl_dist <- scl_traits[cell_spp,]
        scl_dist <- daisy(as.matrix(scl_dist), metric='gower',
                          type=list('factor'='life_form_simp','numeric'=c('height','blade_area','nutlet_volume')))
        tempFD <- dbFD(x=scl_dist, stand.FRic=F, messages=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
        Matnull$func_richness[n] <- tempFD$FRic
        Matnull$func_dispersion[n] <- tempFD$FDis
        
        # calculate average Faith and mpd from imputed trees
        temp_tree <- imputed_trees[[sample(1:length(imputed_trees),1)]]
        Matnull$phylo_richness[n] <- pd(cell_comm, temp_tree, include.root=TRUE)[1,'PD']
        Matnull$phylo_diversity[n] <- mpd(cell_comm, cophenetic(temp_tree), abundance.weighted=FALSE)
        
      }
      
    }
    
    # save
    reg_pool[[scl_vrich[r]]] <- Matnull
    
    # progress
    print(rep('+',scl_vrich[r]))
    
  }
  
  regional_pool_list[[reg1$code]] <- reg_pool

  # progress
  print(paste('--- ', round(e/nrow(regional_pools)*100, 2), '% ---', sep=''))

}

save(regional_pool_list, file='results/regional_pool_list.Rdata')


