

# NULL MODELS FOR ALPHA DIVERSITY INDICES (EXCEPT RICHNESS)


library(tidyverse)
library(cluster)
library(BAT)
library(FD) # Fdisp
library(randtip) # imputation
library(picante) # pd (Faith), mpd


# load data for analyses
source('scripts/import data alpha and beta diversity calc.R')


# Xnull randomization: shuffling the taxon labels on the phylogeny & func.distance matrix + species richness constant

# we are going to define two pools:
  # global pool: any species can occur in any grid (Swenson et al. 2012)
  # regional pool: species ranges are constant (Spalink et al. 2018)



# global pool: the most efficient way to do this is to generate a distribution for each species richness value ####
scl_vrich <- read.csv("results/slc_indices.txt", sep="") %>% dplyr::select(taxon_richness) %>% deframe() %>% unique()
scl_vrich <- scl_vrich[-which(scl_vrich==1)] # at least two species

global_pool_list <- list()

# site x species presence matrix to subset communities from
scl_presences <- matrix(ncol=nrow(scl_traits), nrow=1)
colnames(scl_presences) <- rownames(scl_traits)
scl_presences[1,] <- 1

Xnull=100

Matnull <- matrix(ncol=4, nrow=Xnull) %>% data.frame() # results for every combination
colnames(Matnull) <- c('Frich_','Fdisp_','Faith_','mpd_')

for (r in 1:length(scl_vrich)) { # for every observed richness value
  
  Matnull[] <- NA
  
  for (n in 1:nrow(Matnull)) {  # create X random values
    
    cell_spp <- sample(rownames(scl_traits), scl_vrich[r], replace=FALSE)
    
    # presence matrix
    cell_comm <- t(as.matrix(scl_presences[1,cell_spp]))
    
    # Frich & Fdisp
    Matnull$Frich_[n] <- alpha(cell_comm, scl_dendrogram)
    tempFD <- as.matrix(scl_dist)[cell_spp,cell_spp] %>% as.dist() %>% dbFD(stand.FRic=F, messages=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
    Matnull$Fdisp_[n] <- tempFD$FDis
    
    # calculate average Faith and mpd from imputed trees
    temp_tree <- imputed_trees[[sample(1:length(imputed_trees),1)]]
    Matnull$Faith_[n] <- pd(cell_comm, temp_tree, include.root=TRUE)[1,'PD']
    Matnull$mpd_[n] <- mpd(cell_comm, cophenetic(temp_tree), abundance.weighted=FALSE)
    
  }
  
  # save
  global_pool_list[[as.character(scl_vrich[r])]] <- Matnull
  
  # progress
  print(paste('--- ', round(r/length(scl_vrich)*100, 2), '% ---', sep=''))
  
}

save(global_pool_list, file='results/global_pool_list.Rdata')



# regional pools: identify species that occur within the same BIOME in the same REALM ####
ecoregions <- read.csv("results/df_ecoregions.txt", sep="") %>%
  subset(scientific_name%in%v_spp) %>%
  dplyr::select(scientific_name,REALMxBIOME)

regional_pools <- unique(ecoregions[,'REALMxBIOME']) %>% as.data.frame() # 39 potential pools, 17 species per pool on average
colnames(regional_pools) <- 'REALMxBIOME'
regional_pools$n_species <- NA
for (p in 1:nrow(regional_pools)) {
  regional_pools$n_species[p] <- ecoregions %>%
    subset(REALMxBIOME==regional_pools$REALMxBIOME[p]) %>%
    dplyr::select(scientific_name) %>% deframe() %>% length()
  }
regional_pools <- regional_pools[regional_pools$n_species>2,] # remove combinations with 2 or less species
write.table(regional_pools, 'results/regional_pools.txt')

regional_pool_list <- list()

Matnull <- matrix(ncol=4, nrow=Xnull) %>% data.frame() # results for every combination
colnames(Matnull) <- c('Frich_','Fdisp_','Faith_','mpd_')

for (e in 1:nrow(regional_pools)) {
  
  reg1 <- regional_pools[e,]
  pool1 <- ecoregions %>% subset(REALMxBIOME==reg1$REALMxBIOME)
  pool1 <- pool1$scientific_name
  
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
        Matnull$Frich_[n] <- tempFD$FRic
        Matnull$Fdisp_[n] <- tempFD$FDis
        
        # calculate average Faith and mpd from imputed trees
        temp_tree <- imputed_trees[[sample(1:length(imputed_trees),1)]]
        Matnull$Faith_[n] <- pd(cell_comm, temp_tree, include.root=TRUE)[1,'PD']
        Matnull$mpd_[n] <- mpd(cell_comm, cophenetic(temp_tree), abundance.weighted=FALSE)
        
      }
      
    }
    
    # save
    reg_pool[[as.character(scl_vrich[r])]] <- Matnull
    
    # progress
    print(rep('+',scl_vrich[r]))
    
  }
  
  regional_pool_list[[reg1$REALMxBIOME]] <- reg_pool

  # progress
  print(paste('--- ', round(e/nrow(regional_pools)*100, 2), '% ---', sep=''))

}

save(regional_pool_list, file='results/regional_pool_list.Rdata')


