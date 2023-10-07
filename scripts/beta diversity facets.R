

# CALCULATE BETA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY


library(tidyverse)
library(terra)
library(vegan)
library(BAT)
library(ade4)


# load data for analyses
source('scripts/import data alpha and beta diversity calc.R')



# Presence matrix
scl_comm_mat <- matrix(nrow=length(unique(scl_occurrences$cell)), ncol=length(unique(scl_occurrences$scientific_name)))
rownames(scl_comm_mat) <- unique(scl_occurrences$cell)
colnames(scl_comm_mat) <- unique(scl_occurrences$scientific_name)

for (i in 1:nrow(scl_comm_mat)) { # presence matrix
  
  spp1 <- scl_occurrences %>% filter(cell==rownames(scl_comm_mat)[i]) %>%
    dplyr::select(scientific_name) %>% unique() %>% deframe()
  scl_comm_mat[i,colnames(scl_comm_mat)%in%spp1] <- 1
  scl_comm_mat[i,!(colnames(scl_comm_mat)%in%spp1)] <- 0
  
  print(paste('--- ', round(i/nrow(scl_comm_mat), 2)*100, '% ---', sep='')) 
  
}



# Taxonomic beta diversity: Bray-Curtis
scl_beta_tax <- vegdist(scl_comm_mat, method="bray", binary=FALSE, diag=FALSE, upper=TRUE)

# Functional beta diversity: Btotal (total) = Brepl (replacement) + Brich (species loss/gain)
scl_beta_fun <- BAT::beta(scl_comm_mat, scl_dendrogram, func="jaccard", abund=F, raref=0, comp=F)

# Phylogenetic beta diversity
scl_beta_phy <- BAT::beta(scl_comm_mat, scl_phylogeny, func="jaccard", abund=F, raref=0, comp=F)

# save results
scl_beta_div <- list()
scl_beta_div[['taxonomic']] <- scl_beta_tax
scl_beta_div[['functional']] <- scl_beta_fun
scl_beta_div[['phylogenetic']] <- scl_beta_phy

save(scl_beta_div, file='results/scl_beta_div.Rdata')



# correlation between matrices
mantel.rtest


