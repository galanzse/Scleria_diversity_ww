


# CALCULATE ALPHA TAXONOMIC, FUNCTIONAL AND PHYLOGENETIC DIVERSITY, AND SES


library(tidyverse)
library(rnaturalearth)
library(cluster)
library(terra)
library(FD) # Fdisp
library(BAT)
library(picante) # pd (Faith), mpd
library(SpatialPack) # modified.ttest


# load data for analyses
source('scripts/import data alpha and beta diversity calc.R')


# null models and BIOMExREALM id
# regional_pools <- read.csv("results/regional_pools.txt", sep="")
# rast_ecoregions <- rast('results/rast_ecoregions.tiff')
load("results/global_pool_list.Rdata")
# load("results/regional_pool_list.Rdata")


# compute alpha diversity and trait means ####

# unique cells
scl_occurrences$cell <- terra::extract(world_grid, scl_points, cells=TRUE) %>% dplyr::select(cell) %>% deframe()
# scl_occurrences$BIOMExREALM <- terra::extract(rast_ecoregions, scl_points, cells=TRUE) %>% # info regional pools
#   dplyr::select(REALMxBIOME) %>% deframe()
scl_indices <- scl_occurrences %>% dplyr::select(cell) %>% unique() # BIOMExREALM

# table(is.na(scl_indices$BIOMExREALM)) # NAs in regional pool
# scl_indices$BIOMExREALM[which(!(scl_indices$BIOMExREALM %in% ecoregions$REALMxBIOME))] %>% unique()
# scl_indices$BIOMExREALM[scl_indices$BIOMExREALM %in% c('AT98','AA14','NA2','PA13','NA98')] <- NA

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

# scl_indices$SES_Frich_rg <- NA #    SES REGIONAL
# scl_indices$SES_Fdisp_rg <- NA #
# scl_indices$SES_Faith_rg <- NA #
# scl_indices$SES_mpd_rg <- NA #


# loop
for (cl in 1:nrow(scl_indices)) {
  
  cell_spp <- scl_occurrences %>% subset(cell==scl_indices$cell[cl]) %>%
    dplyr::select(scientific_name) %>% unique() %>% deframe()
  
  # remove species not present in trait/phylogeny matrices
  cell_spp <- cell_spp[which(cell_spp%in%colnames(scl_presences))]
  
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


  if (length(cell_spp)>1) {
    scl_indices$Frich_[cl] <- alpha(cell_comm, scl_dendrogram) # Frich
    tempFD <- as.matrix(scl_dist)[cell_spp,cell_spp] %>% as.dist() %>% dbFD(stand.FRic=F, messages=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
    scl_indices$Fdisp_[cl] <- tempFD$FDis # Fdisp
  }
  
  # # calculate average Faith and mpd from imputed trees
  # # Xi=length(imputed_trees)
  Xi=1000
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
    # if (!(is.na(scl_indices$BIOMExREALM[cl]))) {
    #   pool2 <- regional_pool_list[[as.character(scl_indices$BIOMExREALM[cl])]]
    #   if (!(is.null(pool2))) { # some rare pools are absent
    #     pool2 <- pool2[[as.character(scl_indices$richness[cl])]]
    #     percentile <- ecdf(pool2$Frich_); scl_indices$SES_Frich_rg[cl] <- percentile(scl_indices$Frich_[cl])
    #     percentile <- ecdf(pool2$Fdisp_); scl_indices$SES_Fdisp_rg[cl] <- percentile(scl_indices$Fdisp_[cl])
    #     percentile <- ecdf(pool2$Faith_); scl_indices$SES_Faith_rg[cl] <- percentile(scl_indices$Faith_[cl])
    #     percentile <- ecdf(pool2$mpd_); scl_indices$SES_mpd_rg[cl] <- percentile(scl_indices$mpd_[cl])
    #   }
    # }
  }
  
  

  # progress
  print(paste('--- ', round(cl/nrow(scl_indices)*100, 2), '% ---', sep=''))
  
}


# scl_indices[is.na(scl_indices)] <- 0
# write.table(scl_indices, 'results/scl_indices.txt')



# maps ####


# templates to fill
map_richness <- world_grid; names(map_richness) <- 'richness'

map_ann_per <- world_grid; names(map_ann_per) <- 'ann_per'
map_cwm_height <- world_grid; names(map_cwm_height) <- 'CWM_height'
map_cwm_blade <- world_grid; names(map_cwm_blade) <- 'CWM_bladearea'
map_cwm_nutlet <- world_grid; names(map_cwm_nutlet) <- 'CWM_nutletvolume'

map_Frich <- world_grid; names(map_Frich) <- 'Frich' #   OBSERVED
map_Fdisp <- world_grid; names(map_Fdisp) <- 'Fdisp'
map_Faith <- world_grid; names(map_Faith) <- 'Faith'
map_mpd <- world_grid; names(map_mpd) <- 'mpd'

map_Frich_SESgl <- world_grid; names(map_Frich_SESgl) <- 'SES_Frich_gl' #   SES GLOBAL
map_Fdisp_SESgl <- world_grid; names(map_Fdisp_SESgl) <- 'SES_Fdisp_gl'
map_Faith_SESgl <- world_grid; names(map_Faith_SESgl) <- 'SES_Faith_gl'
map_mpd_SESgl <- world_grid; names(map_mpd_SESgl) <- 'SES_mpd_gl'

# map_Frich_SESrg <- world_grid; names(map_Frich_SESrg) <- 'SES_Frich_rg' #   SES REGIONAL
# map_Fdisp_SESrg <- world_grid; names(map_Fdisp_SESrg) <- 'SES_Fdisp_rg'
# map_Faith_SESrg <- world_grid; names(map_Faith_SESrg) <- 'SES_Faith_rg'
# map_mpd_SESrg <- world_grid; names(map_mpd_SESrg) <- 'SES_mpd_rg'



# loop to fill values
for (cl in 1:nrow(scl_indices)) {
  
  map_richness[scl_indices$cell[cl]] <- scl_indices$richness[cl]
  
  map_ann_per[scl_indices$cell[cl]] <- scl_indices$ann_per[cl]
  map_cwm_height[scl_indices$cell[cl]] <- scl_indices$cwm_height[cl]
  map_cwm_blade[scl_indices$cell[cl]] <- scl_indices$cwm_blade[cl]
  map_cwm_nutlet[scl_indices$cell[cl]] <- scl_indices$cwm_nutlet[cl]
  
  map_Frich[scl_indices$cell[cl]] <- scl_indices$Frich_[cl]
  map_Fdisp[scl_indices$cell[cl]] <- scl_indices$Fdisp_[cl]
  map_Faith[scl_indices$cell[cl]] <- scl_indices$Faith_[cl]
  map_mpd[scl_indices$cell[cl]] <- scl_indices$mpd_[cl]
  
  map_Frich_SESgl[scl_indices$cell[cl]] <- scl_indices$SES_Frich_gl[cl]
  map_Fdisp_SESgl[scl_indices$cell[cl]] <- scl_indices$SES_Fdisp_gl[cl]
  map_Faith_SESgl[scl_indices$cell[cl]] <- scl_indices$SES_Faith_gl[cl]
  map_mpd_SESgl[scl_indices$cell[cl]] <- scl_indices$SES_mpd_gl[cl]
  
  # map_Frich_SESrg[scl_indices$cell[cl]] <- scl_indices$SES_Frich_rg[cl]
  # map_Fdisp_SESrg[scl_indices$cell[cl]] <- scl_indices$SES_Fdisp_rg[cl]
  # map_Faith_SESrg[scl_indices$cell[cl]] <- scl_indices$SES_Faith_rg[cl]
  # map_mpd_SESrg[scl_indices$cell[cl]] <- scl_indices$SES_mpd_rg[cl]
  
  
  # progress
  print(paste('--- ', round(cl/nrow(scl_indices)*100, 2), '% ---', sep=''))
  
}


# save rasters
# alpha_div_rasters <- c(map_richness, map_ann_per, map_cwm_height, map_cwm_blade, map_cwm_nutlet, map_Frich, map_Fdisp, map_Faith, map_mpd, map_Frich_SESgl, map_Fdisp_SESgl, map_Faith_SESgl, map_mpd_SESgl) # , map_Frich_SESrg, map_Fdisp_SESrg, map_Faith_SESrg, map_mpd_SESrg
# writeRaster(alpha_div_rasters, 'results/alpha_div_rasters.tiff', overwrite=TRUE)



# prepare maps for representation

alpha_div_rasters[['richness']] <- log(alpha_div_rasters[['richness']])
names(alpha_div_rasters)[1] <- 'log(richness)'

# matrix to reclassify SES maps
rc_mat <- matrix(nrow=length(seq(0, 1, 0.01)), ncol=2)
colnames(rc_mat) <- c('from','to')
rc_mat[,'from'] <- seq(0, 1, 0.01)
rc_mat[1:6,'to'] <- 1
rc_mat[7:95,'to'] <- 2
rc_mat[96:101,'to'] <- 3

v_SES <- c("SES_Frich_gl","SES_Fdisp_gl","SES_Faith_gl","SES_mpd_gl","SES_Frich_rg","SES_Fdisp_rg","SES_Faith_rg","SES_mpd_rg")

# plot
for (s in 1:dim(alpha_div_rasters)[3]) {

    pdf(file = paste("results/alpha diversity_150/",names(alpha_div_rasters)[s],".pdf",sep=""),
      width = 9.30, # The width of the plot in inches
      height = 5.74) # The height of the plot in inches
  
  alpha_div_rasters[[s]] %>%
    plot(main=names(alpha_div_rasters)[s], col=rev(grDevices::heat.colors(20)[1:15]))
  lines(world_lines)
  
  dev.off()
  
}



# alpha diversity ~ climate ####

# response variables 
alphadiv <- rast('results/alpha_div_rasters.tiff')[[1:9]]
names(alphadiv)

# bioclim
wc_bioclim <- list.files('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m', full.names=T) %>% rast()
# [-which(list.files('C:/Users/user/Desktop/wc2.1_30s_bio', full.names=T)=="C:/Users/user/Desktop/wc2.1_30s_bio/variables.txt")]
names(wc_bioclim) <- substr(names(wc_bioclim), 12, 17) #11, 16

rs_wc_bioclim <- wc_bioclim[[names(wc_bioclim)%in%c('bio_1','bio_6','bio_7','bio_12','bio_13','bio_15')]] %>%
  aggregate(fact=15, fun="mean") %>% # aggregate to reduce computation time
  project('+proj=eqearth') %>% # project
  resample(y=alphadiv[[1]], method='average') # resample

alphadiv_bioclim <- c(alphadiv, rs_wc_bioclim) # group rasters

df_alphadiv_bioclim <- alphadiv_bioclim %>% as.data.frame(xy=TRUE) %>% # convert into dataframe
  na.omit() %>% subset(CWM_height > 0) # remove NAs and zeros

# par(mfrow=c(1,1), mar=c(5,5,5,5))
# plot(CWM_height ~ bio_1, data=df_CWM_bioclim); abline(lm(CWM_height ~ bio_1, data=df_CWM_bioclim), col='blue', lwd=3)

test_alphadiv_bioclim <- matrix(ncol=6, nrow=dim(alphadiv)[3]*6) %>% as.data.frame()
colnames(test_alphadiv_bioclim) <- c('ALPHA', 'BIOCLIM', 'F_', 'DF', 'R', 'P_value')
test_alphadiv_bioclim$ALPHA <- rep(names(alphadiv), 6)
test_alphadiv_bioclim$BIOCLIM <- c( rep('bio_1',9), rep('bio_6',9), rep('bio_7',9),
                                    rep('bio_12',9), rep('bio_13',9), rep('bio_15',9))
                                    
for (i in 1:nrow(test_alphadiv_bioclim)) {
  
  modtest1 <- modified.ttest(df_alphadiv_bioclim[,test_alphadiv_bioclim$ALPHA[i]],
                             df_alphadiv_bioclim[,test_alphadiv_bioclim$BIOCLIM[i]],
                             df_alphadiv_bioclim[,c('x','y')], nclass=13)
  
  test_alphadiv_bioclim$F_[i] <- modtest1$Fstat
  test_alphadiv_bioclim$DF[i] <- modtest1$dof
  test_alphadiv_bioclim$R[i] <- modtest1$corr
  test_alphadiv_bioclim$P_value[i] <- modtest1$p.value
  
  print(paste('--- ', round(i/nrow(test_alphadiv_bioclim),2)*100,'% ---', sep=''))
  
}

write.csv(test_alphadiv_bioclim, 'results/test_alphadiv_bioclim.csv')


