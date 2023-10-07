

library(tidyverse)
library(ggfortify)
library(ggplot2)
library(hypervolume)
library(terra)


# load data for analyses
load("results/data_final.RData")


# bioclim
wc_bioclim <- list.files('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m', full.names=T) %>% rast()
# [-which(list.files('C:/Users/user/Desktop/wc2.1_30s_bio', full.names=T)=="C:/Users/user/Desktop/wc2.1_30s_bio/variables.txt")] %>% rast()
names(wc_bioclim) <- substr(names(wc_bioclim), 12, 17)
wc_bioclim <- wc_bioclim[[names(wc_bioclim)%in%c('bio_1','bio_6','bio_7','bio_12','bio_13','bio_15')]]

# PCA to explore climatic niches
scl_climatic <- wc_bioclim %>% extract(vect(data_final[['occurrences']], geom=c('x','y')), ID=F)
scl_climatic$scientific_name <- data_final[['occurrences']]$scientific_name
scl_climatic <- merge(scl_climatic, data_final[['taxa']][,c('subgenus','section','scientific_name')])
scl_climatic <- na.omit(scl_climatic)

clim_pca <- prcomp(scl_climatic[,c("bio_1","bio_6","bio_7","bio_12","bio_13","bio_15")], scale=T, center=T)
summary(clim_pca)

# autoplot(clim_pca, data=scl_climatic, colour='section', loadings=TRUE, loadings.label=TRUE)
autoplot(clim_pca, data=scl_climatic, colour='subgenus', loadings=TRUE, loadings.label=TRUE) +
  scale_color_manual(values = c('chocolate1','palegreen', 'plum1','lightgoldenrod2')) +
  theme_bw()

# tabla para analysis
scl_climatic <- cbind(scl_climatic, clim_pca$x[,1:3])



# SUBGENERA ####

# determino las caracteristicas del bucle
N_analyses <- 100
table(scl_climatic$subgenus)
N_observations <- 300

# matriz de resultados
hyp_volume <- matrix(nrow=N_analyses, ncol=4)
colnames(hyp_volume) <- c('Browniae', 'Hypoporum', 'Trachylomia', 'Scleria')
hyp_overlap <- matrix(nrow = N_analyses, ncol = 6)
colnames(hyp_overlap) <- c('Brow.Hypo', 'Brow.Trac', 'Brow.Scle', 'Hypo.Trac', 'Hypo.Scle', 'Trac.Scle')

for (i in 1:N_analyses) {
  
  mysample <- filter(scl_climatic, subgenus=='Browniae')
  mysample <- mysample[sample(nrow(mysample), N_observations), c('PC1','PC2','PC3')]
  hyp_browniae<-hypervolume(mysample, method='box')
  hyp_volume[i,'Browniae'] <- hyp_browniae@Volume
  
  mysample <- filter(scl_climatic, subgenus=='Hypoporum')
  mysample <- mysample[sample(nrow(mysample), N_observations), c('PC1','PC2','PC3')]
  hyp_hypoporum<-hypervolume(mysample, method='box')
  hyp_volume[i,'Hypoporum'] <- hyp_hypoporum@Volume
  
  mysample <- filter(scl_climatic, subgenus=='Trachylomia')
  mysample <- mysample[sample(nrow(mysample), N_observations), c('PC1','PC2','PC3')]
  hyp_trachylomia<-hypervolume(mysample, method='box')
  hyp_volume[i,'Trachylomia'] <- hyp_trachylomia@Volume
  
  mysample <- filter(scl_climatic, subgenus=='Scleria')
  mysample <- mysample[sample(nrow(mysample), N_observations), c('PC1','PC2','PC3')]
  hyp_scleria<-hypervolume(mysample, method='box')
  hyp_volume[i,'Scleria'] <- hyp_scleria@Volume
  
  
  # overlaps sorensen
  hyp <- hypervolume_set(hyp_browniae, hyp_hypoporum, check.memory = FALSE)
  hyp_overlap[i,'Brow.Hypo'] <- hypervolume_overlap_statistics(hyp)[2]
  
  hyp <- hypervolume_set(hyp_browniae, hyp_trachylomia, check.memory = FALSE)
  hyp_overlap[i,'Brow.Trac'] <- hypervolume_overlap_statistics(hyp)[2]
  
  hyp <- hypervolume_set(hyp_browniae, hyp_scleria, check.memory = FALSE)
  hyp_overlap[i,'Brow.Scle'] <- hypervolume_overlap_statistics(hyp)[2]
  
  hyp <- hypervolume_set(hyp_hypoporum, hyp_trachylomia, check.memory = FALSE)
  hyp_overlap[i,'Hypo.Trac'] <- hypervolume_overlap_statistics(hyp)[2]
  
  hyp <- hypervolume_set(hyp_hypoporum, hyp_scleria, check.memory = FALSE)
  hyp_overlap[i,'Hypo.Scle'] <- hypervolume_overlap_statistics(hyp)[2]
  
  hyp <- hypervolume_set(hyp_trachylomia, hyp_scleria, check.memory = FALSE)
  hyp_overlap[i,'Trac.Scle'] <- hypervolume_overlap_statistics(hyp)[2]
  
}

rm(mysample, i, N_analyses, N_observations, hyp, hyp_browniae, hyp_hypoporum, hyp_trachylomia, hyp_scleria)

write.table(hyp_volume, 'results/hyp_subgenera_vol.txt')
write.table(hyp_overlap, 'results/hyp_subgenera_ove.txt')


# riqueza del nicho climatico
long_hyp_volume <- hyp_volume %>% as.data.frame() %>% pivot_longer(1:4, names_to='subgenus', values_to='hypervolume')
ggplot(aes(x=subgenus, y=hypervolume, fill=subgenus), data=long_hyp_volume) +
  geom_boxplot()+
  scale_fill_manual(values=c('chocolate1','palegreen', 'plum1','lightgoldenrod2')) +
  theme_bw() +
  theme(legend.position='none') + xlab('')

# overlap entre nichos climaticos
long_hyp_overlap <- hyp_overlap %>% as.data.frame() %>% pivot_longer(1:6, names_to='subgenus', values_to='overlap')
ggplot(aes(x=subgenus, y=overlap, fill=subgenus), data=long_hyp_overlap) +
  geom_boxplot()+
  theme_bw() +
  theme(legend.position='none') + xlab('')



# SECTION ####

N_analyses <- 100
table(datapca$section)
N_observations <- 45

# eliminamos las secciones con menos de 5. observaciones
sections50 <- scl_climatic$section %>% table() %>% as.data.frame()
colnames(sections50) <- c('section','Freq')

# matriz de resultados
hyp_volume_sec <- matrix(nrow = N_analyses, ncol = length(unique(sections50$section)))
colnames(hyp_volume_sec) <- unique(sections50$section)

for (sect in unique(sections50$section)) {
  
  ss_sp <- filter(scl_climatic, section==sect)
  
  for (i in 1:N_analyses) {
    mysample <- ss_sp[sample(nrow(ss_sp), N_observations), c('PC1','PC2','PC3')]
    hyp_sect <- hypervolume(mysample, method="box")
    # hyp_sect <- hypervolume_box(mysample, kde.bandwidth=0.5)
    hyp_volume_sec[i,sect] <- hyp_sect@Volume
  }
  
  rm(hyp_sect, i, sect)
}

rm(N_analyses, N_observations, mysample, ss_sp)

write.table(hyp_volume_sec, 'results/hyp_section_vol.txt')


# riqueza del nicho climatico
long_hyp_volume_sec <- hyp_volume_sec %>% as.data.frame() %>% pivot_longer(1:17, names_to='section', values_to='hypervolume')
ggplot(aes(x=section, y=hypervolume, fill=section), data=long_hyp_volume_sec) +
  geom_boxplot()+
  theme_bw() + xlab('') +
  theme(legend.position='none', axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))



# hay relacion entre la riqueza y la diversidad?

sppxsect <- scl_climatic %>% select(section) %>% unique()
sppxsect$hypervolume <- NA
sppxsect$richness <- NA
for (i in 1:nrow(sppxsect)) {
  sppxsect$hypervolume[i] <- hyp_volume_sec[,sppxsect$section[i]] %>% mean()
  sppxsect$richness[i] <- scl_climatic %>% subset(section==sppxsect$section[i]) %>% select(scientific_name) %>% unique() %>% nrow()
}

plot(sppxsect$hypervolume ~ log(sppxsect$richness)); abline(lm(sppxsect$hypervolume ~ log(sppxsect$richness)))



