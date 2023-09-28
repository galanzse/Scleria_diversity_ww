

# IS NUTLET SIZE CORRELATED TO DISPERSAL CAPACITY AT THE SPECIES LEVEL?


library(readxl)
library(tidyverse)
library(ggpubr)
library(terra)



# import data ####

load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")

# traits
scl_nutlets <- data_final[['traits']] %>% dplyr::select(scientific_name, subgenus, section, nutlet_volume, blade_area, height) %>% na.omit()

# automated assessments that include AOO and EOO
IUCN_results <- read_excel("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/IUCN_results.xlsx")

# ecorregions
scl_ecorregions <- read.table('results/df_ecoregions.txt')



# Predictors: number of realms (aprox. continents) where it is present, Biomes, EEO, AOO, locations, latitudinal and longitudinal range
scl_nutlets <- merge(scl_nutlets, IUCN_results[,c('scientific_name', 'Category_code', 'EOO', 'AOO', 'Nbe_loc')], by='scientific_name', all.x=T)

scl_nutlets$x_range <- NA
scl_nutlets$y_range <- NA
scl_nutlets$n_realm <- NA
scl_nutlets$biome <- NA

for (i in 1:nrow(scl_nutlets)) {

  temp_occ <- data_final[['occurrences']] %>% subset(scientific_name==scl_nutlets$scientific_name[i]) %>%
    dplyr::select(x,y)
  
  scl_nutlets$x_range[i] <-  max(temp_occ$x) - min(temp_occ$x)
    
  scl_nutlets$y_range[i] <- max(temp_occ$y) - min(temp_occ$y)

  scl_nutlets$n_realm[i] <- scl_ecorregions %>% subset(scientific_name==scl_nutlets$scientific_name[i]) %>%
    dplyr::select(REALM) %>% unique() %>% na.omit() %>% nrow()
  
  scl_nutlets$biome[i] <- scl_ecorregions %>% subset(scientific_name==scl_nutlets$scientific_name[i]) %>%
    dplyr::select(BIOME) %>% unique() %>% na.omit() %>% nrow()
  
  # 
  print(scl_nutlets$scientific_name[i])
  
}


# categorise nutlets
my_comparisons <- list( c("Small", "Medium"), c("Small", "Big"), c("Medium", "Big"))
qnt1 <- quantile(scl_nutlets$nutlet_volume,c(0,1/3,2/3,1), na.rm=T)[2:3]
scl_nutlets$nutlet_cat <- cut(scl_nutlets$nutlet_volume, breaks=c(-Inf, qnt1[1], qnt1[2], Inf), labels=c("Small","Medium","Big"))


# save results
write.table(scl_nutlets, 'results/scl_nutlets.txt')



# analyses ####

scl_nutlets <- read.csv("results/scl_nutlets.txt", sep="")
str(scl_nutlets)
scl_nutlets$nutlet_cat <- as.factor(scl_nutlets$nutlet_cat)
scl_nutlets$nutlet_cat <- factor(scl_nutlets$nutlet_cat, c('Small','Medium','Big'))



# my_comparisons <- list( c("Small", "Medium"), c("Small", "Big"), c("Medium", "Big"))
ggplot(aes(y=log(EOO), x=nutlet_cat), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('log(EOO)') + xlab('Nutlet size') +
  theme_classic() 
  # stat_compare_means(method='t.test', paired=F, comparisons=my_comparisons, bracket.size=.1,
  #                    label="p.signif", hide.ns=T, vjust=0.4)

ggplot(aes(y=log(EOO), x=subgenus), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('log(EOO)') + xlab('') +
  theme_classic() 



anova(lm(log(EOO) ~ nutlet_cat, data=scl_nutlets))
ggplot(aes(x=nutlet_cat, y=log(EOO)), data=scl_nutlets) + # AOO
  geom_boxplot() +
  geom_smooth(method='lm') +
  xlab('height') + ylab('log(EOO)') +
  theme_classic()

ggplot(aes(y=y_range, x=nutlet_cat), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('Latitudinal range') + xlab('Nutlet size category') +
  theme_classic()

ggplot(aes(y=x_range, x=nutlet_cat), data=scl_nutlets) + # AOO
  geom_boxplot() +
  ylab('Longitudinal range') + xlab('Nutlet size category') +
  theme_classic()


