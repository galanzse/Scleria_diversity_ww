

# EXPLORE RELATIONSHIP BETWEEN FUNCTIONALITY AND CLIMATIC NICHE IN SCLERIA AT THE SUBGENUS LEVEL


library(tidyverse); library(ggpubr)
library(terra)
library(randomForest)
library(caret); library(caTools)
library(SpatialPack) # modified.ttest

load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData")



# functional differences ####

LHS_data <- data_final[['traits']] %>% dplyr::select(subgenus, section, height, blade_area, nutlet_volume)
v_traits <- c('height','blade_area','nutlet_volume')

par(mfrow=c(1,3)); hist(LHS_data$height); hist(LHS_data$blade_area); hist(LHS_data$nutlet_volume)
LHS_data[,v_traits] <- LHS_data[,v_traits] %>% apply(2, log)
long_LHS_data <- LHS_data %>% pivot_longer(3:5, names_to='trait', values_to='value') %>% na.omit()

long_LHS_data$subgenus <- str_to_title(long_LHS_data$subgenus) # correct names
long_LHS_data$subgenus <- as.factor(long_LHS_data$subgenus)

# there is a clear functional segregation at the subgenus level for leaf and nutlet traits
lbs = setNames(c("'log10 blade area ('*cm^2*')'", 
                 "'log10 maximum height (cm)'", 
                 "'log10 nutlet volume ('*mm^3*')'"),
               c('blade_area','height','nutlet_volume'))

my_comparisons <- list( c("Scleria", "Browniae"), c("Trachylomia", "Browniae"), c("Hypoporum", "Browniae"),
                        c("Hypoporum", "Trachylomia"), c("Scleria", "Hypoporum"),
                        c("Scleria", "Trachylomia"))

ggplot(aes(x=subgenus, y=value), data=long_LHS_data) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) +
  facet_wrap(.~trait, scales="free", ncol=3, labeller=as_labeller(lbs, label_parsed)) +
  theme_classic()
  # theme(axis.text.x=element_text(angle=-45, vjust=1.2, hjust=0)) +
  # stat_compare_means(method='t.test', paired=F, comparisons=my_comparisons, bracket.size=.1,
  #                    label="p.signif", hide.ns=T, vjust=0.4)



# climatic differences ####

# wc_bioclim <- list.files('C:/Users/user/Desktop/wc2.1_30s_bio', full.names=T)[-which(list.files('C:/Users/user/Desktop/wc2.1_30s_bio', full.names=T)=="C:/Users/user/Desktop/wc2.1_30s_bio/variables.txt")] %>% rast()
# names(wc_bioclim) <- substr(names(wc_bioclim), 11, 16)
# df_bioclim <- df_bioclim %>% terra::extract(data_final[['occurrences']][,c('x','y')], ID=F)
# scl_bioclim <- cbind(data_final[['occurrences']][,'scientific_name'], df_bioclim)
# colnames(scl_bioclim)[1] <- 'scientific_name'
# scl_bioclim <- merge(data_final[['traits']][,c('scientific_name','section','subgenus')], scl_bioclim)
# write.table(scl_bioclim, 'results/scl_bioclim.txt')


# We will use the following variables as they represent average, extreme and seasonal indexes of temperature and precipitation:
# lbs = c('BIO1 = Annual Mean Temperature', 'BIO6 = Min Temperature of Coldest Month', 'BIO7 = Temperature Annual Range (BIO5-BIO6)', 'BIO12 = Annual Precipitation', 'BIO13 = Precipitation of Wettest Month', 'BIO15 = Precipitation Seasonality (Coefficient of Variation)')

scl_bioclim <- scl_bioclim %>% dplyr::select(subgenus, bio_1, bio_6, bio_7, bio_12, bio_13, bio_15)

# Do Scleria Subgenera have different average climatic niches?
long_scl_bioclim <- scl_bioclim %>% pivot_longer(2:7, names_to='bioclim', values_to='value') %>% na.omit()
long_scl_bioclim$bioclim <- as.factor(long_scl_bioclim$bioclim)
long_scl_bioclim$bioclim <- factor(long_scl_bioclim$bioclim, levels=c('bio_1','bio_6', 'bio_7', 'bio_12', 'bio_13', 'bio_15'))

long_scl_bioclim$subgenus <- as.factor(long_scl_bioclim$subgenus)

lbs = setNames(c("'Annual Mean Temperature ' (degree*C)",
                 "'Min Temperature of Coldest Month ' (degree*C)",
                 "'Temperature Annual Range ' (degree*C)",
                 "'Annual Precipitation (mm)'",
                 "'Precipitation of Wettest Month (mm)'",
                 "'Precipitation Seasonality (CV)'"),
               levels(long_scl_bioclim$bioclim))

ggplot(aes(x=subgenus, y=value), data=long_scl_bioclim) +
  geom_boxplot() + xlab(NULL) + ylab(NULL) +
  facet_wrap(.~bioclim, scales="free", nrow=2, labeller=as_labeller(lbs, label_parsed)) +
  theme_classic()



# random forest ####

# Can we predict species' infrageneric taxa based on LHS scheme and climatic niches?
LHS_data <- data_final[['traits']] %>% dplyr::select(scientific_name, subgenus, section, height, blade_area, nutlet_volume)
scl_bioclim <- read.csv("results/scl_bioclim.txt", sep="") %>% dplyr::select(scientific_name, bio_1, bio_6, bio_7, bio_12, bio_13, bio_15)
scl_bioclim <- merge(data_final[['traits']][,c('scientific_name','section','subgenus')], scl_bioclim)

scl_bioclim <- scl_bioclim %>% group_by(scientific_name,section,subgenus) %>%
  summarise(bio_1=mean(bio_1, na.rm=T),
            bio_6=mean(bio_6, na.rm=T),
            bio_7=mean(bio_7, na.rm=T),
            bio_12=mean(bio_12, na.rm=T),
            bio_13=mean(bio_13, na.rm=T),
            bio_15=mean(bio_15, na.rm=T))

rf_data <- merge(LHS_data, scl_bioclim, by=c('scientific_name','section','subgenus')) %>% na.omit()

# 20% of observations for external validation
v_ex <- sample(1:nrow(rf_data), nrow(rf_data)/5)
data_ex <- rf_data[v_ex,]; table(data_ex$subgenus)
data_tr <- rf_data[-v_ex,]

v_vars <- c('subgenus', colnames(rf_data)[4:12])

# 5-fold cross validation 
control <- trainControl(method="repeatedcv", number=5, repeats=5, search="grid")
rf_random <- train(subgenus~., data=data_tr[,v_vars], method="rf", metric="Accuracy", trControl=control)
rf_random
plot(rf_random)
varImp(rf_random, scale=F)
plot(varImp(rf_random, scale=F))

# confusion matrix
predicted_class <- predict(rf_random, newdata=data_ex[,v_vars], type='raw')
confusionMatrix(data=predicted_class, reference=as.factor(data_ex$subgenus))



# CWM ~ climate ####

# CWM 
CWM <- rast('results/alpha_div_rasters.tiff')[[c('CWM_height', 'CWM_bladearea', 'CWM_nutletvolume')]]

# bioclim
rs_wc_bioclim <- wc_bioclim[[names(wc_bioclim)%in%c('bio_1','bio_6','bio_7','bio_12','bio_13','bio_15')]] %>%
  aggregate(fact=15, fun="mean") # aggregate to reduce computation time

cellSize(rs_wc_bioclim, unit='km')

rs_wc_bioclim <- rs_wc_bioclim %>%
  project('+proj=eqearth') %>% # project
  resample(y=CWM[[1]], method='average') # resample

CWM_bioclim <- c(CWM, rs_wc_bioclim) # group rasters

df_CWM_bioclim <- CWM_bioclim %>% as.data.frame(xy=TRUE) %>% # convert into dataframe
  na.omit() %>% subset(CWM_height > 0) # remove NAs and zeros

# par(mfrow=c(1,1), mar=c(5,5,5,5))
# plot(CWM_height ~ bio_1, data=df_CWM_bioclim); abline(lm(CWM_height ~ bio_1, data=df_CWM_bioclim), col='blue', lwd=3)

test_CWM_bioclim <- matrix(ncol=6, nrow=3*6) %>% as.data.frame()
colnames(test_CWM_bioclim) <- c('CWM', 'BIOCLIM', 'F_', 'DF', 'R', 'P_value')
test_CWM_bioclim$CWM <- c(rep('CWM_height',6), rep('CWM_bladearea',6), rep('CWM_nutletvolume',6))
test_CWM_bioclim$BIOCLIM <- rep(c('bio_1','bio_6','bio_7','bio_12','bio_13','bio_15'), 3)

for (i in 1:nrow(test_CWM_bioclim)) {
  
  modtest1 <- modified.ttest(df_CWM_bioclim[,test_CWM_bioclim$CWM[i]], df_CWM_bioclim[,test_CWM_bioclim$BIOCLIM[i]],
                             df_CWM_bioclim[,c('x','y')], nclass=13)
  
  test_CWM_bioclim$F_[i] <- modtest1$Fstat
  test_CWM_bioclim$DF[i] <- modtest1$dof
  test_CWM_bioclim$R[i] <- modtest1$corr
  test_CWM_bioclim$P_value[i] <- modtest1$p.value
  
  print(paste('--- ', round(i/nrow(test_CWM_bioclim),2)*100,'% ---', sep=''))
  
}


