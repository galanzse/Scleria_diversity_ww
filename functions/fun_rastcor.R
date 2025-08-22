

# rastcor identifies uncorrelated variables in a multi-layer raster following two methods:
# 1/ calculating a given correlation coefficient from c("pearson", "kendall", "spearman") and a threshold
# 2/ using the variance inflation factor


library(tidyverse)
library(terra)
library(corrplot) # corrplot
library(usdm) # vif


temp <- rast('SDM/agg_predictors_MaxEnt.tiff') %>% crop(vect(ne_countries(country='Spain')[1]))
plot(temp) 

raster=temp
method='pearson'
threshold.cor=0.8
threshold.vif=0.8

rastcor <- function(raster=NULL, method=NULL, threshold.cor=NULL, threshold.vif=NULL) {
  
  # FIRST METHOD: CORRELATION
  
  # calculate correlation matrix
  var.cor <- raster %>% as.data.frame() %>% cor(method=method, use='pairwise.complete.obs')
  
  # graphical display of a correlation matrix
  # print(corrplot(var.cor, type="upper", method="number", tl.cex=1, cl.cex=1, cl.ratio=0.1, col=COL2('RdBu', 10)))
  
  # matrix to dataframe
  cor.df <- as.data.frame(var.cor)

  # keep upper tri
  lower <- var.cor
  lower[lower.tri(var.cor, diag=TRUE)] <- ""
  lower.df <- as.data.frame(lower)

  # correlation matrix to distance matrix
  var.dist <- abs(as.dist(cor.df))

  # build dendrogram: distance is inversely proportional to the value of the correlation coefficient (shorter distance = higher correlation)
  var.cluster <- hclust(1-var.dist)

  # plot
  print(plot(var.cluster))
  abline(h=1-threshold.cor, lty=2, lwd=2, col="red")
  
  
  
  # SECOND METHOD: VIF
  
  # usdm::vif(as.data.frame(raster), size=5000)
  print(usdm::vifstep(as.data.frame(raster), th=threshold.vif))

}
 

