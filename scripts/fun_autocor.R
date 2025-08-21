

# Calculate Moran.I


library(tidyverse)
library(terra)
library(ape) # moran.i
library(geosphere) # distm


# autocor takes a spatraster and dataframe of occurrences and calulates Moran's I
# https://cran.r-project.org/web/packages/ape/vignettes/MoranI.pdf


autocor <- function(raster = NULL, xy = NULL){
  
  # extract values from the predictor variable for the points of presence
  xy.variables <- na.omit(
    data.frame(xy,
               terra::extract(x=raster, y=xy, df=TRUE, cells=FALSE, ID=FALSE)
    )
  )
  
  # rename variable to call it later
  colnames(xy.variables)[!(colnames(xy.variables)%in%c('x','y'))] <- 'var'
  
  # matrix of inverse distances between occurrences
  xy.distancias <- 1/ (geosphere::distm(x=xy.variables[, c("x","y")], fun=distGeo) / 1000)
  xy.distancias[!is.finite(xy.distancias)] <- 0 # replace Inf with zero
  diag(xy.distancias) <- 0 # diagonal to 0
  
  # Moran's I 
  moran.i <- ape::Moran.I(xy.variables$var, xy.distancias, na.rm = TRUE)
  # print(moran.i)
  
}


