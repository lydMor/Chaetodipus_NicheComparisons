### reading in dists and gettting important env variables: 
MouseDist_EVERYTHING <- read_excel("MouseDist_EVERYTHING.xlsx")

### make a column for specific group: 
for (i in 1:nrow(MouseDist_EVERYTHING)){
  query<- MouseDist_EVERYTHING$ssp.[i]
  other<-MouseDist_EVERYTHING$species[i]
  if(!(is.na(query))){
    MouseDist_EVERYTHING$GroupCat[i]<-paste(query)
  }
  else{
    MouseDist_EVERYTHING$GroupCat[i]<-paste(other)
  }
}

### now just get your dist data and the group cat: 
Dist<- MouseDist_EVERYTHING %>% dplyr::select(Long, Lat, GroupCat, species, TEST_LOCALE)

##### NOW GET ENVIRONEMENTAL VARIABLES: 
library(raster)
library(terra)


#### 
######
### BIOCLIM
######
###

#### we want higher res bioclim: 
BIO<-getData('worldclim', var='bio', res=2.5)
alt<- getData("worldclim", var="alt", res=2.5)
BIO_all<- stack(BIO, alt)
BIO2<- rast(BIO)

# save this: 
writeRaster(BIO2, "bioclim_2.4.tif", overwrite=T)

Env<- extract(BIO2, Dist[,1:2])
Env<- cbind(Env, Dist)


####
######
### SOIL 
######
###
setwd("~/Documents/Spiranthes?/NICHE/LANDRAST")

Soil<- readOGR("DSMW.shp")

#### get a spdf: 
wgs1984 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

points<- SpatialPoints(coords=Dist[,1:2], proj4string = wgs1984)
pointsSPDF <- SpatialPointsDataFrame(coords=points, data=Dist)
crs(Soil)<-crs(pointsSPDF)

## extract info: 
Senv<- over(pointsSPDF, Soil)

## add what you need
Senv$GroupCat<- pointsSPDF$GroupCat
Senv$Long<- pointsSPDF$Long
Senv$Lat<- pointsSPDF$Lat

## get just the variables you want: 

Senv<- Senv %>% dplyr::select(SNUM, DOMSOI)


## combine these: 

Env2<- bind_cols(Senv, Env)

### save this: 
write.csv(Env2, "mouse_EnvData_allVars.csv", row.names=F)



#### 
######
###  GETTING IMPORTANT VARIABLES FOR EACH SPECIES: 
#######
###

## we know we're going to keep SNUM/DOMSOI for all of them. we just need to decide which bioclim: 
## we'll just use one set of env variables for all models: 
dat<- Env2 %>% dplyr::select(-DOMSOI, -species, -GroupCat, -Long, -Lat, -TEST_LOCALE)

env_cor <- cor(dat, use = "pairwise.complete.obs")

threshold = 0.7  ## the maximum permissible correlation

dists <- as.dist(1 - abs(env_cor))

clust <- hclust(dists, method = "single")
groups <- cutree(clust, h = 1 - threshold)

## Visualize groups:
plot(clust, hang = -1)
rect.hclust(clust, h = 1 - threshold)

## plotting a PCA (to see which variables capture the most variation)
Env_nona<- na.omit(dat)
comb_pr <- prcomp(Env_nona, cor = T)
summary(comb_pr)
?princomp
comb_pr$loadings
biplot(comb_pr)
biplot(comb_pr, c(3,4))
sumP<-summary(comb_pr)
sumP$importance

scores<-comb_pr$rotation
## important variables: 4*,12*,17,8,9,


### now we need to get everything into one raster stack (that means messing with the shapefile): 
# turn our shapefile into a raster: 

## this takes awhile so we'll crop our bioclim to the extent of our system first:
xrange <- c(min(Dist$Long) - 1, max(Dist$Long) + 1)
yrange <- c(min(Dist$Lat) - 1, max(Dist$Lat) + 1)
crop_bio <- crop(BIO_all, c(xrange, yrange))

## create a raster and then resample to get the same extent and resolution as bioclim: 
ext <- floor(extent(crop_bio))
rr <- raster(ext, res=res(crop_bio))
rr <- rasterize(Soil, rr, field=Soil$SNUM)
extent(crop_bio)<-extent(rr)
soi <- resample(rr, crop_bio, method = "bilinear")
Rasts<-stack(crop_bio, soi)


writeRaster(Rasts, "EnvVarsAll_2.5.tif")
Rasts<-rast(Rasts)


### 1: Bigger option: (more variables): # bio 1, 2, 3, 4, 12, 9, 8, 14, 15, SNUM 

Rasts1<- subset(Rasts, subset=c('bio1', 'bio2', 'bio3', 'bio4', 'bio12', 'bio9', 
                                 'bio8', 'bio14', 'bio15', 'alt', 'layer'))

writeRaster(Rasts1, "Rast_subset_big.tif")

##### 2: Smaller option: bio 2, 4, SNUM, alt, 12, 9, 17

Rasts2<- subset(Rasts, subset=c('bio2', 'bio4', 'layer', 'alt', 'bio12', 'bio9', 'bio17'))

writeRaster(Rasts2, "Rast_subset_small.tif")
      ## SAVE YOUR WORK: 
save.image("~/Documents/PocketMouse/ENM_upto151.RData")

### NOW- get rid of all the shit you don't need. 
rm(list=ls())
## load in the stuff you do need: 
Rasts<- rast("EnvVarsAll_2.5.tif")
Rasts1<- rast("Rast_subset_big.tif")
Rasts2<- rast("Rast_subset_small.tif")
Dat<- read.csv("mouse_EnvData_allVars.csv")

### you can make a loop, or you can individually do each species. There's only like 4 so. IDK. 

### do it manually: (species groups)

collis<- Dat %>% filter(species=="collis")

durangae<- Dat %>% filter(species=="durangae")

lineatus<- Dat %>% filter(species=="lineatus")

nelsoni<- Dat %>% filter(species=="nelsoni")

#subspeciesGroups: 
cc<- Dat %>% filter(GroupCat=="collis")
cm<- Dat %>% filter(GroupCat=="mapimienss")

## save this: 
save.image("Species_DistDataAndEnvVars.RData")


######                   ######
########                  ####
##        COLLIS FULL     #######
#######                  #######

### env variables version one: 
Env2<- terra::extract(Rasts1, collis[,23:24])
remove_na <- which(is.na(Env2), arr.ind=TRUE)[,1]
Env2_nona <- Env2[-remove_na,]
#is empty
comb_pr <- princomp(Env2, cor = T)
source("envSample.R")
points_R2 <- envSample(points[-remove_na, 1:2], filters =
                         list(comb_pr$scores[,1], comb_pr$scores[, 2], comb_pr$scores[, 3]),
                       res = list(0.25, 0.25, 0.25), do.plot = TRUE)
dim(points_R2)
## 610 points



#############################


