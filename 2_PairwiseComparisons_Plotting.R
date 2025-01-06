### Mouse Humboldt: 
rm(list=ls())

library(humboldt)
library(tcltk)
library(tidyverse)
library(terra)
library(raster)
###MOUSE NICHE COMPARISONS: 

## read in dist data:
Dat<- read.csv("mouse_EnvData_allVars.csv")
Dist<- Dat %>% dplyr::select(GroupCat, Long, Lat)

##read in rasters: (just use the small one to make sure that all of the variables are the same across models)
Rasts<- rast("rast_subset_small.tif")
Env_Rast<- stack(Rasts)

## get these env data into humboldt format:
##convert one of the rasters to a point dataframe to sample.  Use any raster input.
env.points_Full<-rasterToPoints(Env_Rast[[1]], fun=NULL, spatial=TRUE)




##Extract values to points from rasters
## climate: 
EnvPoints_humb<-data.frame(raster::extract(Env_Rast, env.points_Full))
EnvFull_humb<-cbind(env.points_Full@coords,EnvPoints_humb)
names(EnvFull_humb)[5]<-"SNUM"


anyNA(EnvFull_humb)

Env_humb_NONA<- na.omit(EnvFull_humb)


##### Get a df of pairwise species comparisons:
names(Dist)<-c("taxon_name", "lon", 'lat')
SpNames<- Dist %>% distinct(taxon_name)
SpecComb<-t(combn(SpNames$taxon_name,2))
SpecComb<- as.data.frame(SpecComb)
names(SpecComb)[1]<-"s1"
names(SpecComb)[2]<-"s2"

## separate your thing in to a list of occurrence records based on species pairs:
PairsList<-list()
d<-1
for(i in 1:nrow(SpecComb)){
  print(i)
  s1<- SpecComb$s1[i]
  s2<-SpecComb$s2[i]
  temp<- Dist %>% filter(taxon_name== s1 | taxon_name==s2)
  PairsList[[d]]<-temp
  d<- d +1
}

### now write a loop to do the humboldt (and indivdually save each comparison)

## FULL EXTENT: 
stats_Full<- data.frame(D=NA, I=NA, p.D=NA, p.I=NA, num=1:10)

for(i in 1:length(PairsList)){
  TEST<-PairsList[[i]]
  names<- unique(TEST$taxon_name)
  T_1<- TEST %>% filter(taxon_name==names[1])
  T_2<- TEST %>% filter(taxon_name==names[2])
  
  full<-humboldt.doitall(inname=paste(names[1],"_",names[2],"Full", sep=""), Env_humb_NONA, 
                         env2=Env_humb_NONA, sp1=T_1, sp2=T_2, rarefy.dist=5, rarefy.units="km", 
                         env.reso=0.04166667, reduce.env=0, reductype="PCA", non.analogous.environments="YES", 
                         correct.env=T, env.trim=T,  env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, 
                         pcx=1, pcy=2, col.env=e.var, e.var=c(3:9), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="YES",
                         thresh.espace.z=0.0001, p.overlap=F, p.boxplot=F, p.scatter=F, run.silent=T, ncores=10)
  
  stats_Full$D[i]<- paste(full$eqiv$obs$D)
  stats_Full$I[i]<-paste(full$eqiv$obs$I)
  stats_Full$p.D[i]<- paste(full$eqiv$p.D)
  stats_Full$p.I[i]<- paste(full$eqiv$p.I)
  stats_Full$s1[i]<- paste(names[1])
  stats_Full$s2[i]<- paste(names[2])
  
  print(i)
}


### shared e space only. but like. they're fully overlapping so. 
stats_Shared<- data.frame(D=NA, I=NA, p.D=NA, p.I=NA, num=1:10)
for(i in 4:length(PairsList)){
  TEST<-PairsList[[i]]
  names<- unique(TEST$taxon_name)
  T_1<- TEST %>% filter(taxon_name==names[1])
  T_2<- TEST %>% filter(taxon_name==names[2])
  
  shared<-humboldt.doitall(inname=paste(names[1],"_",names[2],"Shared", sep=""), Env_humb_NONA, 
                           env2=Env_humb_NONA, sp1=T_1, sp2=T_2, rarefy.dist=5, rarefy.units="km", 
                           env.reso=0.04166667, reduce.env=1, reductype="PCA", non.analogous.environments="NO", 
                           correct.env=T, env.trim=T,  env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, 
                           pcx=1, pcy=2, col.env=e.var, e.var=c(3:9), R=100, kern.smooth=1, e.reps=100, b.reps=100,
                           thresh.espace.z=0.0001, p.overlap=T, p.boxplot=F, p.scatter=F, run.silent=T, ncores=10)
  
  stats_Shared$D[i]<- paste(shared$eqiv$obs$D)
  stats_Shared$I[i]<-paste(shared$eqiv$obs$I)
  stats_Shared$p.D[i]<- paste(shared$eqiv$p.D)
  stats_Shared$p.I[i]<- paste(shared$eqiv$p.I)
  stats_Shared$s1[i]<- paste(names[1])
  stats_Shared$s2[i]<- paste(names[2])
  
  print(i)
}



### save ur stats: 
write.csv(stats_Full, "NicheComp_FullRange.csv", row.names=F)
write.csv(stats_Shared, "NicheComp_SharedRange.csv", row.names=F)


### make a little heatmap:  

##### FULL ####### (NICHE OVERLAP TEST)
## get mirrored values: 
full_2<-stats_Full
names(full_2)[6:7]<-c("s2", "s1")
tax1<- unique(full_2$s1)

## bind this: 
stats_Full2<- bind_rows(stats_Full, full_2)
stats_Full2$pair<- paste(stats_Full2$s1, stats_Full2$s2, sep="_")

## make a matrix and get just the lower diagonal: 
mat_full1<- dcast(stats_Full2, s1~s2, value.var="D")
mat_full1[upper.tri(mat_full1, diag=F)] <- NA 
rownames(mat_full1)<- mat_full1$s1
mat_full1<- mat_full1 %>% dplyr::select(-s1)

### change up some names and shit:
stats_Full3<- melt(as.matrix(mat_full1))
names(stats_Full3)<-c("s1", "s2", "D_1")
stats_Full3$pair<- paste(stats_Full3$s1, stats_Full3$s2, sep="_")

## combine this with your full one: 
PVals<- stats_Full3 %>% filter(is.na(D_1))
PVals2<- stats_Full2 %>% filter(pair %in% PVals$pair) %>% dplyr::select(pair, p.D)
PVals3<- full_join(PVals, PVals2)
PVals3$D_1<- paste(PVals3$p.D)
PVals3<- PVals3 %>% dplyr::select(s1, s2, pair, D_1)

Lower_full<- stats_Full3 %>% filter(!(is.na(D_1)))

stats_Full4<- rbind(Lower_full, PVals3)

## now add back just the D score: 
D_full<- stats_Full3 %>% dplyr::select(pair, D_1)
names(D_full)[2]<-"D"

plot_full<- full_join(D_full, stats_Full4)
plot_full[plot_full=="NA"]<-NA

plot_full$D_1<- round(as.numeric(plot_full$D_1), digits=3)
plot_full$D_1[is.na(plot_full$D_1)]<-""

###### 
ggplot(plot_full, aes(x = s1, y=s2, fill=as.numeric(D))) + 
  geom_tile()+ 
  coord_fixed()+
  geom_text(aes(label = format(D_1, nsmall=2)))+
  scale_fill_gradient(low = "#fcfdbf", high= "#b73779", na.value= "white")+ 
  xlab("Taxon 1")+
  ylab("Taxon 2")+
  theme_light()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())+
  ggtitle("Niche Overlap Test")



###### SHARED: ############# (NICHE DIVERGENCE TEST)

## get mirrored values: 
Shared_2<-stats_Shared
names(Shared_2)[6:7]<-c("s2", "s1")
tax1<- unique(Shared_2$s1)

## bind this: 
stats_Shared2<- bind_rows(stats_Shared, Shared_2)
stats_Shared2$pair<- paste(stats_Shared2$s1, stats_Shared2$s2, sep="_")

## make a matrix and get just the lower diagonal: 
mat_Shared1<- dcast(stats_Shared2, s1~s2, value.var="D")
mat_Shared1[upper.tri(mat_Shared1, diag=F)] <- NA 
rownames(mat_Shared1)<- mat_Shared1$s1
mat_Shared1<- mat_Shared1 %>% dplyr::select(-s1)

### change up some names and stuff:
stats_Shared3<- melt(as.matrix(mat_Shared1))
names(stats_Shared3)<-c("s1", "s2", "D_1")
stats_Shared3$pair<- paste(stats_Shared3$s1, stats_Shared3$s2, sep="_")

## combine this with your Shared one: 
PVals<- stats_Shared3 %>% filter(is.na(D_1))
PVals2<- stats_Shared2 %>% filter(pair %in% PVals$pair) %>% dplyr::select(pair, p.D)
PVals3<- full_join(PVals, PVals2)
PVals3$D_1<- paste(PVals3$p.D)
PVals3<- PVals3 %>% dplyr::select(s1, s2, pair, D_1)

Lower_Shared<- stats_Shared3 %>% filter(!(is.na(D_1)))

stats_Shared4<- rbind(Lower_Shared, PVals3)

## now add back just the D score: 
D_Shared<- stats_Shared3 %>% dplyr::select(pair, D_1)
names(D_Shared)[2]<-"D"

plot_Shared<- full_join(D_Shared, stats_Shared4)
plot_Shared[plot_Shared=="NA"]<-NA

plot_Shared$D_1<- round(as.numeric(plot_Shared$D_1), digits=3)
plot_Shared$D_1[is.na(plot_Shared$D_1)]<-""
###### 
ggplot(plot_Shared, aes(x = s1, y=s2, fill=as.numeric(D))) + 
  geom_tile()+ 
  coord_fixed()+
  geom_text(aes(label = format(D_1, nsmall=2)))+ 
  scale_fill_gradient(low = "#fcfdbf", high= "#b73779", na.value= "white")+ 
  xlab("Taxon 1")+
  ylab("Taxon 2")+
  theme_light()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())+
  ggtitle("Niche Divergence Test")




### Plot nelsoni and lineatus in g-space: 
Rast<- Env_Rast$alt

### make a nice plot of your raster: 
Rast_df<- as.data.frame(Rast, xy=TRUE)
names(Rast_df)[3]<- "elevation"
Rast_plot<- ggplot(data = Rast_df) +
  geom_raster(aes(x = x, y = y, fill = elevation)) +
  scale_fill_viridis_c(name= "elevation", option="cividis") +
  theme_void()
Rast_plot

### get your df all to one place: 
Lin_nels<- rbind(T_1, T_2)

Rast_plot +
  geom_point(aes(x=as.numeric(lon), y=as.numeric(lat), color=taxon_name), data=Lin_nels)