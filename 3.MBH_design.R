###   ###   ###   MBH design clustered BRUVs    ###   ###   ###


# libraries ----
#install.packages("MBHdesign")
library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)
library( raster)
library( rgdal)
library( sp)
library( rgeos)


# clear environment ----
rm(list = ls())

# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
#p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')
o.dir <- paste(w.dir, "outputs", sep='/')

####    SET INFO   ####

study <- "HillRiver"

total.no.deployments <- "50_deployments"

design.version <- "v1"

design.info <- paste(study, total.no.deployments, design.version, sep='.')
small.info <- paste(study, design.version, sep='.')


# Read in the inclusion probs ----
slope.cuts <- raster(paste(r.dir, paste("slope_cuts", small.info, "tif", sep ='.'), sep='/'))
plot(slope.cuts)
inclProbs <- raster(paste(r.dir, paste("inclProbs", design.info, "tif", sep='.'), sep='/'))
plot(inclProbs)
inclProbs <- setValues( inclProbs, values( inclProbs) / sum( values( inclProbs), na.rm=TRUE))
plot(inclProbs)

# check sun of incl probs --
cellStats(inclProbs, 'sum')

rootInclProbs <- inclProbs

rootInclProbs <- setValues( rootInclProbs, sqrt( values( rootInclProbs)))

cellStats(rootInclProbs, 'sum')
plot(rootInclProbs)


# Read data ----
zones <- readRDS(paste(d.dir, "Zones_Jurien_HillRiver.RDS", sep='/')) # this one in different folder
#Deans <- readRDS( "DeansPoints_forinNOutMP-d3.RDS")
rast <- readRDS(paste(d.dir, "HillRiver_rasters.RDS", sep='/'))
#if( class( BRUVS) != "SpatialPointsDataFrame")
#Deans <- SpatialPointsDataFrame( coords=Deans[,c("Longitude","Latitude")], data=Deans, proj4string = CRS( proj4string( zones[[1]])))
#proj4string(Deans) <- proj4string(swrast$bathy)
straw.nums <- readRDS(paste(d.dir, "StrawmanNumbers_hr.50_deployments.v1.RDS", sep ='/'))
names(straw.nums) <- c("in.hr" , "out.hr")
names(straw.nums)
straw.nums


############################
####  Spatial sample of new sites ----
####  from altered incl. probs.
############################

### Here use quasiSamp to get random points ####
## these points will be the center of buffer for transects ###

####  Set the seed for reproducability
#set.seed( 777)
#### HAVE NOT BEEN ABLE TO MAKE THIS FUNCTION WORK ----
newSites <- list(in.hr = NULL, out.hr = NULL)

for( zz in c("in.hr", "out.hr")){
  print( zz)
  #the number of samples to take (specified minus the legacy number)
  #numby <- floor( (straw.nums[zz])/4)  # for clustered cluster - without legacy sites
  numby <- floor( (straw.nums[zz])) # for not clustered sites
  #numby <- floor( (straw.nums[zz] - numRef[zz])/2)
  #numby <- floor( (straw.nums[zz] - numRef[zz])) # with legacy sites 
  #set up spatial domain
  myZone <- zones[[zz]]
  #if( zz == "AMP"){
  # myZone = zones$AMP - zones$IUCN2
  #set.seed( 747)
  #}
  #tmpIP <- mask( rootInclProbs, myZone)
  tmpIP <- mask( inclProbs, myZone)
  tmpIP <- crop( tmpIP, myZone)
  #take the sample of clusters based on root incl probs
  newSites[[zz]] <- quasiSamp( n=numby, potential.sites=coordinates( tmpIP), inclusion.probs=values(tmpIP), nSampsToConsider=5000)
  
  #plotting (maybe remove at a later date?)
  tmpIPFull <- mask( inclProbs, myZone)
  tmpIPFull <- crop( tmpIPFull, myZone)
  plot( tmpIPFull)
  #plot( legacySites, add=TRUE, pch=1, col='red')
  points( newSites[[zz]][,c("x","y")], pch=20, col='black')
}
newSites <- do.call( "rbind", newSites)
head(newSites)

# Give id to sites and zones --
site.names <-row.names(newSites)
newSites$site <- as.factor(site.names)
zone.names <- gsub('.{5}$', '', site.names) # remove last 3 characters
zone.names <- gsub('in.', 'in', zone.names)
zone.names <- gsub('out.', 'out', zone.names)
zone.names
#zone.names <- substr(site.names, 1, 3) # extract first three characters
newSites$zone <- as.factor(zone.names)
newSites$zone
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#some of the spatial balance is not great...  Presumably because the balance of the reference sites is also not great...

# Plot --
plot(inclProbs)
plot(slope.cuts, main = "Slope strata")
plot(rast$slope, main = "Slope")
plot(zones$both, add=T)
plot(newSites, col=newSites$zone, pch = 20, add=T) # 41
newSites$zone



### Make sure the clusters centres are ~ 1 km apart ----


## Get CRS in utm ----
crs1 <- CRS("+init=epsg:32750") # WGS 84 / UTM zone 50S


## transform the points into UTM --
p1u <- spTransform(newSites, crs1)

## calculate if 2 points fall within 1500 m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

dist1 <- gDistance(p1u, byid =T)
dist1
max(dist1)
min(dist1[dist1 > 0]) # minimum distance other than 0

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 50, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
p1u[v1, ] # 98 features left

remaining.sites <- p1u[v1, ]
remaining.sites <- spTransform(remaining.sites, proj4string(inclProbs))

# plot --
plot(inclProbs)
plot(rast$slope, main = "Slope")
plot(zones$both, add=T)
plot(remaining.sites, col=remaining.sites$zone, pch = 20, add=T) # 41
remaining.sites$zone

## Save --
site <- "Jurien"
design <- "50Boss"
version <- "v1"

writeOGR(remaining.sites, o.dir, paste(site, design, version, sep='-'), driver = "ESRI Shapefile", overwrite = T)
