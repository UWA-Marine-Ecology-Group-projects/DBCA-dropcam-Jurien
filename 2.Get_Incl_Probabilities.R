###   ###   ###   Get inclusion probabilities   ###   ###   ###

####   NPZ 6   ####


# libraries ----
library( rgdal)
library( sp)
library( raster)

# clear environment ----
rm(list = ls())

# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "outputs", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


# Read in data ----
hrDat <- readRDS(paste(d.dir, "hrDat.RDS", sep='/'))
hr_rasters <- readRDS(paste(d.dir, "HillRiver_rasters.RDS", sep='/'))
zones <- readRDS(paste(d.dir, "Zones_Jurien_HillRiver.RDS", sep='/'))


# Straw man for numbers of samples in each region ----

study <- "HillRiver"

total.no.deployments <- "50_deployments"

design.version <- "v1"

straw.nums <- c(25, 25)  # numbers of drops in and out
straw.props <- straw.nums / sum( straw.nums) # 0.75 0.25
names( straw.nums) <- names( straw.props) <- c("in.hr", "out.hr")
saveRDS( straw.nums, file=paste(d.dir, paste("StrawmanNumbers_hr", total.no.deployments, design.version, "RDS", sep='.'), sep='/'))


# Get slope cut points ----
# and their numbers of drops
#slope.quant <- c(0,0.5,0.9,0.98,1) # version 1
slope.quant <- c(0,0.1,0.5,0.95,1) # version 2
slope.cuts <- quantile(hr_rasters$slope, slope.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
slope.cuts
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( slope.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4) #floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
slope.targetNums <- rep(floor( ( tmp / sum( tmp))[-1] * sum(straw.nums))) #rep( 40, 5)#c( 20,20,30,65,65)
#slope.targetNums <- rep( floor( 80/8), 2) # 8 8
slope.targetProps <-  slope.targetNums / sum( slope.targetNums) # 0.5 0.5


# Proportion of potential sites in each zone ----
hr_small <- hrDat[!is.na( hrDat$slope),]
tmp <- colSums( hr_small[,c("in.hr", "out.hr")], na.rm=TRUE) # number of cells in npz and out
#tmp[2] <- tmp[2] - tmp[1] # so similar amount of sites in SPZ and MUZ
tmp[1] # 7495
tmp[2] # 2676 
props <- tmp / nrow( hr_small) # inside 0.7222222 - outside 0.2777778
props <- props / sum( props) # inside 0.7222222 - outside 0.2777778  



# To get cut points ----
catB <- cut( hr_rasters$slope, breaks=slope.cuts, na.rm=TRUE)
plot(catB)
plot( zones$in.hr, add=T); plot( zones$out.hr, add=TRUE)


#writeRaster(catB, paste(r.dir, paste('slope_cuts', total.no.deployments, design.version, 'tif', sep='.'), sep='/'), overwrite=TRUE)



# Get inclusion probabilities within each zone ----
# Weight according to straw.props 

inclProbs <- catB
for( zz in c( "in.hr", "out.hr")){
  print( zz)
  #if( zz == "MUZ")
  #zoneID <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
  #zoneID <- extract( x=catB, y=zones$MUZ-zones$NPZ, cellnumbers=TRUE)
  #else
  zoneID <- extract( x=catB, y=zones[[zz]], cellnumbers=TRUE)
  propsOfslope <- table( catB@data@values[zoneID[[1]][,"cell"]])
  propsOfslope <- propsOfslope / sum( propsOfslope)
  tmp <- slope.targetProps / propsOfslope #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfslope)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
}

inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  # cheats way to crop
plot( inclProbs)


# Standardising so that the zone totals are correct according to straw.props | straw.nums ----
cells.in.hr <- extract( x=catB, y=zones$in.hr, cellnumbers=TRUE)
cells.out.hr <- extract( x=catB, y=zones$out.hr, cellnumbers=TRUE)


inclProbs@data@values[cells.in.hr[[1]][,'cell']] <- inclProbs@data@values[cells.in.hr[[1]][,'cell']] * straw.props["in.hr"]
inclProbs@data@values[cells.out.hr[[1]][,'cell']] <- inclProbs@data@values[cells.out.hr[[1]][,'cell']] * straw.props["out.hr"]

#inclProbs@data@values[HPZZone[[1]][,'cell']] <- inclProbs@data@values[HPZZone[[1]][,'cell']] * straw.props["HPZ"]
#inclProbs@data@values[NPZZone[[1]][,'cell']] <- inclProbs@data@values[NPZZone[[1]][,'cell']] * straw.props["NPZ"]

plot(inclProbs)


writeRaster( inclProbs, paste(r.dir, paste('inclProbs', study, total.no.deployments, design.version, 'tif', sep ='.'), sep='/'), overwrite=TRUE)
