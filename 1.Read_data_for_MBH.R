###   ###   ###   Read data for MBH analysis    ###   ###   ###


# clear environment ----
rm(list = ls())


# libraries ----
library( rgdal)
library( sp)
library( raster)


# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
o.dir <- paste(w.dir, "outputs", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')



# Read polygons of NPZs and adjacent areas ----
in.hr <- readOGR(paste(s.dir, "SPZ_HillRiver.shp", sep = '/'))
out.hr <- readOGR(paste(s.dir, "South_HillRiver.shp", sep = '/'))
both.hr <- readOGR(paste(s.dir, "In-Out-HillRiver.shp", sep = '/'))


# check crs --
proj4string(in.hr)
proj4string(out.hr)
proj4string(both.hr)


# Prepare list of polygons ----
zones <- list()
zones$in.hr <- in.hr
zones$out.hr <- out.hr
zones$both.hr <- both.hr


#intial look to see area
plot( zones$both.hr, border='black')
plot( zones$out.hr, add=TRUE, col='orange')
plot( zones$in.hr, add=TRUE, col='green')


# Save zones rds ----

saveRDS(zones, file= paste(d.dir, "Zones_Jurien_HillRiver.RDS", sep='/'))


## Read raster data ----
ders <- stack(paste(r.dir, "HillRiver_ders.tif", sep ='/'))
plot(ders)
# set names of derivatives --
names(ders) <- c("depth", "slope", "tpi", "aspect")

b <- ders$depth
s <- ders$slope


# Save rasters rds ----
hr_rasters <- list()
hr_rasters$depth <- b
hr_rasters$slope <- s

saveRDS(hr_rasters, file= paste(d.dir, "HillRiver_rasters.RDS", sep='/'))



## Converting polygons to a common raster ----

###       ###       ### this takes a while for fine res data  ###      ###       ###
in.hr_raster <- rasterize(x=zones$in.hr, y=b, field=zones$in.hr@data[,1], bkg.value=NA, fun="first")
plot(in.hr_raster)
out.hr_raster <- rasterize(x=zones$out.hr, y=b, field=zones$out.hr@data[,1], bkg.value=NA, fun="first")
plot(out.hr_raster)
both.hr_raster <- rasterize(x=zones$both.hr, y=b, field=zones$both.hr@data[,1], bkg.value=NA, fun="first")
plot(both.hr_raster)



#convert and combine --
tmp1 <- as.data.frame( in.hr_raster, xy=TRUE)
tmp2 <- as.data.frame( out.hr_raster, xy=TRUE)
tmp3 <- as.data.frame( both.hr_raster, xy=TRUE)
tmp4 <- as.data.frame( b, xy=TRUE)
tmp5 <- as.data.frame( s, xy=TRUE)


# Join data for Hill River and adjacent analysis --

hrDat <- cbind( tmp1, tmp2[,3])
#hrDat <- cbind( hrDat, tmp3[,3])
hrDat <- cbind( hrDat, tmp4[,3])
hrDat <- cbind( hrDat, tmp5[,3])
head(hrDat)



# Set column names --
df.names <- c("Eastern", "Northing", "in.hr", "out.hr", "depth", "slope")

names(hrDat) <- df.names
names(hrDat) <- df.names


# Save raster dfs rds ----
saveRDS(hrDat, file= paste(d.dir, "hrDat.RDS", sep='/'))







