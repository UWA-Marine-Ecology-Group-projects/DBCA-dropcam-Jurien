###   ###   ###    Prepare spatial data for analysis    ###   ###   ###

# Libraries ----
library(raster)
library(rgdal)
library(sp)
library(ggplot2)
library(pals)
library(RColorBrewer)

# clear environment ----
rm(list = ls())

# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


## Load bathy and park data ----

# load West coast bathy --
b <- raster(paste(r.dir, "CE2016_Mean_Lidar.tif", sep='/'))
plot(b)

crs.utm <- proj4string(b)

# load part of Abrolhos Marine Park --
hr <- readOGR(paste(s.dir, "SPZ_HillRiver.shp", sep='/'))

crs.ll <- proj4string(hr) # +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 

shr <- readOGR(paste(s.dir, "South_HillRiver.shp", sep='/'))

# Join the polygons
hr.both <- raster::union(hr, shr)

writeOGR(hr.both, s.dir, "In-Out-HillRiver", driver = "ESRI Shapefile")

# Reproject lidar ----
b2 <- projectRaster(b, crs=crs.ll)
plot(b2)
plot(hr, add = T)
plot(shr, add = T)

# e <- drawExtent()
e <- extent(114.9628, 115.0831, -30.51398, -30.39573) #  (xmin, xmax, ymin, ymax)
b.cropped <- crop(b2, e)
plot(b.cropped)
plot(hr, add=T); plot(shr, add=T)

# Remove areas above water ----
b.cropped[b.cropped > 0] <- NA
plot(b.cropped)
plot(hr.both, add=T)

# Save new bathy ----
writeRaster(b.cropped, paste(r.dir, "HillRiver_Lidar.tif", sep='/'), overwrite=T)

## Remove regions shallower than 2.5 m  for BOSS ----
b3 <- b.cropped 
b3[b3 > -2] <- NA
plot(b3)
plot(hr.both, add=T)

# Save new bathy ----
writeRaster(b3, paste(r.dir, "HillRiver_Lidar_deep.tif", sep='/'), overwrite=T)


# Separate polys of NPZ ----
# swabrnpz06 <- ab[ab$POLYGONID=="swabrnpz06",]
# #plot(swabrnpz06)
# writeOGR(swabrnpz06, s.dir, "Ab_NPZ06", driver = "ESRI Shapefile")
# 
# swabrnpz09 <- ab[ab$POLYGONID=="swabrnpz09",]
# #plot(swabrnpz09, add=T)
# writeOGR(swabrnpz09, s.dir, "Ab_NPZ09", driver = "ESRI Shapefile")
# 
# npzs <- raster::union(swabrnpz06, swabrnpz09)
# plot(npzs)



# Save --
#writeOGR(npzs, s.dir, "Ab_NPZs", driver = "ESRI Shapefile")


# Bathy  ----
plot(b3)
plot(hr.both, add=T)

hr.slope <- terrain(b3, "slope")
plot(hr.slope, main = "Slope")
plot(hr.both, add=T)

hr.tpi <- terrain(b3, "tpi")
plot(hr.tpi, main = "tpi")
plot(hr.both, add=T)

hr.rou <- terrain(b3, "roughness")
plot(hr.rou, main = "roughness")
plot(hr.both, add=T)


# Slope
scuts <- c(0,0.1,0.2,0.3,0.4,0.5)
pal <- colorRampPalette(c("light blue", "dark blue",  'red'))

raster::plot(hr.slope,
             breaks = scuts,
             col = pal(5),
             main = "slope")
plot(hr.both, add=T)

# TPI
scuts <- c(-1,0,0.5,1.5)
pal <- colorRampPalette(c("light blue", "dark blue",  'red'))

raster::plot(hr.tpi,
             breaks = scuts,
             col = pal(3),
             main = "Depth")
plot(hr.both, add=T)




# Plotting with ggplot ----

# r <- b.09 #raster object
# #preparing raster object to plot with geom_tile in ggplot2
# r_points = rasterToPoints(r)
# r_df = data.frame(r_points)
# head(r_df) #breaks will be set to column "layer"
# r_df$cuts=cut(r_df$WA_500m_Bathy, breaks=c(-40, -60, -80, -100, -120, -140, -160)) #set breaks
# 
# ggplot(data=r_df) + 
#   geom_tile(aes(x=x,y=y,fill=cuts)) + 
#   scale_fill_brewer("depth",type = "seq", palette = "Blues") +
#   coord_equal() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank()) +
#   xlab("Longitude") + ylab("Latitude")


# Load adjacent polys ----
# p06 <- readOGR(paste(s.dir, "Area_next_npz6.shp", sep='/'))
# p09 <- readOGR(paste(s.dir, "Area_next_npz9.shp", sep='/'))


# Join npz and adjacent area polys ----

# # NPZ  --
# plot(b3)
# plot(npz06, add= T)
# plot(p06, add=T)
# 
# 
# npz06 <- spTransform(npz06, proj4string(b))
# zone6 <- union(npz06, p06)
# plot(zone6)
# 
# # crop bathy for zone 6 --
# b.zone6 <- crop(b.cropped, zone6)
# plot(b.zone6)
# 
# # NPZ 9 --
# plot(b.cropped)
# plot(npz09, add= T)
# plot(p09, add=T)
# 
# 
# npz09 <- spTransform(npz09, proj4string(b))
# zone9 <- union(npz09, p09)
# plot(zone9)
# 
# # crop bathy for zone 9 --
# b.zone9 <- crop(b.cropped, zone9)
# plot(b.zone9)



### Calculate slope, aspect, tpi ----

# NPZ 06
slope <- terrain(b3, 'slope')
plot(slope)

aspect <- terrain(b3, 'aspect', unit = 'degrees')
plot(aspect, col=rainbow(100))

tpi <- terrain(b3, 'TPI')
plot(tpi)

depth <- b3

ders <- stack(depth, slope, tpi, aspect)
names(ders) <- c("depth", "slope", "tpi", "aspect")
plot(ders)

# Save derivatives ----

writeRaster(ders, paste(r.dir, "HillRiver_ders.tif", sep ='/'), overwrite = TRUE)



