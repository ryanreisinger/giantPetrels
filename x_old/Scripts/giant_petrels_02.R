#### Giant Petrels

## Ryan Reisinger

library(ggplot2)
library(argosfilter)
library(geosphere)
library(rgdal)
library(raster)
library(marmap)
library(cowplot)
library(orsifronts)
# library(rgeos)
# library(maptools)
library(pals)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

#-----------------------------------
## Stuff for figures
# Widths (in inches)
single.col <- 84*0.0393701
double.col <- 174*0.0393701
double.col.sup <- 150*0.0393701
fig.scale <- 8/11

source("./Scripts/theme_rr.R")

#-----------------------------------
## Get metadata
met <- read.csv("./Data/Metadata_2015-16_02_v3.csv", stringsAsFactors = F)
met <- met[met$Species == "Northern giant petrel" | met$Species == "Southern giant petrel", ]
met <- met[met$File.name != "", ]
met$File.name <- paste0(met$File.name, ".csv")


# Convert lat and lon to decimal degrees
dd <- function(string){
  d <- strsplit(x = string, split = " ")
  d <- unlist(d)
  deg <- as.numeric(d[1])
  min <- as.numeric(d[2])/60
  p <- deg + min
}

met$Lat <- unlist(lapply(FUN = dd, X = met$Lat))
met$Lat <- met$Lat * -1 #assume all positions are South
met$Lon <- unlist(lapply(FUN = dd, X = met$Lon))

# Remove individual which goes nowhere
rm <- "GPS_SGP_151004_INC_KD_0915.csv"
met <- met[met$File.name != rm, ]

# More accurate deployment locations
met[met$File.name == "GPS_SGP_150925_INC_KD_0908.csv", "Lat"] <- -46.96359
met[met$File.name == "GPS_SGP_150925_INC_KD_0908.csv", "Lon"] <- 37.85200

met[met$File.name == "GPS_NGP_151001_INC_KD_1019.csv", "Lat"] <- -46.95483
met[met$File.name == "GPS_NGP_151001_INC_KD_1019.csv", "Lon"] <- 37.86412

met[met$File.name == "GPS_NGP_150916_INC_KD_1003.csv", "Lat"] <- -46.952148
met[met$File.name == "GPS_NGP_150916_INC_KD_1003.csv", "Lon"] <- 37.858187

met[met$File.name == "GPS_SGP_150925_INC_KD_0906.csv", "Lat"] <- -46.963811
met[met$File.name == "GPS_SGP_150925_INC_KD_0906.csv", "Lon"] <- 37.851975

# KML plot of nest locations
mnt <- met
coordinates(mnt) <- c("Lon", "Lat")
crd <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(mnt) <- crd
writeOGR(mnt,
         dsn = paste0("./Output/KML/Nests.kml"),
         layer = "Nests",
         driver = "KML",
         overwrite_layer = TRUE)


# Sex
# Break at ~96 mm
met$Culmen_Length <- as.numeric(met$Culmen_Length)
met$Culmen_Depth <- as.numeric(met$Culmen_Depth)

# Update culmen depth of one individual
# from metadata in Data_original
met[met$ID == "SGP04_KD_SEP_2015", "Culmen_Depth"] <- 37.6

# Culmen Length threshold
v <- 96

for (i in 1:nrow(met)) {
  if (is.na(met$Culmen_Length[i])) {
    met$sex[i] <- "Unknown"
  } else if (met$Culmen_Length[i] < v) {
    met$sex[i] <- "Female"
  } else {
    met$sex[i] <- "Male"
  }
}

# Add results from molecular sexing
met$sex.m <- rep(NA, nrow(met))
met[met$ID == "NGP02_KD_SEP_2015", "sex.m"] <- "Female"
met[met$ID == "NGP14_KD_SEP_2015", "sex.m"] <- "Male"
met[met$ID == "NGP18_KD_SEP_2015", "sex.m"] <- "Female"
met[met$ID == "NGP19_KD_SEP_2015", "sex.m"] <- "Male"

met[met$ID == "SGP01_KD_SEP_2015", "sex.m"] <- "Male"
met[met$ID == "SGP02_KD_SEP_2015", "sex.m"] <- "Male"
met[met$ID == "SGP04_KD_SEP_2015", "sex.m"] <- "Female"
met[met$ID == "SGP18_KD_SEP_2015", "sex.m"] <- "Female"
met[met$ID == "SGP20_KD_SEP_2015", "sex.m"] <- "Female"

# All known sexes match genetics results, so replace "Unknown" sex
met[met$ID == "NGP19_KD_SEP_2015", "sex"] <- "Male"
met[met$ID == "SGP02_KD_SEP_2015", "sex"] <- "Male"


#-----------------------------------
# Read in files
files <- met$File.name

# Check structure
# Some GPS files were still semicolon delimited
t <- matrix(ncol = 2, nrow = length(files))
t <- as.data.frame(t)
names(t) <- c("File.name", "Ncol")

for (i in 1:length(files)) {
  a <- read.csv(paste0("./Data/GPS/", files[i]), stringsAsFactors = F)
  t$File.name[i] = files[i]
  t$Ncol[i] = ncol(a)
}

rm(t)

dt <- matrix(ncol = 11)
dt <- as.data.frame(dt)
names(dt) <- c("Date", "Time", "Latitude", "Longitude", "Altitude", "Speed", "Course", "Type", "Distance", "Essential", "File.name")

for (i in 1:length(files)) {
  a <- read.csv(paste0("./Data/GPS/", files[i]), stringsAsFactors = F)
  a$File.name <- rep(files[i], nrow(a))
  dt <- rbind(dt, a)
}

dt <- dt[-1, ]

#-----------------------------------
## Speed filter
# Adjust speed (m/s) with vmax argument

ids <- unique(dt$File.name) #'ids' should be same as 'files' created above

fn <- matrix(ncol = 13)
fn <- as.data.frame(fn)
names(fn) <- c("Date", "Time", "Latitude", "Longitude", "Altitude", "Speed", "Course", "Type", "Distance", "Essential", "File.name", "date.time", "Speed.flag")
fn$date.time <- as.POSIXct(fn$date.time, format = "%Y/%m/%d %H:%M:%S", tz = "GMT") # Set class to avoid conversion in rbind later

for (i in 1:length(ids)) {
  f <- dt[dt$File.name == ids[i], ]
  f$date.time <- paste(f$Date, f$Time, sep = " ")
  f$date.time <- as.POSIXct(f$date.time, format = "%Y/%m/%d %H:%M:%S", tz = "GMT")
  f$Speed.flag <- vmask(lat = f$Latitude,
                        lon = f$Longitude,
                        dtime = f$date.time,
                        vmax = 30)
  f <- f[f$Speed.flag != "removed", ]
  fn <- rbind.data.frame(fn, f)
}

fn <- fn[-1, ]

#-----------------------------------
## Merge metadata and tracking data

dat <- merge(fn, met, by = "File.name")

#-----------------------------------
## Quick plot as sanity check

# Straight ggplot
world <- borders("world", colour = "gray50", fill = "gray50", xlim = c(min(dat$Longitude), max(dat$Longitude)), ylim = c(min(dat$Latitude), max(dat$Latitude)))
ggplot(data = dat, aes(x = dat$Longitude, y = dat$Latitude, colour = Species)) + geom_point() + coord_quickmap() + world

#-----------------------------------
## Distance from deployment

dat$Distance <- distGeo(p1 = cbind(dat$Lon, dat$Lat), p2 = cbind(dat$Longitude, dat$Latitude))
dat$Distance <- dat$Distance/1000

for (i in 1:length(ids)) {
  d <- dat[dat$File.name == ids[i], ]
  p <- ggplot(data = d, aes(x = date.time, y = Distance, colour = sex)) + geom_point() + labs(title = ids[i]) 
  print(p)
}

#-----------------------------------
## Maximum distance

l <- dat
l <- l[ -(1:nrow(l)), ]
for (i in 1:length(ids)) {
  tm <- dat[dat$File.name == ids[i], ]
  m <- max(tm$Distance)
  tm$Maxdist <- rep(m, nrow(tm))
  l <- rbind.data.frame(l, tm)
}

dat <- l # Replace dat

#-----------------------------------
## KML files for closer inspection

for (i in 1:length(ids)) {
  mnt <- dat[dat$File.name == ids[i], ]
  coordinates(mnt) <- c("Longitude", "Latitude")
  crd <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  proj4string(mnt) <- crd
  writeOGR(mnt, dsn = paste0("./Output/KML/", ids[i], ".kml"), layer = ids[i], driver = "KML")
}

#-----------------------------------
## Trim start and end portions of tracks using locator

# No need to re-run this,
# pre-run output can be read in

# tr <- dat[1, ]
# tr$trim <- rep(NA, 1)
# tr<- tr[-1, ]
# 
# for (i in 1:length(ids)){
#   test <- dat[dat$File.name == ids[i], ]
#   test$d <- as.numeric(test$date.time)
#   plot(test$date.time, test$Distance, pch = 19)
#   t <- locator(2)$x
#   for (j in 1:nrow(test)){
#     if (test$d[j] > t[1] && test$d[j] < t[2]) {
#       test$trim[j] <- "Keep"
#     } else {
#       test$trim[j] <- "Trim"
#     }
#   }
#   test$d <- NULL
#   tr <- rbind.data.frame(tr, test)
# }
# 
# tr <- tr[tr$trim == "Keep", ] #Remove the trimmed sections
# 
# # Save intermediate output
# saveRDS(tr, file = "./Output/tr_01.RDS")

# Read in pre-run
tr <- readRDS("./Output/tr_01.RDS")

# Plot to check
for (i in 1:length(ids)) {
  d <- tr[tr$File.name == ids[i], ]
  p <- ggplot(data = d, aes(x = date.time, y = Distance, colour = sex)) + geom_point() + geom_line() + labs(title = ids[i])
  print(p)
}

# Overwrite KMLs
for (i in 1:length(ids)) {
  mnt <- tr[tr$File.name == ids[i], ]
  coordinates(mnt) <- c("Longitude", "Latitude")
  crd <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  proj4string(mnt) <- crd
  writeOGR(mnt, dsn = paste0("./Output/KML/", ids[i], ".kml"), layer = paste0(ids[i], "_trim"), driver = "KML")
}

## KMLs for males, by species
mnt <- tr[tr$Species == "Northern giant petrel" & tr$sex == "Male", ]
coordinates(mnt) <- c("Longitude", "Latitude")
crd <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(mnt) <- crd
writeOGR(mnt, dsn = paste0("./Output/KML/", "NGP_males", ".kml"),
         layer = "NGP_males", driver = "KML",
         overwrite_layer = TRUE)

mnt <- tr[tr$Species == "Southern giant petrel" & tr$sex == "Male", ]
coordinates(mnt) <- c("Longitude", "Latitude")
crd <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(mnt) <- crd
writeOGR(mnt, dsn = paste0("./Output/KML/", "SGP_males", ".kml"),
         layer = "SGP_males", driver = "KML",
         overwrite_layer = TRUE)

#-----------------------------------
## Define trips

## When at nest
tr$nest <- rep("nest", nrow(tr))
tr[tr$Distance > 0.15, "nest"] <- "trip"

# No need to re-run this,
# pre-run output can be read in

#Keep trips together for the moment, later discard positions near nest
# tr2 <- tr[1, ]
# tr2$trip <- rep(NA, 1)
# tr2<- tr2[-1, ]
# 
# 
# for (j in 1:length(ids)){
#   test <- tr[tr$File.name == ids[j], ]
#   test$d <- as.numeric(test$date.time)
#   plot(test$date.time, test$Distance, type = "l") +
#   points(test$date.time, test$Distance, pch = 19, col = as.factor(test$nest))
#   abline(h=0.1)
#   t <- locator(9)$x
#   # t <- t[!duplicated(t)]
#   for (i in 1:nrow(test)){
#     print(i)
#     if (test$d[i] > t[1] && test$d[i] < t[2]){
#       test$trip[i] <- 1
#     } else if (test$d[i] > t[2] && test$d[i] < t[3]) {
#       test$trip[i] <- 2
#     } else if (test$d[i] > t[3] && test$d[i] < t[4]) {
#       test$trip[i] <- 3
#     } else if (test$d[i] > t[4] && test$d[i] < t[5]) {
#       test$trip[i] <- 4
#     } else if (test$d[i] > t[5] && test$d[i] < t[6]) {
#       test$trip[i] <- 5
#     } else if (test$d[i] > t[6] && test$d[i] < t[7]) {
#       test$trip[i] <- 6
#     } else if (test$d[i] > t[7] && test$d[i] < t[8]) {
#       test$trip[i] <- 7
#     } else if (test$d[i] > t[8] && test$d[i] < t[9]) {
#       test$trip[i] <- 8
#     } else {
#       test$trip[i] <- NA
#     }
#   }
#   tr2 <- rbind.data.frame(tr2, test)
# }
# 
# tr2$trip2 <- paste0(tr2$File.name, "_", tr2$trip)
# 
# saveRDS(tr2, file = "./Output/tr_02.RDS")

# Read in pre-run
tr2 <- readRDS("./Output/tr_02.RDS")

#-----------------------------------
## Trim again
tr3 <-tr2[tr2$Distance > 0.15,] # Remove locations less than 100 m from putative nest
saveRDS(tr3, file = "./Output/tr_03.RDS")

#-----------------------------------
## Displacement plots

tr4 <- tr3
tr4$lag <- rep(NA, nrow(tr4))
tr4 <- tr4[1, ]
tr4 <- tr4[-1, ]

trips <- unique(tr3$trip2)

for (i in 1:length(trips)) {
  trp <- tr3[tr3$trip2 == trips[i], ]
  tms <- diff(trp$date.time)
  units(tms) <- "hours"
  tms <- c(0, tms)
  
  tms <- cumsum(tms)
  trp$lag <- tms
  tr4 <- rbind.data.frame(tr4, trp)
}

tr4$ND <- tr4$Distance^2

tr4sub <- tr4

# Long trips only
tr4long <- tr4[tr4$Maxdist > 100, ]
disp <- ggplot(tr4long, aes(x = lag/24, y = Distance, group = trip2, colour = sex)) +
  geom_line() +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"),
                      name = "Sex") + 
  theme_rr() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(colour = "black"),
                     plot.title = element_text(size = rel(1))) +
  labs(title = "Distant trips",
       x = "Days since deployment", y = "Distance from nest (km)") #+
# scale_x_continuous(limits = c(0, 17))
disp


# Short trips only
tr4short <- tr4[tr4$Maxdist < 100, ]

disp2 <- ggplot(tr4short, aes(x = lag/24, y = Distance, group = trip2, colour = sex)) +
  geom_line() +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"),
                      name = "Sex") + 
  theme_rr() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(colour = "black"),
                     plot.title = element_text(size = rel(1))) +
  labs(title = "Nearby trips",
       x = "Days since deployment", y = "Distance from nest (km)") +
  scale_x_continuous(limits = c(0, 17))
disp2

pdf("./Plots/displacementPlot.pdf",
    width = double.col/fig.scale, height = (double.col.sup/1.6)/fig.scale,
    useDingbats = FALSE)
cowplot::plot_grid(disp, disp2, nrow = 2)
dev.off()

# Combined plot

tr4sub.copy <- tr4sub
tr4sub.copy[tr4sub.copy$Maxdist > 100, "which"] <- "Long trips"
tr4sub.copy[tr4sub.copy$Maxdist < 100, "which"] <- "Short trips"

dispBoth <- ggplot(tr4sub.copy, aes(x = lag/24, y = Distance, group = trip2, colour = sex)) +
  geom_line() +
  facet_wrap(~which, nrow = 2, scales = "free") +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"),
                      name = "Sex") + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(colour = "black")) +
  labs(x = "Days since deployment", y = "Distance from nest (km)") +
  scale_x_continuous(limits = c(0, 17))

pdf("./Plots/displacementPlotFacet.pdf",
    width = double.col/fig.scale, height = (double.col.sup/1.6)/fig.scale,
    useDingbats = FALSE)
dispBoth
dev.off()

# Numbers in each category
length(unique(tr4sub.copy[tr4sub.copy$sex == "Male" & tr4sub.copy$which == "Short trips", "ID"]))
length(unique(tr4sub.copy[tr4sub.copy$sex == "Female" & tr4sub.copy$which == "Short trips", "ID"]))

length(unique(tr4sub.copy[tr4sub.copy$sex == "Male" & tr4sub.copy$which == "Long trips", "ID"]))
length(unique(tr4sub.copy[tr4sub.copy$sex == "Female" & tr4sub.copy$which == "Long trips", "ID"]))

#-----------------------------------
# Differences in trip distance or bill dims?

test <- aggregate(cbind(lag, Distance) ~ trip2, FUN = max, data = tr4sub)
test$lag <- test$lag/24
test$id <- substr(start = 1, stop = 30, x = test$trip2)
hld <- unique(tr4sub[ , c("ID", "trip2", "sex", "Species", "Culmen_Length")])
test <- merge(x = test, y = hld)

mean(test[test$Distance > 100 & test$sex == "Male", "lag"])
mean(test[test$Distance > 100 & test$sex == "Female", "lag"])

t.test(test[test$Distance > 100 & test$sex == "Female", "lag"], test[test$Distance > 100 & test$sex == "Male", "lag"])

# Culmen length?
test2 <- test[ , -which(names(test) == "trip2")]
test2 <- unique(test2)
mean(test2[test2$Distance > 100 & test2$sex == "Male", "Culmen_Length"], na.rm = T)
mean(test2[test2$Distance < 100 & test2$sex == "Male", "Culmen_Length"], na.rm = T)

#-----------------------------------
## Utilization distributions

library(adehabitatHR)

s <- tr
s <- s[s$sex != "Unknown", ]
s <- s[s$Distance > 0.4, ]
s.males <- s[s$sex == "Male", ]
s.females <- s[s$sex == "Female", ]

# Grid for UD
lms <- c(min(s$Longitude) - 5, max(s$Longitude) + 5, min(s$Latitude) - 5, max(s$Latitude) + 5)
rt <- raster(ext = extent(lms), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), res = 0.05)
rt.sp <- as(rt, "SpatialPixelsDataFrame")

# Group sex x species
s$group <- paste0(s$Species, " ", s$sex)
s <- s[ , c("Longitude", "Latitude", "group")]
coordinates(s) <- c("Longitude", "Latitude")
proj4string(s) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud <- kernelUD(s, h = 1, grid = rt.sp)
image(kud)

# Create rasters of each
vud <- getvolumeUD(kud)
vud.males.raster.ngp <- raster(as(vud$'Northern giant petrel Male',"SpatialPixelsDataFrame"))
vud.males.raster.sgp <- raster(as(vud$'Southern giant petrel Male',"SpatialPixelsDataFrame"))
vud.females.raster.ngp <- raster(as(vud$'Northern giant petrel Female',"SpatialPixelsDataFrame"))
vud.females.raster.sgp <- raster(as(vud$'Southern giant petrel Female',"SpatialPixelsDataFrame"))

# Write the rasters
writeRaster(vud.males.raster.ngp, "./Output/vud.males.raster.ngp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.males.raster.sgp, "./Output/vud.males.raster.sgp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.females.raster.ngp, "./Output/vud.females.raster.ngp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.females.raster.sgp, "./Output/vud.females.raster.sgp.grd", format = "raster",
            overwrite = T)


# By sex
s.males <- s.males[ , c("Longitude", "Latitude", "Species")]
coordinates(s.males) <- c("Longitude", "Latitude")
proj4string(s.males) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
#kud.males <- kernelUD(s.males, h="LSCV", grid = rt.sp) # Does not converge
kud.males <- kernelUD(s.males, h = 1, grid = rt.sp)
image(kud.males)

overlap.males <- kerneloverlaphr(kud.males, method = "BA")
overlap.males
vud.males <- getvolumeUD(kud.males)
image(vud.males)

vud.males.raster.ngp <- raster(as(vud.males$'Northern giant petrel',"SpatialPixelsDataFrame"))
vud.males.raster.sgp <- raster(as(vud.males$'Southern giant petrel',"SpatialPixelsDataFrame"))
vud.males.raster <- stack(vud.males.raster.ngp, vud.males.raster.sgp)
names(vud.males.raster) <- c("Northern giant petrel", "Southern giant petrel")
plot(vud.males.raster, col = terrain.colors(100))

s.females <- s.females[ , c("Longitude", "Latitude", "Species")]
coordinates(s.females) <- c("Longitude", "Latitude")
proj4string(s.females) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
#kud.females <- kernelUD(s.females, h="LSCV", grid = rt.sp) # Does not converge
kud.females <- kernelUD(s.females, h = "href", grid = rt.sp)
image(kud.females)

overlap.females <- kerneloverlaphr(kud.females, method = "BA")
overlap.females
vud.females <- getvolumeUD(kud.females)
image(vud.females)

vud.females.raster.ngp <- raster(as(vud.females$'Northern giant petrel',"SpatialPixelsDataFrame"))
vud.females.raster.sgp <- raster(as(vud.females$'Southern giant petrel',"SpatialPixelsDataFrame"))
vud.females.raster <- stack(vud.females.raster.ngp, vud.females.raster.sgp)
names(vud.females.raster) <- c("Northern giant petrel", "Southern giant petrel")
plot(vud.females.raster, col = terrain.colors(100))

#-----------------------------------
# Extract depth and SST

# Define extent
minx <- min(tr$Longitude)
maxx <- max(tr$Longitude)
miny <- min(tr$Latitude)
maxy <- max(tr$Latitude)

# For presentation
minx <- 8
maxx <- 60
miny <- -62
maxy <- -30

# Get bathymetry
# b <- getNOAA.bathy(minx,maxx, miny, maxy, resolution = 1)
# b2 <- as.raster(b)
# writeRaster(b2, filename = "./Data/Mapping/bathymetry_large.grd", format = "raster", overwrite = T)
# b2 <- raster("./Data/Mapping/bathymetry_large.grd")

# Smaller version for large plot
# b3 <- getNOAA.bathy(minx,maxx, miny, maxy, resolution = 4)
# b4 <- as.raster(b3)
# b4[b4>0] <- NA # Set land
# writeRaster(b4, filename = "./Data/Mapping/bathymetry.grd", format = "raster", overwrite = T)

# Read in saved bathymetry
b4 <- raster("./Data/Mapping/bathymetry.grd")

b4.p <- rasterToPoints(b4)
b4.p <- data.frame(b4.p)
colnames(b4.p) <- c("lon", "lat", "val")

# Read in saved SST
sst <- raster("./Data/Mapping/ghrsst.grd")
sst.large <- aggregate(sst, 10)

sst.p <- rasterToPoints(sst.large)
sst.p <- data.frame(sst.p)
colnames(sst.p) <- c("lon", "lat", "val")

# Add the Orsifronts

myOrsi <- fortify(orsifronts)
myOrsi <- myOrsi[myOrsi$id %in% c("stf", "saf", "pf"), ]

# Create a dummy variable for plotting
dat$Species2 <- paste(dat$Species, dat$sex, sep = " - ")

# # Bathymetry
# # Coloured
# p1 <- ggplot(data = b4.p, aes(x = lon, y = lat))
# p1 <- p1 + geom_tile(aes(fill = val))
# p1 <- p1 + scale_x_continuous(expand = c(0,0))
# p1 <- p1 + scale_y_continuous(expand = c(0,0))
# p1 <- p1 + scale_fill_gradient2(low = "dodgerblue4", mid = "gainsboro", high = "burlywood4", midpoint = 0)
# p1 <- p1 + geom_point(data = dat, aes(x = dat$Longitude, y = dat$Latitude, colour = Species)) + facet_wrap(Species ~ sex) +
#   scale_colour_manual(values = c("#4daf4a", "#984ea3"))
# p1 <- p1 + coord_quickmap()
# p1 <- p1 + theme_bw() + theme(panel.grid.minor = element_blank(),
#                               panel.grid.major = element_blank())
# p1

# # Greyscale
# p2 <- ggplot(data = b4.p, aes(x = lon, y = lat))
# p2 <- p2 + geom_tile(aes(fill = val))
# # p2 <- p2 + scale_fill_gradient(low = "gray60", high = "gray99", name = "Depth (m)", limits = c(-7000, 0)) #if land is NA
# # p2 <- p2 + geom_path(data = myOrsi, aes(x = long, y = lat, group = group), linetype = 2)
# p2 <- p2 + geom_point(data = dat, aes(x = dat$Longitude, y = dat$Latitude, colour = Species), size = 0.6) + facet_wrap(Species ~ sex) +
#   scale_colour_manual(values = c("#4daf4a", "#984ea3"), guide = FALSE)
# p2 <- p2 + coord_quickmap()
# p2 <- p2 + scale_x_continuous(expand = c(0,0), limits = c(minx, maxx))
# p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(miny, maxy))
# p2 <- p2 + theme_bw() + theme(panel.grid.minor = element_blank(),
#                               panel.grid.major = element_blank(),
#                               strip.background = element_blank(),
#                               panel.border = element_rect(colour = "black"),
#                               axis.title.x = element_blank(),
#                               axis.title.y = element_blank())
# world <- borders("world", colour = "black", fill = "grey", xlim = c(min(dat$Longitude), max(dat$Longitude)), ylim = c(min(dat$Latitude), max(dat$Latitude)))
# p2 <- p2 + world #+ xlim(minx, maxx) + ylim(miny, maxy)
# p2

# With SST
p2 <- ggplot(data = sst.p, aes(x = lon, y = lat))
p2 <- p2 + geom_tile(aes(fill = val))
p2 <- p2 + scale_fill_gradientn(colours=parula(100),
                                guide = FALSE,
                                name = "SST",
                                limits = c(-2, 24))
# p2 <- p2 + geom_path(data = myOrsi, aes(x = long, y = lat, group = group), colour = "white")
p2 <- p2 + geom_point(data = dat, aes(x = dat$Longitude, y = dat$Latitude, colour = Species), size = 0.6) +
  facet_wrap(~ Species2, ncol = 2) +
  # scale_colour_manual(values = c("#4daf4a", "#984ea3"), guide = FALSE) # coloured points
  scale_colour_manual(values = c("black", "black"), guide = FALSE)
p2 <- p2 + coord_quickmap()
p2 <- p2 + labs(title = "a")
p2 <- p2 + scale_x_continuous(expand = c(0,0), limits = c(minx, maxx))
p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(miny, maxy))
p2 <- p2 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              axis.text = element_text(colour = "black"),
                              panel.border = element_rect(colour = "black"),
                              plot.title = element_text(size = rel(1)),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank())
world <- borders("world", colour = "black", fill = "grey", xlim = c(min(dat$Longitude), max(dat$Longitude)), ylim = c(min(dat$Latitude), max(dat$Latitude)))
p2 <- p2 + world #+ xlim(minx, maxx) + ylim(miny, maxy)
p2 <- p2 + annotate("point", x = 37.74, y = -46.91, size = 3, colour = "red")
# p2

#-----------------------------------
## Close-up map of island

library(rgeos)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)
library(rgdal)
library(maptools)

# Get island shapefile
island = readOGR(dsn = "./Data/Mapping", layer = "Islands_Polygonizer")
island@data$id = rownames(island@data)
island.points = fortify(island, region = "id")
island.df = join(island.points, island@data, by = "id")

# Crop to Marion only
island.df <- island.df[island.df$lat < -46.7, ]

# Quick check
ggplot(island.df) + aes(long, lat, group = group) + geom_polygon(fill = "grey") + geom_path(color = "black")

# Crop raster
# Bathymetry
b5 <- getNOAA.bathy(37.45, 38.2, -47.2, -46.8, resolution = 1)
b5 <- marmap::as.raster(b5)
b5[b5 > 0] <- NA # Set land
b5.p <- rasterToPoints(b5)
b5.p <- data.frame(b5.p)
colnames(b5.p) <- c("lon", "lat", "val")

# SST
sst.small <- crop(sst, extent(37.45, 38.2, -47.2, -46.8))
sst.small.p <- rasterToPoints(sst.small)
sst.small.p <- data.frame(sst.small.p)
colnames(sst.small.p) <- c("lon", "lat", "val")

# Only near trips
dat.short <- dat[dat$Maxdist < 100, ]

# Bathymetry
# p3 <- ggplot(data = b5.p, aes(x = lon, y = lat))
# p3 <- p3 + geom_tile(aes(fill = val))
# p3 <- p3 + scale_fill_gradient(low = "gray60", high = "gray99", name = "Depth (m)", limits = c(-7000, 0)) #if land is NA
# p3 <- p3 + geom_polygon(data = island.df, aes(x = long, y = lat, group = group), fill = "grey")
# p3 <- p3 + geom_path(data = island.df, aes(x = long, y = lat, group = group), colour = "black")
# p3 <- p3 + geom_point(data = dat.short, aes(x = dat.short$Longitude, y = dat.short$Latitude, colour = Species), size = 0.6) + facet_wrap(~Species, nrow = 2) +
#   scale_colour_manual(values = c("#4daf4a", "#984ea3"))
# p3 <- p3 + coord_quickmap()
# p3 <- p3 + scale_x_continuous(expand = c(0,0), limits = c(37.45, 38.2))
# p3 <- p3 + scale_y_continuous(expand = c(0,0), limits = c(-47.2, -46.8))
# p3 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),
#                               panel.grid.major = element_blank(),
#                               strip.background = element_blank(),
#                               panel.border = element_rect(colour = "black"),
#                               axis.title.x = element_blank(),
#                               axis.title.y = element_blank())
# p3

# SST
p3 <- ggplot(data = sst.small.p, aes(x = lon, y = lat))
p3 <- p3 + geom_tile(aes(fill = val))
p3 <- p3 + scale_fill_gradientn(colours=parula(100),
                                guide = "colourbar",
                                name = "SST",
                                limits = c(-2, 24))
p3 <- p3 + geom_polygon(data = island.df, aes(x = long, y = lat, group = group), fill = "grey")
p3 <- p3 + geom_path(data = island.df, aes(x = long, y = lat, group = group), colour = "black")
p3 <- p3 + geom_point(data = dat.short, aes(x = dat.short$Longitude, y = dat.short$Latitude, colour = Species), size = 0.6) +
  facet_wrap(~Species, nrow = 2) +
  # scale_colour_manual(values = c("#4daf4a", "#984ea3"), guide = FALSE)
  scale_colour_manual(values = c("black", "black"), guide = FALSE)
p3 <- p3 + coord_quickmap()
p3 <- p3 + labs(title = "b")
p3 <- p3 + scale_x_continuous(expand = c(0,0), limits = c(37.45, 38.2))
p3 <- p3 + scale_y_continuous(expand = c(0,0), limits = c(-47.2, -46.795))
p3 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              axis.text = element_text(colour = "black"),
                              panel.border = element_rect(colour = "black"),
                              plot.title = element_text(size = rel(1)),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank())
# p3

# grid.arrange(p2, p3, ncol = 2)

# Plot to file

# Distant trips
pdf("./Plots/mapDistant.pdf",
    width = double.col/fig.scale, height = (double.col/1.6)/fig.scale,
    useDingbats = FALSE)
print(p2)
dev.off()

# Nearby trips
pdf("./Plots/mapNearby.pdf",
    width = double.col/fig.scale, height = (double.col/1.6)/fig.scale,
    useDingbats = FALSE)
print(p3)
dev.off()

# Combined
pdf("./Plots/mapCombined.pdf",
    width = double.col/fig.scale, height = double.col/fig.scale,
    useDingbats = FALSE)
cowplot::plot_grid(p2, p3, axis = "t", ncol = 2, rel_widths = c(2, 1.5))
dev.off()

tiff("./Plots/mapCombined.tiff",
    width = double.col/fig.scale, height = double.col/fig.scale,
    res = 600,
    units = "in")
cowplot::plot_grid(p2, p3, axis = "t", ncol = 2, rel_widths = c(2, 1.5))
dev.off()

# Remove dummy variable
dat$Species2 <- NULL

#-----------------------------------------------------
## Extract depth

tr3$DEPTH <- extract(x = b2, y = tr3[ , c("Longitude", "Latitude")])

# ## Extract SST
# library(xtractomatic)
# 
# # Hack to extract SST, connection seems to time out for too many lines
# tr.sub1 <- tr3[1:1000, ]
# tr.sub2 <- tr3[1001:2000,]
# tr.sub3 <- tr3[2001:3000,]
# tr.sub4 <- tr3[3001:4000,]
# tr.sub5 <- tr3[4001:5000,]
# tr.sub6 <- tr3[5001:6000,]
# tr.sub7 <- tr3[6001:nrow(tr3),]
# 
# SST1 <- xtracto(xpos = tr.sub1$Longitude, ypos = tr.sub1$Latitude, tpos = tr.sub1$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# SST2 <- xtracto(xpos = tr.sub2$Longitude, ypos = tr.sub2$Latitude, tpos = tr.sub2$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# SST3 <- xtracto(xpos = tr.sub3$Longitude, ypos = tr.sub3$Latitude, tpos = tr.sub3$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# SST4 <- xtracto(xpos = tr.sub4$Longitude, ypos = tr.sub4$Latitude, tpos = tr.sub4$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# SST5 <- xtracto(xpos = tr.sub5$Longitude, ypos = tr.sub5$Latitude, tpos = tr.sub5$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# SST6 <- xtracto(xpos = tr.sub6$Longitude, ypos = tr.sub6$Latitude, tpos = tr.sub6$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# SST7 <- xtracto(xpos = tr.sub7$Longitude, ypos = tr.sub7$Latitude, tpos = tr.sub7$date.time,
#                 xlen = .01, ylen = .01, dtype = "jplG1SST")
# 
# SST <- rbind.data.frame(SST1, SST2, SST3, SST4, SST5, SST6, SST7)
# 
# tr3$SST <- SST$mean
# tr3[tr3$SST == "NaN", "SST"] <- NA
# 
# # Save
# saveRDS(tr3, "./Output/tr_03-2.RDS")

# Get pre-run
tr3 <- readRDS("./Output/tr_03-2.RDS")


#-----------------------------------
# Continue extractions on  R Server
# using script giant_petrels_extract_envars.R

tr3_env <- readRDS("./Output/tr_03-2_env.RDS")

# Check temp is in C:
if (max(tr3_env$SST) > 100) {
  tr3_env$SST <- tr3_env$SST-273.15
}

#-----------------------------------
# Visualization

# Bill dimensions
# Tweak sex.m label for display
met[is.na(met$sex.m), "sex.m"] <- "Not tested"

p3 <- ggplot(data = met, aes(x = Culmen_Length, y = Culmen_Depth, colour = Species, shape = sex.m))
p3 <- p3 + geom_point(size = 3)
p3 <- p3 + scale_colour_manual(values = c("#4daf4a", "#984ea3"),
                               name = "Species")
p3 <- p3 + scale_shape_manual(values = c(15, 16, 17),
                              name = "Sex\n(molecular testing)")
p3 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text = element_text(colour = "black"))
p3 <- p3 + labs(x = "Culmen length (mm)", y = "Bill depth (mm)")
p3 <- p3 + geom_vline(xintercept = 96, lty = 2)
# p3

pdf("./Plots/billDimensions.pdf",
    width = single.col/fig.scale, height = (single.col/1.6)/fig.scale,
    useDingbats = FALSE)
print(p3)
dev.off()

ngp.culmen <- met[met$Species == "Northern giant petrel", "Culmen_Length"]
sgp.culmen <- met[met$Species == "Southern giant petrel", "Culmen_Length"]
wilcox.test(ngp.culmen, sgp.culmen)

ngp.depth <- met[met$Species == "Northern giant petrel", "Culmen_Depth"]
sgp.depth <- met[met$Species == "Southern giant petrel", "Culmen_Depth"]
wilcox.test(ngp.depth, sgp.depth)

# Lat & Lon
library(beanplot)

# Latitude
beanplot(Latitude ~ Species, data = tr3[tr$sex == "Male", ],
         main = NULL, ylab = "Latitude", side = "both",
         ylim = c(min(tr3[tr3$sex == "Male", "Latitude"], na.rm = T), max(tr3[tr3$sex == "Male", "Latitude"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

beanplot(Latitude ~ Species, data = tr3[tr3$sex == "Female", ],
         main = NULL, ylab = "Latitude", side = "both",
         ylim = c(min(tr3[tr3$sex == "Female", "Latitude"]), max(tr3[tr3$sex == "Female", "Latitude"])),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

# SST
beanplot(SST ~ Species, data = tr3[tr3$sex == "Male", ],
         main = NULL, ylab = "SST", side = "both",
         ylim = c(min(tr3[tr3$sex == "Male", "SST"], na.rm = T), max(tr3[tr3$sex == "Male", "SST"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

beanplot(SST ~ Species, data = tr3[tr3$sex == "Female", ],
         main = NULL, ylab = "SST", side = "both",
         ylim = c(min(tr3[tr3$sex == "Female", "SST"], na.rm = T), max(tr3[tr3$sex == "Female", "SST"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

# Depth
beanplot(DEPTH ~ Species, data = tr3[tr3$sex == "Male", ],
         main = NULL, ylab = "DEPTH", side = "both",
         ylim = c(min(tr3[tr3$sex == "Male", "DEPTH"], na.rm = T), max(tr3[tr3$sex == "Male", "DEPTH"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

beanplot(DEPTH ~ Species, data = tr3[tr3$sex == "Female", ],
         main = NULL, ylab = "DEPTH", side = "both",
         ylim = c(min(tr3[tr3$sex == "Female", "DEPTH"], na.rm = T), max(tr3[tr3$sex == "Female", "DEPTH"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

#-----------------------------------
# Try split violins with GGPlot
source("./Scripts/function_splitViolin.R")
library(tidyr)

datV <- tr3_env
datV$CHL <- log10(datV$CHL)

datV <- gather(data = datV,
              key = "Covariate",
              value = "Value",
              "Latitude",
              "SST",
              "CHL",
              "DEPTH2")

datV[datV$Covariate == "DEPTH2", "Covariate"] <- "DEPTH"

# Get in the correct order
datV$Covariate <- factor(datV$Covariate, levels = c("Latitude", "SST", "DEPTH", "CHL"))

pdf("./Plots/envarDensity.pdf",
    width = (single.col*1.6/fig.scale),
    height = (single.col*1.6/fig.scale),
    useDingbats = FALSE)
beans <- ggplot(data = datV, aes(x = sex, y = Value, fill = Species)) +
  geom_split_violin() +
  scale_fill_manual(values = c("#4daf4a", "#984ea3"), name = "Species") +
  facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Sex", y = "Covariate value") +
  theme_rr()
print(beans)
dev.off()

#-----------------------------------
# Summary table

summary <- met[ , c("Species", "ID", "Culmen_Length", "sex", "Deployment.date", "Deployment.time..start.","Retrieval.date", "Retrieval.time")]
summary <- unique(summary)
tmp <- tr3[ , c("ID", "Maxdist")]
tmp <- unique(tmp)
summary <- merge(x = summary, y = tmp, by = "ID", all.y = F)
summary <- summary[ , c("Species", "ID", "Culmen_Length", "sex", "Deployment.date", "Deployment.time..start.","Retrieval.date", "Retrieval.time", "Maxdist")]

summary$date.start <- paste0(summary$Deployment.date, " ", summary$Deployment.time..start.)
summary$date.start <- strptime(summary$date.start, format = "%d-%m-%Y %H:%M", tz = "GMT")

summary$date.end <- paste0(summary$Retrieval.date, " ", summary$Retrieval.time)
summary$date.end <- strptime(summary$date.end, format = "%d-%m-%Y %H:%M", tz = "GMT")

summary$duration <- rep(NA, nrow(summary))
for (i in 1:nrow(summary)) {
  t <- difftime(summary$date.end[i], summary$date.start[i], units = "days")
  summary$duration[i] <- t
}

summary$duration <- round(summary$duration, 1)
summary$Maxdist <- round(summary$Maxdist, 1)

write.csv(summary, file = "./Output/SummaryTable.csv", row.names = F)

# Read back in if neccessary
# summary <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)

summary$date.start <- as.POSIXct(summary$date.start, format = "%Y-%m-%d %H:%M:%S")
summary$date.end <- as.POSIXct(summary$date.end, format = "%Y-%m-%d %H:%M:%S")

pdf("./Plots/trackingPeriods.pdf",
    width = double.col.sup/fig.scale,
    height = double.col.sup/fig.scale,
    useDingbats = FALSE)
p1 <- ggplot(data = summary) +
  geom_linerange(aes(x = ID, ymin = date.start, ymax = date.end, colour = Species)) +
  scale_y_datetime(date_breaks = "1 week") +
  scale_colour_manual(values = c("#4daf4a", "#984ea3")) +
  coord_flip() +
  theme_rr() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
print(p1)
dev.off()

#-----------------------------------
## Random forests

library(randomForest)

dat.f <- tr3_env[ , c("ID", "Species", "sex", "Distance", "DEPTH2", "SST", "DIST.200", "SSTgrad",
                      "SSHA", "SSHgrad", "EKE", "MLD", "CHL", "PROD", "WINDU", "WINDV")]
# Select only complete cases
dat.f <- dat.f[complete.cases(dat.f), ]

# Log Chlorophyll
dat.f$CHL <- log10(dat.f$CHL)

# Only at-sea
dat.f <- dat.f[dat.f$DEPTH2 <= 0, ]

# First check for colinearity
cor(x = dat.f[ , c(4:16)], method = "spearman")

library(corrplot)
cr <- cor(x = dat.f[ , c(4:16)], method = "pearson")
corrplot(cr, method = "number", type = "lower")

library(GGally)
ggpairs(dat.f[ , c(4:16)])

# Leave out Distance, DIST.200, PROD, EKE

# Add species abreviations
dat.f$Species.code <- rep(NA, nrow(dat.f))
dat.f[dat.f$Species == "Northern giant petrel", "Species.code"] <- "NGP"
dat.f[dat.f$Species == "Southern giant petrel", "Species.code"] <- "SGP"

dat.f$Group <- rep(NA, nrow(dat.f))
for (i in 1:nrow(dat.f)) {
  dat.f$Group[i] <- paste0(dat.f$Species.code[i], " ", dat.f$sex[i])
}

# Build RF
rf <- randomForest(data = dat.f, as.factor(Group) ~ DEPTH2 + SST + SSTgrad +
                     SSHA + SSHgrad + CHL + WINDU + WINDV + MLD, ntree = 1000,
                   na.action = na.omit, proximity = T, importance = F)
rf

varImpPlot(rf, pch = 16)

# Plots
pdf("./Plots/variableImportance.pdf",
    width = single.col/fig.scale,
    height = single.col/fig.scale,
    useDingbats = FALSE)
# par(ps = 9) #set font size
varImpPlot(rf, pch = 16)
dev.off()


plot(dat.f$SST, dat.f$DEPTH2, col = as.factor(dat.f$Group)) # Plot two most important vars



pdf("./Plots/envars2D.pdf",
    width = single.col/fig.scale,
    height = single.col/fig.scale,
    useDingbats = FALSE)
p <- ggplot(dat.f, aes(x = SST, y = DEPTH2, col = Group)) +
  geom_point() +
  scale_colour_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a")) +
  theme_rr() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
print(p)
dev.off()

MDSplot(rf, fac = as.factor(dat.f$Group), k = 5)

# MDS with rfPermute
library(rfPermute)
proximityPlot(rf, legend.loc = "right", circle.size = NULL)


# Make own plot - proximity as MDS
rf.mds <- cmdscale(1 - rf$proximity, k = 2)
mdsframe <- cbind.data.frame(dat.f$Group, rf.mds[ , 1], rf.mds[ ,2])
names(mdsframe) <- c("Group", "Dim1", "Dim2")

p.mds <- ggplot(mdsframe, aes(x = Dim1, y = Dim2, colour = Group))
p.mds <- p.mds + geom_point(size = 1.5, alpha = 1) + scale_colour_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"))
p.mds <- p.mds + theme_bw() + theme(panel.grid.minor = element_blank(),
                                    panel.grid.major = element_blank())
p.mds <- p.mds + labs(x = "Dimension 2", y = "Dimension 1")
p.mds

# 3D plots
library(plot3D)
Groups <- unique(dat.f$Group)
dat.f$g <- rep(NA, nrow(dat.f))
dat.f[dat.f$Group == Groups[1], "g"] <- 1
dat.f[dat.f$Group == Groups[2], "g"] <- 2
dat.f[dat.f$Group == Groups[3], "g"] <- 3
dat.f[dat.f$Group == Groups[4], "g"] <- 4

scatter3D(x = dat.f$DEPTH2, y = dat.f$SST, z = dat.f$CHL, colvar = dat.f$g,
          col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"),
          bty = "f",
          pch = 20,
          phi = 30, theta = 50,
          ticktype = "detailed",
          xlab = "DEPTH", ylab = "SST", zlab = "log(CHL)",
          colkey = list(at = c(1,2,3,4),labels = c("NGP Male", "NGP Female", "SGP Male", "SGP Female")))