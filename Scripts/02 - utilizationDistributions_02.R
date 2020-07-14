# Giant Petrels - Utilization distributions

## Ryan Reisinger

library(adehabitatHR)
library(raster)
library(viridis)
library(ggplot2)
library(orsifronts)

library(rgeos)
library(rgdal)
library(plyr)

library(mnormt)
library(EMbC)

library(classInt)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

# --------------------------------
## Plotting stuff

## Figure widths in mm
single.col <- 84*0.0393701
double.col <- 117*0.0393701
double.col.sup <- 150*0.0393701

## Scaling for font size
fig.scale <- 8/11

## Custom theme for plots
source("./Scripts/theme_rr.R")

# --------------------------------
## Get data
dat <- readRDS("./Output/tracks_trips.RDS")

#-----------------------------------
## Utilization distributions

## Drop missing coordinates
dat <- dat[!is.na(dat$decimal_longitude) & !is.na(dat$decimal_latitude), ]

## Expectation maximization binary clustering
ids <- unique(dat$trip_id)
expth <- list()
for (i in ids) {
this_dat <- dat[dat$trip_id == i, c("date.time", "decimal_longitude", "decimal_latitude")]
expth <- append(expth, list(this_dat))
rm(this_dat)
}

mybcp <- stbc(expth, info=-1)

# Look at selected individuals
sctr(slct(mybcp, 100))
view(slct(mybcp, 100))

# Get output
dat$embc_velocity <- mybcp@bC@X[,1]
dat$embc_turnrad <- mybcp@bC@X[,2]
dat$embc_label <- mybcp@bC@A

# Check the labels
ggplot(data = dat, aes(x = embc_velocity, y = embc_turnrad, colour = as.factor(embc_label), shape = sex)) +
  geom_point()

# Plot spatially
ggplot(data = dat, aes(x = decimal_longitude, y = decimal_latitude, colour = embc_velocity)) +
  geom_point() +
  scale_color_viridis_c()

# Look at velocity distribution
ggplot(data = dat, aes(x = embc_velocity, y = sex)) +
  geom_boxplot()

# Fisher-Jenks breaks?
# The result is very similar to k-means
v_cut <- classInt::classIntervals(dat$embc_velocity, n = 2, style = "fisher")$brks[2]

# Quantile threshold?
if (FALSE) {
v_cut <- quantile(dat$embc_velocity, probs = c(0.75))
}

# Group 3 = high velocity & low turning angle - shouldn't be used in KUD analyses
dat$mode <- "transit"
dat[dat$embc_velocity < v_cut, "mode"] <- "restricted"

# Plot spatially again
ggplot(data = dat, aes(x = decimal_longitude, y = decimal_latitude, colour = mode)) +
  geom_point() +
  scale_color_viridis_d()

s <- dat[dat$mode == "restricted", ]
s.males <- dat[dat$sex == "Male", ]
s.females <- dat[dat$sex == "Female", ]

# Short trips only
s.short <- dat[dat$trip.Maxdist < 50, ]

# Grid for UD
lms <- c(min(s$decimal_longitude, na.rm = T) - 5,
         max(s$decimal_longitude, na.rm = T) + 5,
         min(s$decimal_latitude, na.rm = T) - 5,
         max(s$decimal_latitude, na.rm = T) + 5)
rt <- raster(ext = extent(lms), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
             res = 0.05)
rt.sp <- as(rt, "SpatialPixelsDataFrame")

# Small grid
lms.short <- c(min(s.short$decimal_longitude, na.rm = T) - 0.05,
         max(s.short$decimal_longitude, na.rm = T) + 0.05,
         min(s.short$decimal_latitude, na.rm = T) - 0.05,
         max(s.short$decimal_latitude, na.rm = T) + 0.05)
rt.short <- raster(ext = extent(lms.short),
                   crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                   res = 0.0005)
rt.sp.short <- as(rt.short, "SpatialPixelsDataFrame")

# Group sex x species
s$group <- paste0(s$scientific_name, " ", s$sex)
s <- s[ , c("decimal_longitude", "decimal_latitude", "group")]
coordinates(s) <- c("decimal_longitude", "decimal_latitude")
proj4string(s) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud <- kernelUD(s, h = 0.8, grid = rt.sp)
image(kud)

# KUD of short trips only
s.short <- s.short[ , c("decimal_longitude", "decimal_latitude", "scientific_name")]
coordinates(s.short) <- c("decimal_longitude", "decimal_latitude")
proj4string(s.short) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.short <- kernelUD(s.short, h = 0.005, grid = rt.sp.short)
image(kud.short)

# Create rasters of each
vud <- getvolumeUD(kud)
vud.males.raster.ngp <- raster(as(vud$'Northern Giant Petrel Male',"SpatialPixelsDataFrame"))
vud.males.raster.sgp <- raster(as(vud$'Southern Giant Petrel Male',"SpatialPixelsDataFrame"))
vud.females.raster.ngp <- raster(as(vud$'Northern Giant Petrel Female',"SpatialPixelsDataFrame"))
vud.females.raster.sgp <- raster(as(vud$'Southern Giant Petrel Female',"SpatialPixelsDataFrame"))

vud.short <- getvolumeUD(kud.short)
vud.short.raster.ngp <- raster(as(vud.short$'Northern Giant Petrel',"SpatialPixelsDataFrame"))
vud.short.raster.sgp <- raster(as(vud.short$'Southern Giant Petrel',"SpatialPixelsDataFrame"))

# Get contours
all.contour <- getverticeshr(kud, percent = 95)
short.contour <- getverticeshr(kud.short, percent = 95)

all.contour.50 <- getverticeshr(kud, percent = 50)
short.contour.50 <- getverticeshr(kud.short, percent = 50)

# Write the rasters
writeRaster(vud.males.raster.ngp, "./Output/vud.males.raster.ngp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.males.raster.sgp, "./Output/vud.males.raster.sgp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.females.raster.ngp, "./Output/vud.females.raster.ngp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.females.raster.sgp, "./Output/vud.females.raster.sgp.grd", format = "raster",
            overwrite = T)

writeRaster(vud.short.raster.sgp, "./Output/vud.short.raster.sgp.grd", format = "raster",
            overwrite = T)
writeRaster(vud.short.raster.ngp, "./Output/vud.short.raster.ngp.grd", format = "raster",
            overwrite = T)


# By sex
s.males <- s.males[ , c("decimal_longitude", "decimal_latitude", "scientific_name")]
coordinates(s.males) <- c("decimal_longitude", "decimal_latitude")
proj4string(s.males) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.males <- kernelUD(s.males, h = 0.8, grid = rt.sp)
image(kud.males)

overlap.males <- kerneloverlaphr(kud.males, method = "BA")
overlap.males
vud.males <- getvolumeUD(kud.males)
image(vud.males)

vud.males.raster.ngp <- raster(as(vud.males$'Northern Giant Petrel',"SpatialPixelsDataFrame"))
vud.males.raster.sgp <- raster(as(vud.males$'Southern Giant Petrel',"SpatialPixelsDataFrame"))
vud.males.raster <- stack(vud.males.raster.ngp, vud.males.raster.sgp)
names(vud.males.raster) <- c("Northern Giant Petrel", "Southern Giant Petrel")
plot(vud.males.raster, col = terrain.colors(100))

s.females <- s.females[ , c("decimal_longitude", "decimal_latitude", "scientific_name")]
coordinates(s.females) <- c("decimal_longitude", "decimal_latitude")
proj4string(s.females) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
#kud.females <- kernelUD(s.females, h="LSCV", grid = rt.sp) # Does not converge
kud.females <- kernelUD(s.females, h = 0.8, grid = rt.sp)
image(kud.females)

overlap.females <- kerneloverlaphr(kud.females, method = "BA")
overlap.females
vud.females <- getvolumeUD(kud.females)
image(vud.females)

vud.females.raster.ngp <- raster(as(vud.females$'Northern Giant Petrel',"SpatialPixelsDataFrame"))
vud.females.raster.sgp <- raster(as(vud.females$'Southern Giant Petrel',"SpatialPixelsDataFrame"))
vud.females.raster <- stack(vud.females.raster.ngp, vud.females.raster.sgp)
names(vud.females.raster) <- c("Northern Giant Petrel", "Southern Giant Petrel")
plot(vud.females.raster, col = terrain.colors(100))

#### Maps of giant petrel utilization distributions

vud.males.ngp <- raster("./Output/vud.males.raster.ngp.grd")
vud.males.sgp <- raster("./Output/vud.males.raster.sgp.grd")
vud.females.ngp <- raster("./Output/vud.females.raster.ngp.grd")
vud.females.sgp <- raster("./Output/vud.females.raster.sgp.grd")

vud.short.ngp <- raster("./Output/vud.short.raster.ngp.grd")
vud.short.sgp <- raster("./Output/vud.short.raster.sgp.grd")

# NGP Males
ngpm <- rasterToPoints(vud.males.ngp)
ngpm <- data.frame(ngpm)
colnames(ngpm) <- c("lon", "lat", "val")
ngpm$sp <- "Northern giant petrel"
ngpm$sex <- "Male"

# NGP Females
ngpf <- rasterToPoints(vud.females.ngp)
ngpf <- data.frame(ngpf)
colnames(ngpf) <- c("lon", "lat", "val")
ngpf$sp <- "Northern giant petrel"
ngpf$sex <- "Female"

# SGP Males
sgpm <- rasterToPoints(vud.males.sgp)
sgpm <- data.frame(sgpm)
colnames(sgpm) <- c("lon", "lat", "val")
sgpm$sp <- "Southern giant petrel"
sgpm$sex <- "Male"

# SGP Females
sgpf <- rasterToPoints(vud.females.sgp)
sgpf <- data.frame(sgpf)
colnames(sgpf) <- c("lon", "lat", "val")
sgpf$sp <- "Southern giant petrel"
sgpf$sex <- "Female"

# Bind
kern <- rbind(ngpf, ngpm, sgpf, sgpm)
kern <- kern[kern$val <= 95, ]

## Short trips
# SGP
sgps <- rasterToPoints(vud.short.sgp)
sgps <- data.frame(sgps)
colnames(sgps) <- c("lon", "lat", "val")
sgps$sp <- "Southern giant petrel"
# NGP
ngps <- rasterToPoints(vud.short.ngp)
ngps <- data.frame(ngps)
colnames(ngps) <- c("lon", "lat", "val")
ngps$sp <- "Northern giant petrel"

kern.short <- rbind(sgps, ngps)
kern.short <- kern.short[kern.short$val <= 95, ]

# Get Southern Ocean fronts from Park & Durand 2019
# https://doi.org/10.17882/59800

library(ncdf4)
frnts <- nc_open("./Data/Mapping/62985.nc")
NB <- data.frame(
  "lat" = ncvar_get(frnts, "LatNB"),
  "lon" = ncvar_get(frnts, "LonNB"),
  "name" = "NB"
)
SAF <- data.frame(
  "lat" = ncvar_get(frnts, "LatSAF"),
  "lon" = ncvar_get(frnts, "LonSAF"),
  "name" = "SAF"
)
PF <- data.frame(
  "lat" = ncvar_get(frnts, "LatPF"),
  "lon" = ncvar_get(frnts, "LonPF"),
  "name" = "PF"
)
SACCF <- data.frame(
  "lat" = ncvar_get(frnts, "LatSACCF"),
  "lon" = ncvar_get(frnts, "LonSACCF"),
  "name" = "SACCF"
)
SB <- data.frame(
  "lat" = ncvar_get(frnts, "LatSB"),
  "lon" = ncvar_get(frnts, "LonSB"),
  "name" = "SB"
)
nc_close(frnts)

frnts <- rbind(NB, SAF, PF, SACCF, SB)

## Create contours df
# 95
# All
contourdf_all <- fortify(all.contour)
contourdf_all$sp <- NA
contourdf_all$sex <- NA
contourdf_all[grep("Northern Giant Petrel", contourdf_all$id), "sp"] <- "Northern giant petrel"
contourdf_all[grep("Southern Giant Petrel", contourdf_all$id), "sp"] <- "Southern giant petrel"
contourdf_all[grep("Male", contourdf_all$id), "sex"] <- "Male"
contourdf_all[grep("Female", contourdf_all$id), "sex"] <- "Female"

# Short
contourdf_short <- fortify(short.contour)
contourdf_short$sp <- NA
contourdf_short$sex <- NA
contourdf_short[grep("Northern Giant Petrel", contourdf_short$id), "sp"] <- "Northern giant petrel"
contourdf_short[grep("Southern Giant Petrel", contourdf_short$id), "sp"] <- "Southern giant petrel"
contourdf_short[grep("Male", contourdf_short$id), "sex"] <- "Male"
contourdf_short[grep("Female", contourdf_short$id), "sex"] <- "Female"

# 50
# All
contourdf_all_50 <- fortify(all.contour.50)
contourdf_all_50$sp <- NA
contourdf_all_50$sex <- NA
contourdf_all_50[grep("Northern Giant Petrel", contourdf_all_50$id), "sp"] <- "Northern giant petrel"
contourdf_all_50[grep("Southern Giant Petrel", contourdf_all_50$id), "sp"] <- "Southern giant petrel"
contourdf_all_50[grep("Male", contourdf_all_50$id), "sex"] <- "Male"
contourdf_all_50[grep("Female", contourdf_all_50$id), "sex"] <- "Female"

# Short
contourdf_short_50 <- fortify(short.contour.50)
contourdf_short_50$sp <- NA
contourdf_short_50$sex <- NA
contourdf_short_50[grep("Northern Giant Petrel", contourdf_short_50$id), "sp"] <- "Northern giant petrel"
contourdf_short_50[grep("Southern Giant Petrel", contourdf_short_50$id), "sp"] <- "Southern giant petrel"
contourdf_short_50[grep("Male", contourdf_short_50$id), "sex"] <- "Male"
contourdf_short_50[grep("Female", contourdf_short_50$id), "sex"] <- "Female"

# Map
minx <- 8
maxx <- 60
miny <- -64
maxy <- -30

# Plot
p2 <- ggplot(data = kern, aes(x = lon, y = lat))
p2 <- p2 + geom_tile(aes(fill = val)) + facet_grid(sp ~ sex)
p2 <- p2 + scale_fill_viridis(direction = 1, option = "plasma",
                              name = "Utilisation \ndistribution\n(%)")
# p2 <- p2 + geom_path(data = myOrsi, aes(x = long, y = lat, group = group), linetype = 2)
p2 <- p2 + geom_path(data = frnts, aes(x = lon, y = lat, group = name), linetype = 1, colour = "grey60")
p2 <- p2 + geom_path(data = contourdf_all, aes(x = long, y = lat, group = group), linetype = 1, colour = "black")
p2 <- p2 + geom_path(data = contourdf_all_50, aes(x = long, y = lat, group = group), linetype = 1, colour = "black")
p2 <- p2 + coord_quickmap()
p2 <- p2 + scale_x_continuous(expand = c(0,0), limits = c(minx, maxx))
p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(miny, maxy))
p2 <- p2 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank())
world <- borders("world", colour = "black", fill = "grey", xlim = c(minx, maxx), ylim = c(miny, maxy))
p2 <- p2 +
  world +
  theme_rr() +
  labs(x = "Longitude", y = "Latitude") +
  annotate(geom = "point",
           x = mean(dat$deployment_decimal_longitude),
           y = mean(dat$deployment_decimal_latitude),
           shape = 18,
           size = 2.5,
           colour = "white") #+
  # annotate(geom = "text",
  #          x = 22,
  #          y = -32,
  #          label = "SOUTH\nAFRICA",
  #          colour = "black") +

pdf("./Plots/utilizationDistributions.pdf",
    width = double.col/fig.scale,
    height = double.col/fig.scale,
    useDingbats = FALSE)
print(p2)
dev.off()


# Plot of short trips only...

## New limits
minx <- 37.55
maxx <- 38.00
maxy <- -46.80
miny <- -47.05

## Get island shape
# Get island shapefile
island = readOGR(dsn = "./Data/Mapping", layer = "Islands_Polygonizer")
island@data$id = rownames(island@data)
island.points = fortify(island, region = "id")
island.df = plyr::join(island.points, island@data, by = "id")

# Crop to Marion only
island.df <- island.df[island.df$lat < -46.7, ]

# Create a copy of deployment locations
deploy <- dat[, c("track_id", "sp_code", "deployment_decimal_latitude",
                  "deployment_decimal_longitude")]
deploy <- deploy[!duplicated(deploy$track_id), ]

deploy.ngp <- deploy[deploy$sp_code == "NGP", ]
deploy.ngp$sp <- factor("Northern giant petrel",levels = c("Northern giant petrel",
                                                           "Southern giant petrel"))
deploy.sgp <- deploy[deploy$sp_code == "SGP", ]
deploy.sgp$sp <- factor("Southern giant petrel",levels = c("Northern giant petrel",
                                                           "Southern giant petrel"))

# # And the UDs
# kern.short <- kern.short[kern.short$lon > minx & kern.short$lon < maxx, ]
# kern.short <- kern.short[kern.short$lat > miny & kern.short$lat < maxy, ]

p3 <- ggplot() +
  geom_polygon(data = island.df, aes(x = long, y = lat, group = group), fill = "grey") +
  geom_tile(data = kern.short, aes(x = lon, y = lat, fill = val)) +
  geom_path(data = island.df, aes(x = long, y = lat, group = group), colour = "black") +
  # geom_point(data = island.df, aes(x = long, y = lat, group = group), colour = "black") +
  geom_path(data = contourdf_short, aes(x = long, y = lat, group = group), linetype = 1, colour = "black") +
  geom_path(data = contourdf_short_50, aes(x = long, y = lat, group = group), linetype = 1, colour = "black") +
  facet_grid(sp ~ .) +
  scale_fill_viridis(direction = 1, option = "plasma",
                              name = "Utilisation \ndistribution\n(%)") +
  coord_quickmap() +
  scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
  scale_y_continuous(expand = c(0,0), limits = c(miny, maxy)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank()) +
  theme_rr() +
  labs(x = "Longitude", y = "Latitude") +
  geom_point(data = deploy.ngp, aes(x = deployment_decimal_longitude,
                                    y = deployment_decimal_latitude,
                                    group = sp),
             shape = 18, size = 2.5, colour = "white") +
  geom_point(data = deploy.sgp, aes(x = deployment_decimal_longitude,
                                    y = deployment_decimal_latitude,
                                    group = sp),
             shape = 18, size = 2.5, colour = "white") +
  theme(legend.position = "bottom")

pdf("./Plots/utilizationDistributionsShort.pdf",
    width = single.col/fig.scale,
    height = 4/fig.scale,
    useDingbats = FALSE)
print(p3)
dev.off()
