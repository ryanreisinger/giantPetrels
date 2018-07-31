#### Permutation tests for distribution overlap

## Ryan Reisinger

# The script assumes that 'giant_petrels.R' has been run previously,
# or the output of that script is present in './Output/'.

library(ggplot2)
library(argosfilter)
library(ggmap)
library(geosphere)
library(rgdal)
library(raster)
library(marmap)
# library(rgeos)
# library(maptools)
library(adehabitatHR)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

tr3 <- readRDS("./Output/tr_03.RDS")
met <- readRDS("./Output/met.RDS")

# Make sure sexes are updated if reading in RDS
tr3[tr3$ID == "NGP19_KD_SEP_2015", "sex"] <- "Male"
tr3[tr3$ID == "SGP02_KD_SEP_2015", "sex"] <- "Male"

met[met$ID == "NGP19_KD_SEP_2015", "sex"] <- "Male"
met[met$ID == "SGP02_KD_SEP_2015", "sex"] <- "Male"

#--------------------------------------------------------
ref <- met[ , c("ID", "Species", "sex")]

# Count NGPS and Males (will use them as the reference levels)
NGPS <- length(ref$Species[ref$Species == "Northern giant petrel"])
Males <- length(ref$sex[ref$sex == "Male"])

#---------------------------------------------------
# 1. Swap species labels

# Set up the frame
nperm = 1000
comb <- matrix(nrow = nrow(ref), ncol = nperm)
comb <- as.data.frame(comb)
comb <- cbind.data.frame(ref, comb)

for (i in 1:1000) {
comb[ , i + 3] <- rep("Southern giant petrel", nrow(comb))
one <- sample(x = nrow(ref), size = NGPS, replace = F)
comb[one, i + 3] <- "Northern giant petrel"
}

# Run the tests based on the permutation frame
# First, setup

lms <- c(min(tr3$Longitude) - 5, max(tr3$Longitude) + 5, min(tr3$Latitude) - 5, max(tr3$Latitude) + 5)
rt <- raster(ext = extent(lms), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), res = 0.05)
rt.sp <- as(rt, "SpatialPixelsDataFrame")

s1 <- tr3
s1 <- s1[s1$sex != "Unknown", ]
s1 <- s1[s1$Distance > 0.4, ]
s1.males <- s1[s1$sex == "Male", ]
s1.females <- s1[s1$sex == "Female", ]

##########
# Males
# Real value
s1.males.sp <- s1.males[ , c("Longitude", "Latitude", "Species")]
coordinates(s1.males.sp) <- c("Longitude", "Latitude")
proj4string(s1.males.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.males <- kernelUD(s1.males.sp, h = 1.0, grid = rt.sp)
overlap.males <- kerneloverlaphr(kud.males, method = "BA", percent = 95, conditional = T)
o1 <- overlap.males[2, 1]
o2 <- overlap.males[1, 2]
overlap.males.real <- mean(o1, o2)
overlap.males.real

# Plot hr
# vud.males <- getvolumeUD(kud.males)
# vud.males.raster.ngp <- raster(as(vud.males$'Northern giant petrel',"SpatialPixelsDataFrame"))
# vud.males.raster.sgp <- raster(as(vud.males$'Southern giant petrel',"SpatialPixelsDataFrame"))
# vud.males.raster <- stack(vud.males.raster.ngp, vud.males.raster.sgp)
# names(vud.males.raster) <- c("Northern giant petrel", "Southern giant petrel")
# plot(vud.males.raster, col = terrain.colors(100))
# writeRaster(vud.males.raster, file="vud.males.raster.grd", format = "raster")

# Permute
md <- s1.males[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
m <- comb[ , c(1, i + 3)]
names(m) <- c("ID", "perm.sp")
ml <- merge(x = m, y = md, by = "ID", all = F)
ml$Species <- NULL
ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sp")]
coordinates(ml.sp) <- c("Longitude", "Latitude")
proj4string(ml.sp) <- projt
kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = 95, conditional = T)
ov1 <- overlap.ml[2, 1]
ov2 <- overlap.ml[1, 2]
ovr[i] <- mean(ov1, ov2)
}

ovr.males <- ovr

hist(ovr.males, breaks = seq(0, 1, 0.01))
abline(v = overlap.males.real, col = "red")

##########
# Females
# Real value
s1.females.sp <- s1.females[ , c("Longitude", "Latitude", "Species")]
coordinates(s1.females.sp) <- c("Longitude", "Latitude")
proj4string(s1.females.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.females <- kernelUD(s1.females.sp, h = 1.0, grid = rt.sp)
overlap.females <- kerneloverlaphr(kud.females, method = "BA", percent = 95, conditional = T)
o1 <- overlap.females[2, 1]
o2 <- overlap.females[1, 2]
overlap.females.real <- mean(o1, o2)
overlap.females.real

# Plot hr
# vud.females <- getvolumeUD(kud.females)
# vud.females.raster.ngp <- raster(as(vud.females$'Northern giant petrel',"SpatialPixelsDataFrame"))
# vud.females.raster.sgp <- raster(as(vud.females$'Southern giant petrel',"SpatialPixelsDataFrame"))
# vud.females.raster <- stack(vud.females.raster.ngp, vud.females.raster.sgp)
# names(vud.females.raster) <- c("Northern giant petrel", "Southern giant petrel")
# plot(vud.females.raster, col = terrain.colors(100))
# writeRaster(vud.females.raster, file="vud.females.raster.grd", format = "raster")

md <- s1.females[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  m <- comb[ , c(1, i + 3)]
  names(m) <- c("ID", "perm.sp")
  ml <- merge(x = m, y = md, by = "ID", all = F)
  ml$Species <- NULL
  ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sp")]
  coordinates(ml.sp) <- c("Longitude", "Latitude")
  proj4string(ml.sp) <- projt
  kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
  overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = 95, conditional = T)
  ov1 <- overlap.ml[2, 1]
  ov2 <- overlap.ml[1, 2]
  ovr[i] <- mean(ov1, ov2)
}

ovr.females <- ovr

hist(ovr.females, breaks = seq(0, 1, 0.01))
abline(v = overlap.females.real, col = "red")


#---------------------------------------------------
# 2. Swap sex labels

# Set up the frame
nperm = 1000
comb <- matrix(nrow = nrow(ref), ncol = nperm)
comb <- as.data.frame(comb)
comb <- cbind.data.frame(ref, comb)

for (i in 1:1000) {
  comb[ , i + 3] <- rep("Female", nrow(comb))
  one <- sample(x = nrow(ref), size = Males, replace = F)
  comb[one, i + 3] <- "Male"
}

# Run the tests based on the permutation frame
# First, setup

lms <- c(min(tr3$Longitude), max(tr3$Longitude), min(tr3$Latitude), max(tr3$Latitude))
rt <- raster(ext = extent(lms), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), res = 0.1)
rt.sp <- as(rt, "SpatialPixelsDataFrame")

s1 <- tr3
s1 <- s1[s1$sex != "Unknown", ]
s1 <- s1[s1$Distance > 0.4, ]
s1.ngp <- s1[s1$Species == "Northern giant petrel", ]
s1.sgp <- s1[s1$Species == "Southern giant petrel", ]

##########
# NGPs
# Real value
s1.ngp.sp <- s1.ngp[ , c("Longitude", "Latitude", "sex")]
coordinates(s1.ngp.sp) <- c("Longitude", "Latitude")
proj4string(s1.ngp.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.ngp <- kernelUD(s1.ngp.sp, h = 1.0, grid = rt.sp)
overlap.ngp <- kerneloverlaphr(kud.ngp, method = "BA", percent = 95, conditional = T)
o1 <- overlap.ngp[2, 1]
o2 <- overlap.ngp[1, 2]
overlap.ngp.real <- mean(o1, o2)
overlap.ngp.real

# Plot hr
# vud.ngp <- getvolumeUD(kud.ngp)
# vud.ngp.raster.female <- raster(as(vud.ngp$Female,"SpatialPixelsDataFrame"))
# vud.ngp.raster.male<- raster(as(vud.ngp$Male,"SpatialPixelsDataFrame"))
# vud.ngp.raster <- stack(vud.ngp.raster.female, vud.ngp.raster.male)
# names(vud.ngp.raster) <- c("Female", "Male")
# writeRaster(vud.ngp.raster, file="vud.ngp.raster.grd", format = "raster")
# plot(vud.ngp.raster, col = terrain.colors(100))

# Permute
md <- s1.ngp[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  m <- comb[ , c(1, i + 3)]
  names(m) <- c("ID", "perm.sex")
  ml <- merge(x = m, y = md, by = "ID", all = F)
  ml$sex <- NULL
  ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sex")]
  coordinates(ml.sp) <- c("Longitude", "Latitude")
  proj4string(ml.sp) <- projt
  kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
  overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = 95, conditional = T)
  ov1 <- overlap.ml[2, 1]
  ov2 <- overlap.ml[1, 2]
  ovr[i] <- mean(ov1, ov2)
}

ovr.ngp <- ovr

hist(ovr.ngp, breaks = seq(0, 1, 0.01))
abline(v = overlap.ngp.real, col = "red")

##########
# SGP
# Real value
s1.sgp.sp <- s1.sgp[ , c("Longitude", "Latitude", "sex")]
coordinates(s1.sgp.sp) <- c("Longitude", "Latitude")
proj4string(s1.sgp.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.sgp <- kernelUD(s1.sgp.sp, h = 1.0, grid = rt.sp)
overlap.sgp <- kerneloverlaphr(kud.sgp, method = "BA", percent = 95, conditional = T)
o1 <- overlap.sgp[2, 1]
o2 <- overlap.sgp[1, 2]
overlap.sgp.real <- mean(o1, o2)
overlap.sgp.real

# Plot hr
# vud.sgp <- getvolumeUD(kud.sgp)
# vud.sgp.raster.female <- raster(as(vud.sgp$Female,"SpatialPixelsDataFrame"))
# vud.sgp.raster.male<- raster(as(vud.sgp$Male,"SpatialPixelsDataFrame"))
# vud.sgp.raster <- stack(vud.sgp.raster.female, vud.sgp.raster.male)
# names(vud.sgp.raster) <- c("Female", "Male")
# plot(vud.sgp.raster, col = terrain.colors(100))
# writeRaster(vud.sgp.raster, file="vud.sgp.raster.grd", format = "raster")

md <- s1.sgp[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  m <- comb[ , c(1, i + 3)]
  names(m) <- c("ID", "perm.sex")
  ml <- merge(x = m, y = md, by = "ID", all = F)
  ml$sex <- NULL
  ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sex")]
  coordinates(ml.sp) <- c("Longitude", "Latitude")
  proj4string(ml.sp) <- projt
  kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
  overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = 95, conditional = T)
  ov1 <- overlap.ml[2, 1]
  ov2 <- overlap.ml[1, 2]
  ovr[i] <- mean(ov1, ov2)
}

ovr.sgp <- ovr

hist(ovr.sgp, breaks = seq(0, 1, 0.01))
abline(v = overlap.sgp.real, col = "red")

#---------------------------------------------------
#3. Stats

#ovr.males #overlap.males.real
overlap.males.real
mean(ovr.males)
sd(ovr.males)
  #ovr.females #overlap.females.real
overlap.females.real
mean(ovr.females)
sd(ovr.females)
#ovr.ngp #overlap.ngp.real
overlap.ngp.real
mean(ovr.ngp)
sd(ovr.ngp)
#ovr.sgp #overlap.sgp.real
overlap.sgp.real
mean(ovr.sgp)
sd(ovr.sgp)

# One tailed tests
library(Hmisc)
males.p <- sum(overlap.males.real <= ovr.males)
zapsmall(binconf(males.p, nperm, method = 'exact'))

females.p <- sum(overlap.females.real >= ovr.females)
zapsmall(binconf(females.p, nperm, method = 'exact'))

ngp.p <- sum(overlap.ngp.real >= ovr.ngp)
zapsmall(binconf(ngp.p, nperm, method = 'exact'))

sgp.p <- sum(overlap.sgp.real >= ovr.sgp)
zapsmall(binconf(sgp.p, nperm, method = 'exact'))


saveRDS(ovr.males, "./Output/ovr95.males.RDS")
saveRDS(ovr.females, "./Output/ovr95.females.RDS")
saveRDS(ovr.ngp, "./Output/ovr95.ngp.RDS")
saveRDS(ovr.sgp, "./Output/ovr95.sgp.RDS")