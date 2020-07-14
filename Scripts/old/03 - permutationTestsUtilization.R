#### Permutation tests for distribution overlap

## Ryan Reisinger

library(EMbC)
# library(ggplot2)
# library(argosfilter)
# library(ggmap)
# library(geosphere)
# library(rgdal)
library(raster)
# library(marmap)
# library(rgeos)
# library(maptools)
library(adehabitatHR)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels/")

# Which KUD contour to use?
which.over <- 50
# which.over <- 90

#--------------------------------------------------------
# Get data
tr3 <- readRDS("./Output/tracks_trips.RDS")
tr3 <- tr3[!is.na(tr3$individual_id), ]
met <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)

#--------------------------------------------------------
# Restrict data to foraging locations only
dat <- tr3

## Drop missing coordinates
dat <- dat[!is.na(dat$decimal_longitude) & !is.na(dat$decimal_latitude), ]

## EMbC
ids <- unique(dat$trip_id)
expth <- list()
for (i in ids) {
  this_dat <- dat[dat$trip_id == i, c("date.time", "decimal_longitude", "decimal_latitude")]
  expth <- append(expth, list(this_dat))
  rm(this_dat)
}
mybcp <- stbc(expth, info=-1)

## Get output
dat$embc_velocity <- mybcp@bC@X[,1]
dat$embc_turnrad <- mybcp@bC@X[,2]
dat$embc_label <- mybcp@bC@A

## Fisher-Jenks breaks
v_cut <- classInt::classIntervals(dat$embc_velocity, n = 2, style = "fisher")$brks[2]

## Asign mode
dat$mode <- "transit"
dat[dat$embc_velocity < v_cut, "mode"] <- "restricted"

## Filter
dat <- dat[dat$mode == "restricted", ]

## Return to 'tr3'
tr3 <- dat

#--------------------------------------------------------
# Names for consistency with older code
tr3 <- tr3[ , c("individual_id", "sp_code", "decimal_longitude", "decimal_latitude", "sex", "date.time")]
names(tr3) <- c("ID", "sp_code", "Longitude", "Latitude", "sex", "Date")

tr3$Species <- NA
tr3[tr3$sp_code == "NGP", "Species"] <- "Northern giant petrel"
tr3[tr3$sp_code == "SGP", "Species"] <- "Southern giant petrel"

met <- met[ , c("individual_id", "scientific_name", "sex")]
names(met) <- c("ID", "Species", "sex")

# Add year information
met$year <- NA
met[grep(pattern = "2015", met$ID), "year"] <- 2015
met[grep(pattern = "2016", met$ID), "year"] <- 2016
met[grep(pattern = "2017", met$ID), "year"] <- 2017


# Create ref
ref <- met

# per year
ref_2015 <- ref[ref$year == 2015, ]
ref_2016 <- ref[ref$year == 2016, ]
ref_2017 <- ref[ref$year == 2017, ]

#--------------------------------------------------------
# Count NGPS and Males (will use them as the reference levels)
NGPS_2015 <- length(ref$Species[ref$Species == "Northern Giant Petrel" & ref$year == 2015])
Males_2015 <- length(ref$sex[ref$sex == "Male" & ref$year == 2015])

NGPS_2016 <- length(ref$Species[ref$Species == "Northern Giant Petrel" & ref$year == 2016])
Males_2016 <- length(ref$sex[ref$sex == "Male" & ref$year == 2016])

NGPS_2017 <- length(ref$Species[ref$Species == "Northern Giant Petrel" & ref$year == 2017])
Males_2017 <- length(ref$sex[ref$sex == "Male" & ref$year == 2017])

#---------------------------------------------------
# 1. Swap species labels

# Constrain permutations by year

# Set up the frame
nperm = 1000

#---------
# 2015
comb <- matrix(nrow = nrow(ref), ncol = nperm)
comb <- as.data.frame(comb)
comb <- cbind.data.frame(ref, comb)

for (i in 1:1000) {
comb[ , i + 3] <- rep("Southern giant petrel", nrow(comb))
one <- sample(x = nrow(ref), size = NGPS, replace = F)
comb[one, i + 3] <- "Northern giant petrel"
}

#---------


# Run the tests based on the permutation frame
# First, setup
lms <- c(min(tr3$Longitude) - 5, max(tr3$Longitude) + 5, min(tr3$Latitude) - 5, max(tr3$Latitude) + 5)
rt <- raster(ext = extent(lms), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), res = 0.05)
rt.sp <- as(rt, "SpatialPixelsDataFrame")

s1 <- tr3
s1 <- s1[s1$sex != "Unknown", ]
s1.males <- s1[s1$sex == "Male", ]
s1.females <- s1[s1$sex == "Female", ]

##########
# Males
# Real value
s1.males.sp <- s1.males[ , c("Longitude", "Latitude", "Species")]
coordinates(s1.males.sp) <- c("Longitude", "Latitude")
proj4string(s1.males.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.males <- kernelUD(s1.males.sp, h = 1.0, grid = rt.sp)
overlap.males <- kerneloverlaphr(kud.males, method = "BA", percent = which.over, conditional = T)
o1 <- overlap.males[2, 1]
o2 <- overlap.males[1, 2]
overlap.males.real <- mean(o1, o2)
overlap.males.real

# Permute
md <- s1.males[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  print(i)
m <- comb[ , c(1, i + 3)]
names(m) <- c("ID", "perm.sp")
ml <- merge(x = m, y = md, by = "ID", all = F)
ml$Species <- NULL
ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sp")]
coordinates(ml.sp) <- c("Longitude", "Latitude")
proj4string(ml.sp) <- projt
kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = which.over, conditional = T)
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
overlap.females <- kerneloverlaphr(kud.females, method = "BA", percent = which.over, conditional = T)
o1 <- overlap.females[2, 1]
o2 <- overlap.females[1, 2]
overlap.females.real <- mean(o1, o2)
overlap.females.real


md <- s1.females[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  print(i)
  m <- comb[ , c(1, i + 3)]
  names(m) <- c("ID", "perm.sp")
  ml <- merge(x = m, y = md, by = "ID", all = F)
  ml$Species <- NULL
  ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sp")]
  coordinates(ml.sp) <- c("Longitude", "Latitude")
  proj4string(ml.sp) <- projt
  kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
  overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = which.over, conditional = T)
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
# s1 <- s1[s1$Distance > 0.4, ]
s1.ngp <- s1[s1$Species == "Northern giant petrel", ]
s1.sgp <- s1[s1$Species == "Southern giant petrel", ]

##########
# NGPs
# Real value
s1.ngp.sp <- s1.ngp[ , c("Longitude", "Latitude", "sex")]
coordinates(s1.ngp.sp) <- c("Longitude", "Latitude")
proj4string(s1.ngp.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
kud.ngp <- kernelUD(s1.ngp.sp, h = 1.0, grid = rt.sp)
overlap.ngp <- kerneloverlaphr(kud.ngp, method = "BA", percent = which.over, conditional = T)
o1 <- overlap.ngp[2, 1]
o2 <- overlap.ngp[1, 2]
overlap.ngp.real <- mean(o1, o2)
overlap.ngp.real


# Permute
md <- s1.ngp[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  print(i)
  m <- comb[ , c(1, i + 3)]
  names(m) <- c("ID", "perm.sex")
  ml <- merge(x = m, y = md, by = "ID", all = F)
  ml$sex <- NULL
  ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sex")]
  coordinates(ml.sp) <- c("Longitude", "Latitude")
  proj4string(ml.sp) <- projt
  kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
  overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = which.over, conditional = T)
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
overlap.sgp <- kerneloverlaphr(kud.sgp, method = "BA", percent = which.over, conditional = T)
o1 <- overlap.sgp[2, 1]
o2 <- overlap.sgp[1, 2]
overlap.sgp.real <- mean(o1, o2)
overlap.sgp.real


md <- s1.sgp[ , c("Longitude", "Latitude", "Species", "sex", "ID")]
projt <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
ovr <- rep(NA, nperm)

for (i in 1:nperm) {
  print(i)
  m <- comb[ , c(1, i + 3)]
  names(m) <- c("ID", "perm.sex")
  ml <- merge(x = m, y = md, by = "ID", all = F)
  ml$sex <- NULL
  ml.sp <- ml[ , c("Longitude", "Latitude", "perm.sex")]
  coordinates(ml.sp) <- c("Longitude", "Latitude")
  proj4string(ml.sp) <- projt
  kudm <- kernelUD(ml.sp, h = 1.0, grid = rt.sp)
  overlap.ml <- kerneloverlaphr(kudm, method = "BA", percent = which.over, conditional = T)
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


saveRDS(ovr.males, paste0("./Output/ovr", which.over, ".males.RDS"))
saveRDS(ovr.females, paste0("./Output/ovr", which.over, ".females.RDS"))
saveRDS(ovr.ngp, paste0("./Output/ovr", which.over, ".ngp.RDS"))
saveRDS(ovr.sgp, paste0("./Output/ovr", which.over, ".sgp.RDS"))