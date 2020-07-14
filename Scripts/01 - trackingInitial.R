# Giant Petrels - Initial tracking analyses

## Ryan Reisinger

# --------------------------------
library(ggplot2)
library(maps)
library(argosfilter)
library(geosphere)
library(dplyr)
# library(rgdal)
library(raster)
library(marmap)
library(cowplot)
# library(orsifronts)
# library(rgeos)
# library(maptools)
# library(pals)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

# --------------------------------
## Plotting parameters

## Figure widths in mm
single.col <- 84*0.0393701
double.col <- 140*0.0393701
double.col.sup <- 150*0.0393701

## Scaling for font size
fig.scale <- 8/11

## Custom theme for plots
source("./Scripts/theme_rr.R")

# --------------------------------
## Get data
dat <- read.csv("./Data/GP_tracks_2019-07-01.csv", stringsAsFactors = F)

## Some housekeeping
dat$Culmen_length <- as.numeric(dat$Culmen_length)
dat$Culmen_depth <- as.numeric(dat$Culmen_depth)

# --------------------------------
## Restrictions

## Keep only incubating indiduals
dat <- dat[dat$breeding_stage == "Incubation" | dat$breeding_stage == "incubating", ]

# Remove individual which goes nowhere
dat <- dat[dat$individual_id != "SGP15_KD_SEP_2015", ]

# Updated species codes
dat[dat$sp_code == "NG", "sp_code"] <- "NGP"
dat[dat$sp_code == "SG", "sp_code"] <- "SGP"

# Some locations are missing
dat <- dat[!is.na(dat$decimal_latitude) & !is.na(dat$decimal_longitude), ]

# Some deployment lats and lons
# are text
dat$deployment_decimal_latitude[grepl("S46 56.339", dat$deployment_decimal_latitude)] <- -46.938983
dat$deployment_decimal_longitude[grepl("E37 51.941", dat$deployment_decimal_longitude)] <- 37.865683

dat$deployment_decimal_latitude <- as.numeric(dat$deployment_decimal_latitude)
dat$deployment_decimal_longitude <- as.numeric(dat$deployment_decimal_longitude)

for (i in 1:nrow(dat)) {
  foo <- dat[i, ]
  lon <- foo$deployment_decimal_longitude
  lat <- foo$deployment_decimal_latitude
  if (lon < 0 & lat > 0) {
    lon.new <- lat
    lat.new <- lon
  } else {
    lon.new <- lon
    lat.new <- lat
  }
  if (lat.new > 0) {
    lat.new <- lat.new * -1
  }
  dat[i, "deployment_decimal_longitude"] <- lon.new
  dat[i, "deployment_decimal_latitude"] <- lat.new
}

# --------------------------------
## More accurate deployment locations
dat[dat$individual_id == "SGP07_KD_SEP_2015", "deployment_decimal_latitude"] <- -46.96361
dat[dat$individual_id == "SGP07_KD_SEP_2015", "deployment_decimal_longitude"] <- 37.85207

dat[dat$individual_id == "SGP18_03102016", "deployment_decimal_latitude"] <- -46.96060
dat[dat$individual_id == "SGP18_03102016", "deployment_decimal_longitude"] <- 37.85486

dat[dat$individual_id == "NGP02_092017", "deployment_decimal_latitude"] <- -46.93776
dat[dat$individual_id == "NGP02_092017", "deployment_decimal_longitude"] <- 37.86327

dat[dat$individual_id == "NGP13_26102016", "deployment_decimal_latitude"] <- -46.95383
dat[dat$individual_id == "NGP13_26102016", "deployment_decimal_longitude"] <- 37.86304

dat[dat$individual_id == "NGP05_KD_SEP_2015", "deployment_decimal_latitude"] <- -46.94197 
dat[dat$individual_id == "NGP05_KD_SEP_2015", "deployment_decimal_longitude"] <- 37.87146

dat[dat$individual_id == "NGP17_KD_SEP_2015", "deployment_decimal_latitude"] <- -46.94041
dat[dat$individual_id == "NGP17_KD_SEP_2015", "deployment_decimal_longitude"] <- 37.86775

dat[dat$individual_id == "NGP19_KD_SEP_2015", "deployment_decimal_latitude"] <- -46.95502
dat[dat$individual_id == "NGP19_KD_SEP_2015", "deployment_decimal_longitude"] <- 37.86432

dat[dat$individual_id == "NGP15_26102016", "deployment_decimal_latitude"] <- -46.96903
dat[dat$individual_id == "NGP15_26102016", "deployment_decimal_longitude"] <- 37.85523

dat[dat$individual_id == "NGP01_KD_SEP_2015", "deployment_decimal_latitude"] <- -46.95285
dat[dat$individual_id == "NGP01_KD_SEP_2015", "deployment_decimal_longitude"] <- 37.85796

# --------------------------------
## Create metadata
met <- dat[!duplicated(dat$individual_id), ]

# --------------------------------
## Sex and bill dimensions
met$sex.m <- "Not.sexed"
met[met$track_id == "NGP02_KD_SEP_2015", "sex.m"] <- "Female"
met[met$track_id == "NGP14_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "NGP18_KD_SEP_2015", "sex.m"] <- "Female"
met[met$track_id == "NGP19_KD_SEP_2015", "sex.m"] <- "Male"

met[met$track_id == "SGP01_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "SGP02_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "SGP04_KD_SEP_2015", "sex.m"] <- "Female"
met[met$track_id == "SGP18_KD_SEP_2015", "sex.m"] <- "Female"
met[met$track_id == "SGP20_KD_SEP_2015", "sex.m"] <- "Female"

## Additional sexing
met[met$track_id == "NGP03_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "NGP05_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "NGP011_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "NGP02_26092016", "sex.m"] <- "Female"
met[met$track_id == "NGP05_26102016", "sex.m"] <- "Male"
met[met$track_id == "NGP08_26102016", "sex.m"] <- "Male"
met[met$track_id == "NGP12_26102016", "sex.m"] <- "Female"
met[met$track_id == "NGP13_26102016", "sex.m"] <- "Female"
met[met$track_id == "NGP16_26102016", "sex.m"] <- "Female"
met[met$track_id == "NGP17_26102016", "sex.m"] <- "Female"
met[met$track_id == "NGP21_27102016", "sex.m"] <- "Female"
met[met$track_id == "NGP05_092017", "sex.m"] <- "Male"
met[met$track_id == "NGP19_092017", "sex.m"] <- "Male"

met[met$track_id == "SGP01_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "SGP17_KD_SEP_2015", "sex.m"] <- "Male"
met[met$track_id == "SGP06_03102016", "sex.m"] <- "Female"
met[met$track_id == "SGP07_03102016", "sex.m"] <- "Female"
met[met$track_id == "SGP02_102017", "sex.m"] <- "Male"
met[met$track_id == "SGP09_102017", "sex.m"] <- "Male"
met[met$track_id == "SGP12_102017", "sex.m"] <- "Male"
met[met$track_id == "SGP13_102017", "sex.m"] <- "Male"
met[met$track_id == "SGP23_102017", "sex.m"] <- "Female"
met[met$track_id == "SGP24_102017", "sex.m"] <- "Female"

# Number of individuals sexed molecularly
nrow(met[met$sex.m == "Male" | met$sex.m == "Female", ])

## Plot - Bill dimensions
pdf("./Plots/billDimensions.pdf",
    useDingbats = FALSE,
    width = double.col / fig.scale,
    height = single.col / fig.scale)
ggplot(data = met, aes(x = Culmen_length,
                       y = Culmen_depth,
                       shape = sp_code,
                       colour = sex.m)) +
  scale_shape_manual(values = c(16, 17), name = "Species") +
  scale_colour_manual(values = c("#4daf4a", "#984ea3", "#999999"), name = "Molecular sex") +
  labs(x = "Culmen length (mm)", y = "Culmen depth (mm)") +
  geom_point(size = 1.5) +
  theme_rr() +
  geom_vline(xintercept = 97, colour = "grey")
dev.off()

## Sex
dat$sex <- NA
dat[!is.na(dat$Culmen_length) & dat$Culmen_length < 97, "sex"] <- "Female"
dat[!is.na(dat$Culmen_length) & dat$Culmen_length > 97, "sex"] <- "Male"

## Add molecular sex to the data too
dat <- merge(x = dat, y = met[, c("track_id", "sex.m")],
             all.x = T)

## Check for differences

min(met[!is.na(met$Culmen_length) & met$sex == "Male", "Culmen_length"])
max(met[!is.na(met$Culmen_length) & met$sex == "Male", "Culmen_length"])

min(met[!is.na(met$Culmen_length) & met$sex == "Female", "Culmen_length"])
max(met[!is.na(met$Culmen_length) & met$sex == "Female", "Culmen_length"])

min(met$Culmen_depth, na.rm = T)
max(met$Culmen_depth, na.rm = T)

## Overall
t.test(met[met$sp_code == "NGP",]$Culmen_length,
       met[met$sp_code == "SGP",]$Culmen_length)
t.test(met[met$sp_code == "NGP",]$Culmen_depth,
       met[met$sp_code == "SGP",]$Culmen_depth)

## In males
t.test(met[met$sp_code == "NGP" & met$sex == "Male",]$Culmen_length,
       met[met$sp_code == "SGP" & met$sex == "Male",]$Culmen_length)

## And females
t.test(met[met$sp_code == "NGP" & met$sex == "Female",]$Culmen_length,
       met[met$sp_code == "SGP" & met$sex == "Female",]$Culmen_length)


## Assign sex from molecular results
dat[which(is.na(dat$sex)), "sex"] <-
  dat[which(is.na(dat$sex)), "sex.m"]

## Update metadata
met <- dat[!duplicated(dat$individual_id), ]

## And now drop any individuals where sex is still unknown (should be none)
dat <- dat[!is.na(dat$sex), ]
met <- dat[!duplicated(dat$individual_id), ]

# --------------------------------
## Speed filter
ids <- unique(dat$individual_id)

fn <- dat[0, ]
for (i in 1:length(ids)) {
  print(ids[i])
  f <- dat[dat$individual_id == ids[i], ]
  f$date.time <- paste(f$date, f$time, sep = " ")
  f$date.time <- as.POSIXct(f$date.time, format = "%Y/%m/%d %H:%M:%S", tz = "GMT")
  if (nrow(f) > 5) {
  f$Speed.flag <- vmask(lat = f$decimal_latitude,
                        lon = f$decimal_longitude,
                        dtime = f$date.time,
                        vmax = 30)
  f <- f[f$Speed.flag != "removed", ]
  fn <- rbind.data.frame(fn, f)
  }
}

# What proportion of data are removed?
foo <- function(x, y) {
  df = nrow(x) - nrow(y)
  df = df/nrow(x)
  df = df*100
  return(df)
}
foo(x = dat, y = fn)
rm(foo)

## Replace dat
dat <- fn
rm(fn)

#-----------------------------------
## Distance from deployment

## First update ids
ids <- unique(dat$individual_id)

dat$Distance <- distGeo(p1 = cbind(dat$deployment_decimal_longitude, dat$deployment_decimal_latitude),
                        p2 = cbind(dat$decimal_longitude, dat$decimal_latitude))
dat$Distance <- dat$Distance/1000

# for (i in 1:length(ids)) {
  # d <- dat[dat$individual_id == ids[i], ]
  # p <- ggplot(data = d, aes(x = date.time, y = Distance, colour = sex)) +
  #   geom_point() +
  #   labs(title = ids[i])
  # print(p)
# }

#-----------------------------------
## Maximum distance
l <- dat
l <- l[0, ]
for (i in 1:length(ids)) {
  tm <- dat[dat$individual_id == ids[i], ]
  m <- max(tm$Distance)
  tm$Maxdist <- rep(m, nrow(tm))
  l <- rbind.data.frame(l, tm)
}

## Replace dat
dat <- l
rm(l)

#-----------------------------------
## Update metadata
met <- dat[!duplicated(dat$individual_id), ]

#-----------------------------------
## Automatically assign trips
## based on runs of successive home locations

## Some of this runs code is from:
## https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
## by GUANGMING LANG

## First update ids
ids <- unique(dat$individual_id)

tr2 <- dat
tr2$nest <- "nest"
tr2[tr2$Distance > 0.200, "nest"] <- "trip"

tr2$trip <- NA
tr2[tr2$nest == "nest", "trip"] <- 0
tr2[tr2$nest == "trip", "trip"] <- 1

tr3 <- tr2[0, ]

for (i in 1:length(ids)) {
  print(ids[i])
  
  d <- tr2[tr2$individual_id == ids[i], ]
  
  # Finds runs of locations away from the nest,
  # call these 'realhome'
  d$realhome <- 1 # at nest
  
  runs <- rle(d$trip)
  myruns <- which(runs$values == 1 & runs$lengths >= 3) # Number of locs away
  
  # Find the end position of the runs
  runs.lengths.cumsum <- cumsum(runs$lengths)
  ends <- runs.lengths.cumsum[myruns]
  
  # Find the end position of the runs
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  
  # Store the indices
  dx <- data.frame("starts" = starts, "ends" = ends)
  
  for(j in 1:nrow(dx)) {
    d[dx$starts[j]:dx$ends[j], "realhome"] <- 0 # away from nest
  }
  
  # Plot to check
  # p <- ggplot(data = d,
  #             aes(x = date.time, y = Distance, colour = as.factor(realhome))) +
  #   geom_point() +
  #   labs(title = ids[i])
  # print(p)
  
  ## Drop initial locations away from nest,
  ## GPS setup?
  ## "realhome" == 0
  if (d[1, "realhome"] == 0) {
    dx <- first(which(d$realhome == 1))
    d <- d[dx:nrow(d), ]
  }
  
  ## Then drop trailing locations away from nest,
  ## walking back to base?
  ## Except "SGP11_KD_SEP_2015", which ends at sea
  if (d[nrow(d), "realhome"] == 0 & ids[i] != "SGP11_KD_SEP_2015") {
    dx <- last(which(d$realhome == 1))
    d <- d[1:dx, ]
  }
  
  ## Now trim the leading and trailing at nest locations
  dx.start <- first(which(d$realhome == 0))-1
  dx.end <- last(which(d$realhome == 0))+1
  d <- d[dx.start:dx.end, ]
  
  # Identify trips
  #-------------------
  d$tripno <- NA
  
  runs <- rle(d$realhome)
  myruns <- which(runs$values == 0) # Number of locs away
  
  # Find the end position of the runs
  runs.lengths.cumsum <- cumsum(runs$lengths)
  ends <- runs.lengths.cumsum[myruns]
  
  # Find the end position of the runs
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  
  # Store the indices
  dx <- data.frame("starts" = starts, "ends" = ends)
  
  for(k in 1:nrow(dx)) {
    d[dx$starts[k]:dx$ends[k], "tripno"] <- k # away from nest
  }

  #-------------------
  # Plot to check
  p <- ggplot(data = d,
              aes(x = date.time, y = Distance, colour = as.factor(tripno))) +
    geom_point() +
    labs(title = ids[i])
  print(p)

  ## Add the output
  tr3 <- rbind(tr3, d)
  
}

#-----------------------------------
## Get bathymetry to calculate proportion of locs on land
# Define extent
# minx <- min(tr3$decimal_longitude, na.rm = T) - 1
# maxx <- max(tr3$decimal_longitude, na.rm = T) + 1
# miny <- min(tr3$decimal_latitude, na.rm = T) - 1
# maxy <- max(tr3$decimal_latitude, na.rm = T) + 1
# 
# bath <- getNOAA.bathy(minx,maxx, miny, maxy, resolution = 1)
# bath <- as.raster(bath)
# writeRaster(bath,
#             filename = "./Data/Mapping/bathymetry_large.grd",
#             format = "raster",
#             overwrite = T)
bath <- raster("./Data/Mapping/bathymetry_large.grd")

tr3$depth <- raster::extract(bath, tr3[,c("decimal_longitude", "decimal_latitude")])
rm(bath)

#-----------------------------------
## Displacement plots

## Keep a copy for later summary
tr_all <- tr3

tr3 <- tr3[tr3$realhome == 0, ]
tr3$trip_id <- paste0(tr3$individual_id, "_", tr3$tripno)

trips <- unique(tr3$trip_id)

## Maximum distance per trip
l <- tr3
l <- l[0, ]
for (i in 1:length(trips)) {
  tm <- tr3[tr3$trip_id == trips[i], ]
  m <- max(tm$Distance, na.rm = T)
  tm$trip.Maxdist <- rep(m, nrow(tm))
  l <- rbind.data.frame(l, tm)
}
tr3 <- l
rm(l)

## Look at trip distance distribution
tripdists <- unique(tr3$trip.Maxdist)
hist(tripdists[tripdists < 500], breaks = 100)
abline(v = c(50, 150), col = "red")
plot(sort(tripdists[tripdists < 500]))
abline(h = c(50, 150), col = "red")


## Plots
tr4 <- tr3
tr4$lag <- NA
tr4$Duration <- NA
tr4 <- tr4[0, ]

for (i in 1:length(trips)) {
  trp <- tr3[tr3$trip_id == trips[i], ]
  tms <- diff(trp$date.time)
  units(tms) <- "hours"
  tms <- c(0, tms)
  
  tms <- cumsum(tms)
  trp$lag <- tms
  trp$Duration <- max(tms, na.rm = T)
  tr4 <- rbind.data.frame(tr4, trp)
}

# Look at trip duration and distance together to classify long and short trips
# this produces two groups where distance distinction is very high (~ 800 km)
foo <- tr4[ , c("track_id", "trip.Maxdist", "Duration")]
foo <- unique(foo)
foo <- foo[complete.cases(foo), ]

foo_k <- kmeans(foo[, 2:3], centers = 2)

foo$group <- foo_k$cluster
plot(foo$trip.Maxdist, foo$Duration/24, col = as.factor(foo$group), xlab = "Trip distance (km)", ylab = "Trip duration (days)")

# Try univariate breaks
library(classInt)
classIntervals(var = foo$trip.Maxdist, n = 2, style = "quantile")$brks # Quantile
classIntervals(var = foo$trip.Maxdist, n = 2, style = "kmeans")$brks # univariate k-means
classIntervals(var = foo$trip.Maxdist, n = 2, style = "jenks")$brks # Jenks

# Set threshold to 50 km, according to quantile breaks
thresh <- 50

rm(foo, foo_k)

## Long trips only
tr4long <- tr4[tr4$trip.Maxdist > thresh, ]
disp.long <- ggplot(tr4long, aes(x = lag/24, y = Distance, group = trip_id, colour = sex)) +
  geom_line() +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"),
                      name = "Sex") + 
  theme_rr() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(colour = "black"),
                     plot.title = element_text(size = rel(1))) +
  labs(title = "Distant trips",
       x = "Days since deployment", y = "Distance from nest (km)") +
  scale_x_continuous(limits = c(0, 17))
disp.long


## Short trips only
tr4short <- tr4[tr4$trip.Maxdist < thresh, ]

disp.short <- ggplot(tr4short, aes(x = lag/24, y = Distance, group = trip_id, colour = sex)) +
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
disp.short

pdf("./Plots/displacementPlot.pdf",
    width = double.col/fig.scale, height = (double.col.sup/1.6)/fig.scale,
    useDingbats = FALSE)
cowplot::plot_grid(disp.long, disp.short, nrow = 2)
dev.off()


## Combined plot
tr4sub <- tr4
tr4sub[tr4sub$trip.Maxdist > thresh, "which"] <- "Long trips"
tr4sub[tr4sub$trip.Maxdist < thresh, "which"] <- "Short trips"

dispBoth <- ggplot(tr4sub, aes(x = lag/24, y = Distance, group = trip_id, colour = sex)) +
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

## Numbers in each category
length(unique(tr4sub[tr4sub$sex == "Male" & tr4sub$which == "Short trips", "individual_id"]))
length(unique(tr4sub[tr4sub$sex == "Female" & tr4sub$which == "Short trips", "individual_id"]))

length(unique(tr4sub[tr4sub$sex == "Male" & tr4sub$which == "Long trips", "individual_id"]))
length(unique(tr4sub[tr4sub$sex == "Female" & tr4sub$which == "Long trips", "individual_id"]))

#-----------------------------------
# Summary table

summary <- tr_all[ , c("scientific_name", "individual_id", "Culmen_length", "sex", "Maxdist")]
summary <- summary[!duplicated(summary$individual_id), ]
summary <- summary[!is.na(summary$individual_id), ]

ids <- unique(summary$individual_id)

## Get start and end dates, and calculate mean latititude
## also indicate trip types and proportion of time on land

summary$date.start <- NA
summary$date.end <- NA
summary$duration <- NA
summary$mean.lat <- NA

summary$long.trip <- NA
summary$short.trip <- NA

summary$land.percent

hold <- summary[0, ]
for (i in 1:nrow(summary)) {
  d1 <- summary[summary$individual_id == ids[i], ]
  d2 <- tr_all[tr_all$individual_id == ids[i], ]
  # Start and end date
  d1$date.start <- min(d2$date.time, na.rm = T)
  d1$date.end <- max(d2$date.time, na.rm = T)
  # Duration
  d1$duration <- difftime(d1$date.end, d1$date.start, units = "day")
  # Mean latitude
  d1$mean.lat <- mean(d2$decimal_latitude, na.rm = T)
  
  # Keep only trips, and then calculate
  d2 <- dplyr::filter(d2, !is.na(tripno))
  
  # Percentage of locations on land
  d1$land.percent <- (nrow(d2[d2$depth > -100, ]) / nrow(d2)) * 100
  
  # Calculate Maxdist by trip
  d2 <- d2 %>%
    group_by(tripno) %>%
    summarise(tripmax = max(Distance, na.rm = T))
  
  # Long trips?
  if(any(d2$tripmax > thresh)) {
    d1$long.trip <- "Y"
  } else {
    d1$long.trip <- "N"
  }
  
  # Short trips?
  if(any(d2$tripmax < thresh)) {
    d1$short.trip <- "Y"
  } else {
    d1$short.trip <- "N"
  }
  
  hold <- rbind(hold, d1)
}

summary <- hold

## Write to file
write.csv(summary, file = "./Output/SummaryTable.csv", row.names = F)

# Do males that forage at sea have shorter bills?
t.test(summary[summary$long.trip == "Y" & summary$sex == "Male", "Culmen_length"],
       summary[summary$long.trip == "N" & summary$sex == "Male", "Culmen_length"])


## t-test of trip durations
t.test(summary[summary$long.trip == "Y" & summary$sex == "Male", "duration"],
       summary[summary$long.trip == "Y" & summary$sex == "Female", "duration"])

## Summary of individuals in different groups
foo <- summary
foo <- group_by(summary, scientific_name, long.trip, sex)
foo <- count(foo)
foo
rm(foo)

## Look at relationship between distance, time on land and sex
pdf("./Plots/proportionLand.pdf",
    useDingbats = FALSE,
    width = double.col.sup / fig.scale,
    height = single.col / fig.scale)
ggplot(data = summary, aes(y = land.percent, x = Maxdist, colour = sex, shape = scientific_name)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = thresh, colour = "grey") +
  scale_shape_manual(values = c(16, 17), name = "Species") +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"), name = "Sex") +
  labs(x = "Maximum distance from nest (km)", y = "Proportion of locations on land (%)") +
  theme_rr() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text = element_text(colour = "black"),
                     plot.title = element_text(size = rel(1)))
dev.off()

# Summary of time on land
min(summary[summary$Maxdist < thresh, ]$land.percent, na.rm = T)
max(summary[summary$Maxdist < thresh, ]$land.percent, na.rm = T)
mean(summary[summary$Maxdist < thresh, ]$land.percent, na.rm = T)

#-----------------------------------
## Plot of tracking periods

summary$date.start2 <- format(summary$date.start, format = "%m-%d")
summary$date.end2 <- format(summary$date.end, format = "%m-%d")
summary$year <- format(summary$date.end, format = "%Y")

summary$individual_id <- factor(summary$individual_id, levels = summary$individual_id[order(summary$date.start)])

# Set the year so dates can have the same axis
library(lubridate)
summary$date.start3 <- summary$date.start
year(summary$date.start3) <- 2015
summary$date.end3 <- summary$date.end
year(summary$date.end3) <- 2015

pdf("./Plots/trackingPeriods.pdf",
    width = double.col.sup/fig.scale,
    height = 9/fig.scale,
    useDingbats = FALSE)
p1 <- ggplot(data = summary) +
  geom_linerange(aes(x = individual_id, ymin = date.start3, ymax = date.end3, colour = scientific_name)) +
  scale_colour_manual(values = c("#4daf4a", "#984ea3"), name = "Species") +
  facet_wrap(~year, ncol = 1, scales = "free_y", drop = T) +
  coord_flip() +
  theme_rr() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
print(p1)
dev.off()

# Summarise tracking durations
# This is once trips have been delineated
range(as.numeric(summary$duration))
mean(as.numeric(summary$duration))

# Get the track start and end before trip delineation
foo <- dplyr::select(tr2, individual_id, date.time) %>%
  group_by(individual_id) %>%
  filter(row_number()==1 | row_number()==n())

# Calculate duration
foo$diff <- c(NA, diff(foo$date.time))

# Select every second row
foo <- foo[seq(2, nrow(foo), 2), ]

# Merge with the trip long trip information
foo <- merge(x = foo, y = summary[ , c("individual_id", "long.trip")], by = "individual_id")

# Now check durations
range(as.numeric(foo$diff))
mean(as.numeric(foo$diff))

# Were long.trip tags attached longer?
t.test(
  as.numeric(foo[foo$long.trip == "Y", ]$diff),
as.numeric(foo[foo$long.trip == "N", ]$diff)
)
ggplot(data = foo, aes(x = long.trip, y = diff)) + geom_boxplot() +
  labs(x = "Long trip?", y = "Tracking duration (days)") +
  theme_bw()

# Yes, but this is an artefact of the fieldwork (tags can only be retrieved
# Birds that did not make long trips had enough time to do so
range(as.numeric(foo[foo$long.trip == "N", ]$diff)) # Tracking duration range
mean(as.numeric(foo[foo$long.trip == "N", ]$diff)) # And mean
rm(foo)

# Look at duration of long trip durations
foo <- dplyr::filter(tr3, tr3$trip.Maxdist > 50) %>%
  dplyr::select(., trip_id, date.time)  %>%
  group_by(., trip_id) %>%
  filter(row_number()==1 | row_number()==n())
foo$diff <- c(NA, diff(foo$date.time))/24
foo <- foo[seq(2, nrow(foo), 2), ]

mean(foo$diff)
range(foo$diff)

detach(package:lubridate)
#-----------------------------------
# Save outputs
saveRDS(tr3, "./Output/tracks_trips.RDS")

#-----------------------------------
## Plot bill dimensions again to check for any borderline cases
# pdf("./Plots/billDimensions_maxdist.pdf",
#     width = double.col.sup / fig.scale,
#     height = single.col / fig.scale)
ggplot(data = met, aes(x = Culmen_length,
                       y = Maxdist,
                       shape = sp_code,
                       colour = sex.m)) +
  scale_shape_manual(values = c(16, 17), name = "Species") +
  scale_colour_manual(values = c("#ff7f00", "#a65628", "#999999"), name = "Molecular sex") +
  labs(x = "Culmen length (mm)", y = "Maximum distance (km)") +
  geom_point() +
  theme_rr() +
  geom_vline(xintercept = 97, colour = "grey")
# dev.off()

write.csv(
met[order(met$Culmen_length) , c("individual_id", "Culmen_length", "sex", "sex.m", "Maxdist")],
"./Output/summary_sexing.csv")
