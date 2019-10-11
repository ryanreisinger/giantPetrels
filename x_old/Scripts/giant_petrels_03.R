# Giant Petrels

## Ryan Reisinger

# --------------------------------
library(ggplot2)
library(maps)
library(argosfilter)
library(geosphere)
library(dplyr)
# library(rgdal)
# library(raster)
# library(marmap)
# library(cowplot)
# library(orsifronts)
# library(rgeos)
# library(maptools)
# library(pals)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

# --------------------------------
## Plotting stuff

## Figure widths in mm
single.col <- 84*0.0393701
double.col <- 174*0.0393701
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

## Remove individuals without bill measurements
dat <- dat[!is.na(dat$Culmen_length), ]
dat <- dat[!is.na(dat$Culmen_depth), ]

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
## Better deployment locations
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

## Plot - Bill dimensions
pdf("./Plots/02/billDimensions.pdf",
    width = double.col.sup / fig.scale,
    height = single.col / fig.scale)
ggplot(data = met, aes(x = Culmen_length,
                       y = Culmen_depth,
                       shape = sp_code,
                       colour = sex.m)) +
  scale_shape_manual(values = c(16, 17), name = "Species") +
  scale_colour_manual(values = c("#ff7f00", "#a65628", "#999999"), name = "Molecular sex") +
  labs(x = "Culmen length (mm)", y = "Culmen depth (mm)") +
  geom_point() +
  theme_rr() +
  geom_vline(xintercept = 97, colour = "grey")
dev.off()

## Sex
dat$sex <- NA
dat[dat$Culmen_length < 97, "sex"] <- "Female"
dat[dat$Culmen_length > 97, "sex"] <- "Male"

## Update metadata
met <- dat[!duplicated(dat$individual_id), ]

## Check for differences

## Overall
t.test(met[met$sp_code == "NGP",]$Culmen_length,
       met[met$sp_code == "SGP",]$Culmen_length)

## In males
t.test(met[met$sp_code == "NGP" & met$sex == "Male",]$Culmen_length,
       met[met$sp_code == "SGP" & met$sex == "Male",]$Culmen_length)

## And females
t.test(met[met$sp_code == "NGP" & met$sex == "Female",]$Culmen_length,
       met[met$sp_code == "SGP" & met$sex == "Female",]$Culmen_length)

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
#   d <- dat[dat$individual_id == ids[i], ]
#   p <- ggplot(data = d, aes(x = date.time, y = Distance, colour = sex)) +
#     geom_point() +
#     labs(title = ids[i]) 
#   print(p)
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

## Some runs code from:
## https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
## By GUANGMING LANG

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
  # p <- ggplot(data = d,
  #             aes(x = date.time, y = Distance, colour = as.factor(tripno))) +
  #   geom_point() +
  #   labs(title = ids[i])
  # print(p)

  ## Add the output
  tr3 <- rbind(tr3, d)
  
}

#-----------------------------------
#-----------------------------------
## Displacement plots

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
tr4 <- tr4[0, ]

for (i in 1:length(trips)) {
  trp <- tr3[tr3$trip_id == trips[i], ]
  tms <- diff(trp$date.time)
  units(tms) <- "hours"
  tms <- c(0, tms)
  
  tms <- cumsum(tms)
  trp$lag <- tms
  tr4 <- rbind.data.frame(tr4, trp)
}

## Long trips only
tr4long <- tr4[tr4$trip.Maxdist > 50, ]
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
tr4short <- tr4[tr4$trip.Maxdist < 50, ]

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
