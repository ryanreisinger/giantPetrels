# Giant Petrels - Autocorrelated Kernel Density Estimates

## Ryan Reisinger

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

library(ctmm)

# --------------------------------
## Plotting stuff

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
dat <- readRDS("./Output/tracks_trips.RDS")

## Drop missing coordinates
dat <- dat[!is.na(dat$decimal_longitude) & !is.na(dat$decimal_latitude), ]

# Create a 'telemetry' object
dat_telem <- dat[ , c("track_id", "date.time", "decimal_longitude", "decimal_latitude")]
names(dat_telem) <- c("individual.local.identifier", "timestamp", "location.long", "location.lat")
telem <- as.telemetry(object = dat_telem, timeformat = "%Y-%m-%d %H:%M:%S")

# Calculate AKDE for each track

kde_list <- list()

for(i in 1:length(telem)) {
  
  this_dat <- telem[[i]]
  
  # Variograms
  this_vario <- variogram(data = this_dat)
  this_vario_fit <- variogram.fit(this_vario, interactive = FALSE)
  
  # Guess and then fit a ctmm
  ctmm_guess <- ctmm.guess(this_dat,interactive=FALSE) # automated model guess
  ctmm_fit <- ctmm.fit(this_dat, ctmm_guess) # fit the guessed model
  
  # Calculate AKDE and plot
  kde <- akde(this_dat, CTMM = ctmm_fit, weights = FALSE) # calculate akde
  plot(this_dat, UD = kde)
  
  # Add to list
  kde_list <- append(kde_list, list(kde))
}
