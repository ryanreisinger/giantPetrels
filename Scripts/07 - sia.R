# Stable isotopes - plots and SIBER

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

library(ggplot2)
library(dplyr)
library(data.table)
library(SIBER)

met <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)
dat <- read.csv("./Data/GP_sia_2019-09-20.csv", stringsAsFactors = F)

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
## Initial look
ggplot(data = dat, aes(x = d13C, y = d15N, colour = Species)) +
  geom_point() +
  facet_wrap(~Sample_Type, ncol = 1)

# ----------------------------------
## Match tracking info
d <- merge(x = dat, y = met, by.x = "Individual_id", by.y = "individual_id",
           all.x = T, all.y = F)

## Merge Culmen_length and Culmen_Length
d$Culmen_Length <- apply(cbind(d$Culmen_Length, d$Culmen_length), 1, mean, na.rm = T)
d$Culmen_length <- NULL

## Where sex is missing, assing using culmen length
d[is.na(d$sex) & !is.na(d$Culmen_Length) & d$Culmen_Length > 97, "sex"] <- "Male"
d[is.na(d$sex) & !is.na(d$Culmen_Length) & d$Culmen_Length < 97, "sex"] <- "Female"

## Wide table for easier plotting
d <- dcast(data = setDT(d), Individual_id + Species + sex ~ Sample_Type, value.var = c("d13C", "d15N", "CN.ratio"))

# -----------------------------------------------
## Lipid normalization of d13C where there are no delipidated plasma samples

# Which are these? (13 values)
d[is.na(d$d13C_DP), ]

# Approach described in the supplement
d13C_delta <- d$d13C_Plas - d$d13C_DP # Equation 1
CN_delta <- d$CN.ratio_DP - d$CN.ratio_Plas # Equation 2
m1 <- lm(d13C_delta ~ CN_delta) # Linear regression

# Summary of lm
summary(m1)

# Predict

# First, we need to create a dummy df with average CN
foo <- d
foo$CN_delta <- mean(d$CN.ratio_DP, na.rm = T) - d$CN.ratio_Plas

# Then we can predict
preds <- predict.lm(m1, newdata = foo)
d$Predicted_d13C_DP <- d$d13C_Plas + preds

# Add a flag for tables
d$Normalized <- "N"
d[is.na(d$d13C_DP), "Normalized"] <- "Y"

# Replace missing values with predictions
d[is.na(d$d13C_DP), "d13C_DP"] <- d[is.na(d$d13C_DP), "Predicted_d13C_DP"]

# Tidy up
rm(d13C_delta, CN_delta, m1, preds, foo)

## Merge with metadata again
d <- merge(x = d,
           y = met[, c("individual_id", "Maxdist", "mean.lat", "long.trip", "short.trip")],
           by.x = "Individual_id",
           by.y = "individual_id",
           all.x = TRUE,
           all.y = FALSE)

# The linear model is not great, so proceed without the normalized samples
d <- d[d$Normalized == "N", ]

## How many birds sampled
nrow(d)

# Remove birds were we still don't have sex
d <- filter(d, !is.na(sex))
nrow(d)

## How many of these were tracked
nrow(d[!is.na(d$Maxdist), ])

#----------------------------------------------
## Basic plotting and metrics in SIBER

# Format the dataframe for SIBER
datSiber <- d[ , c("d13C_DP", "d15N_Plas", "Species", "sex")]
names(datSiber) <- c("iso1", "iso2", "group", "community")

# Create a new community field
datSiber$community <- paste0(datSiber$group, "-", datSiber$community)
#Copy the community field back to group
datSiber$group <- datSiber$community

# Create siber object
siber.example <- createSiberObject(datSiber)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 500, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# Palette
palette(viridis::viridis(4))

# Plot
par(mfrow = c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

# Standard Ellipses, 40% of data
par(mfrow = c(1,1))

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2, small.sample = F) # Does not plot if small sample = T
group.hull.args      <- list(lty = 2, col = "grey20")

# This time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber.example,
                ax.pad = 0.5, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "o",
                iso.order = c(1,2),
                points.order = c(16, 16, 16, 16), #specify symbols
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
                cex = 0.5
)

# Or add the elipses separately, if turned off in the above plot
# Add ellipses contained X% of data, specified by p.interval
plotGroupEllipses(siber.example, n = 500, p.interval = 0.4,
                  lty = 1, lwd = 2, small.sample = T) # Only plots one ellipse... sample size problem?

# Check the group names
siber.example$group.names

# Calculate summary statistics
group.ML <- groupMetricsML(siber.example)
print(group.ML)

#----------------------------------------------
## Bayesian ellipses

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)

# Calculate ellipses
# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col = "red", pch = "x", lwd = 2)

# NGPs have much larger ellipses than SGPs

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp <- hdrcde::hdr(x)$hdr},
  prob = cr.p)

# Do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp <- hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes = T)

#----------------------------------------------
## Overlap of Bayesian ellipses

# Returns first 10 posterior draws
# compare sexes
bayesianOverlap("NGP-Female.NGP-Female", "NGP-Male.NGP-Male", ellipses.posterior = ellipses.posterior, do.plot = F)
bayesianOverlap("SGP-Female.SGP-Female", "SGP-Male.SGP-Male", ellipses.posterior = ellipses.posterior, do.plot = F)