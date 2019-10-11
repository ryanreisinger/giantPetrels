#### Stable isotope analyeses using SIBER

## Ryan Reisinger, modified from the SIBER vignettes by Andrew Jackson

# library(devtools)
# devtools::install_github("andrewljackson/SIBER",
#                          build_vignettes = TRUE, force = TRUE)

library(SIBER)

# Set wd
setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

#----------------------------------------------
## Data

# Load stable isotope data
dat <- read.csv("./Data/Giant petrels_isotopes & deployments_2015 data.csv",
                stringsAsFactors = F)

# Assign sex where missing
dat[dat$Ind == "SGP20", "Sex"] <- "F"

# Re-organize dataframe
dat <- dat[ , c("Ind", "dC", "dN", "Sp", "Sex")]

#----------------------------------------------
## Summary stats
s.ave <- aggregate(dat[ , c("dC", "dN")], by = list(dat$Sp, dat$Sex), FUN = mean)
s.sd <- aggregate(dat[ , c("dC", "dN")], by = list(dat$Sp, dat$Sex), FUN = sd)


#----------------------------------------------
## Basic plotting and metrics in SIBER

# Format the dataframe for SIBER
datSiber <- dat[ , 2:ncol(dat)]
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

# Retunrs first 10 posterior draws
# compare sexes
bayesianOverlap("NGP-F.NGP-F", "NGP-M.NGP-M", ellipses.posterior = ellipses.posterior, do.plot = F)
bayesianOverlap("SGP-F.SGP-F", "SGP-M.SGP-M", ellipses.posterior = ellipses.posterior, do.plot = F)
