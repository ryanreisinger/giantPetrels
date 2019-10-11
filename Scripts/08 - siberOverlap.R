#### Further stable isotope analyses using SIBER

## Ryan Reisinger, with sections based on the SIBER vignettes by Andrew Jackson

## Install SIBER if neccessary
# library(devtools)
# devtools::install_github("andrewljackson/SIBER")

#----------------------------------------------
library(MVN)
library(vegan)
library(SIBER)
library(ggplot2)
library(Hmisc)
library(viridis)
library(plot3D)

library(dplyr)
library(data.table)

# Set wd
setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

# Set a palette
palette(viridis::viridis(4))

#----------------------------------------------
met <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)
dat <- read.csv("./Data/GP_sia_2019-09-20.csv", stringsAsFactors = F)

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

## Merge with metadata again
d <- merge(x = d,
           y = met[, c("individual_id", "Maxdist", "mean.lat", "long.trip", "short.trip")],
           by.x = "Individual_id",
           by.y = "individual_id",
           all.x = TRUE,
           all.y = FALSE)

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

# The linear model is not great, so proceed without the normalized samples
d <- d[d$Normalized == "N", ]

# -----------------------------------------------
## How many birds sampled
nrow(d)

# Remove birds were we still don't have sex
d <- filter(d, !is.na(sex))
nrow(d)

## How many of these were tracked
nrow(d[!is.na(d$Maxdist), ])

#----------------------------------------------
## Format the dataframe for SIBER
datSiber <- d[ , c("d13C_DP", "d15N_Plas", "sex", "Species")]
names(datSiber) <- c("iso1", "iso2", "group", "community")

# Create a new community field
datSiber$community <- paste0(datSiber$community, "-", datSiber$group)

# Copy the community field back to group
datSiber$group <- datSiber$community

# Create SIBER object
siber.example <- createSiberObject(datSiber)

#----------------------------------------------
## Plot

# Create lists of plotting arguments to be passed onwards to the  
# plotting functions. With p.interval = NULL, these are SEA. NB not SEAc, though
# which is what we will base our overlap calculations on. This implementation 
# needs to be added in a future update. For now, the best way to plot SEAc is to
# add the ellipses manually following the vignette on this topic.
group.ellipses.args  <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2)

par(mfrow = c(1,1))
plotSiberObject(siber.example,
                ax.pad = 0.5, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "o",
                points.order = c(16, 16, 16, 16), #specify symbols
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

#add SEAc
#-------------------------------------------------------------------------------------
## Or plot in GGplot

first.plot <- ggplot(data = datSiber, aes(iso1, iso2)) +
  geom_point(aes(color = factor(group)), size = 2) +
  scale_x_continuous(limits = c(-24, -17.5)) +
  scale_y_continuous(limits = c(11, 16)) +
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a")) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme_bw() +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank())
print(first.plot)

# second.plot <- first.plot + facet_wrap(~factor(group))
# print(second.plot)

# Add frontal locations from Jaeger et al
# Don't add these, since these are all birds,
# not only those foraging at sea

# second.plot <- first.plot + geom_vline(xintercept = -22.9, colour = "gray") + #PF
#   geom_vline(xintercept = -20.1, colour = "gray") #STF
# 
# print(second.plot)

#---------
# Loop over data to create ellipses
grps <- siber.example$all.groups
list.ell <- list()

for (i in 1:length(grps)) {
c.id <- grps[i] # specify the community ID
g.id <- grps[i] # specify the group ID within the community

coords <- addEllipse(siber.example$ML.mu[[c.id]][ , , g.id],
                     siber.example$ML.cov[[c.id]][ , , g.id],
                     n = 500,
                     m = siber.example$sample.sizes[c.id, c.id],
                     p.interval = 0.4,
                     ci.mean = FALSE,
                     small.sample = T)
coords <- as.data.frame(coords)
names(coords) <- c("x", "y")
list.ell[i] <- list(coords)
}

#---------

second.plot <- first.plot +
  geom_path(data = list.ell[[1]], aes(x = x, y = y), colour = "#377eb8", size = 1) +
  geom_path(data = list.ell[[2]], aes(x = x, y = y), colour = "#e41a1c", size = 1) +
  geom_path(data = list.ell[[3]], aes(x = x, y = y), colour = "#984ea3", size = 1) +
  geom_path(data = list.ell[[4]], aes(x = x, y = y), colour = "#4daf4a", size = 1)

print(second.plot)

third.plot <- second.plot +
  geom_polygon(data = list.ell[[1]], aes(x = x, y = y), fill = "#377eb8", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[2]], aes(x = x, y = y), fill = "#e41a1c", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[3]], aes(x = x, y = y), fill = "#984ea3", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[4]], aes(x = x, y = y), fill = "#4daf4a", alpha = 0.3, size = 1)

print(third.plot)

pdf("./Plots/SIAbiplotAll.pdf",
    width = single.col/fig.scale + 1, height = single.col/fig.scale,
    useDingbats = FALSE)
print(third.plot)
dev.off()


#----------------------------------------------
## With respect to other species from Reisinger et al. 2016
# prey <- read.csv("./Data/Prey_SIA.csv", stringsAsFactors = F)
# 
# prey.plot <- third.plot +
#   geom_point(data = prey, aes(x = d13C, y = d15N, shape = Species), inherit.aes = F) +
#   geom_errorbar(data = prey, aes(x = d13C, ymin = d15N - d15N.SD, ymax = d15N + d15N.SD),
#                 inherit.aes = F) + 
#   geom_errorbarh(data = prey, aes(y = d15N, xmin = d13C - d13C.SD, xmax = d13C + d13C.SD),
#                  inherit.aes = F)
# 
# print(prey.plot)

#----------------------------------------------
## SEAc sizes
groupMetricsML(siber.example)


#----------------------------------------------
## Calculate overlap

# NGPs
ellipse1 <- "NGP-Female.NGP-Female"
ellipse2 <- "NGP-Male.NGP-Male"

# Overlap of standard ellipses
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                             p.interval = NULL, n = 100)

# Overlap as a proportion of the non-overlapping area of 
# the two ellipses
prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                          sea.overlap[1] -
                                          sea.overlap[3])
prop.over.NGP <- prop.over
#prop.over.NGP <- sea.overlap[3] #use this instead if not testing  proportional overlap

# SGPs
ellipse1 <- "SGP-Female.SGP-Female"
ellipse2 <- "SGP-Male.SGP-Male"

# Overlap of standard ellipses
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                             p.interval = NULL, n = 100)

# Overlap as a proportion of the non-overlapping area of 
# the two ellipses
prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                 sea.overlap[1] -
                                 sea.overlap[3])
prop.over.SGP <- prop.over
# prop.over.SGP <- sea.overlap[3] #use this instead if not testing  proportional overlap

# Males
ellipse1 <- "NGP-Male.NGP-Male"
ellipse2 <- "SGP-Male.SGP-Male"

# Overlap of standard ellipses
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                             p.interval = NULL, n = 100)

# Overlap as a proportion of the non-overlapping area of 
# the two ellipses
prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                 sea.overlap[1] -
                                 sea.overlap[3])
prop.over.M <- prop.over
# prop.over.M <- sea.overlap[3] #use this instead if not testing  proportional overlap


# Females
ellipse1 <- "NGP-Female.NGP-Female"
ellipse2 <- "SGP-Female.SGP-Female"

# Overlap of standard ellipses
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                             p.interval = NULL, n = 100)

# Overlap as a proportion of the non-overlapping area of 
# the two ellipses
prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                 sea.overlap[1] -
                                 sea.overlap[3])
prop.over.F <- prop.over
# prop.over.F <- sea.overlap[3] #use this instead if not testing  proportional overlap


#----------------------------------------------
## Overlap permutations

#---------------------------------
#1. NGPS
nperm <- 1000

datPerm <- datSiber
datPerm <- datPerm[datPerm$group == "NGP-Male" | datPerm$group == "NGP-Female", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
#Permute
datPermuted <- datPerm
datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
#calculate
siber.perm <- createSiberObject(datPermuted)
ellipse1 <- "NGP-Female.NGP-Female"
ellipse2 <- "NGP-Male.NGP-Male"
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.perm, 
                             p.interval = NULL, n = 100)
prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                 sea.overlap[1] -
                                 sea.overlap[3])
# prop.over <- sea.overlap[3] #use for non-proportional overlap
perm.overlaps[i] <- list(prop.over)
}

overlaps.NGP <- unlist(perm.overlaps)

mean(overlaps.NGP) #mean of permutations
sd(overlaps.NGP) #sd of permuations
prop.over.NGP #real value

# Plot
hist(overlaps.NGP, breaks = seq(0, 1, 0.01))
abline(v = prop.over.NGP, col = "red")

ngp.p <- sum(prop.over.NGP >= overlaps.NGP)
zapsmall(binconf(ngp.p, nperm, method = 'exact'))
#---------------------------------

#2. SGPS
nperm <- 1000

datPerm <- datSiber
datPerm <- datPerm[datPerm$group == "SGP-Male" | datPerm$group == "SGP-Female", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
  #Permute
  datPermuted <- datPerm
  datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
  #calculate
  siber.perm <- createSiberObject(datPermuted)
  ellipse1 <- "SGP-Female.SGP-Female"
  ellipse2 <- "SGP-Male.SGP-Male"
  sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.perm, 
                               p.interval = NULL, n = 100)
  prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                   sea.overlap[1] -
                                   sea.overlap[3])
  # prop.over <- sea.overlap[3] #use for non-proportional overlap
  perm.overlaps[i] <- list(prop.over)
}

overlaps.SGP <- unlist(perm.overlaps)

mean(overlaps.SGP) #mean of permutations
sd(overlaps.SGP) #sd of permuations
prop.over.SGP #real value

# Plot
hist(overlaps.SGP, breaks = seq(0, 1, 0.01))
abline(v = prop.over.SGP, col = "red")

sgp.p <- sum(prop.over.SGP >= overlaps.SGP)
zapsmall(binconf(sgp.p, nperm, method = 'exact'))

#---------------------------------

#3. Females
nperm <- 1000

datPerm <- datSiber
datPerm <- datPerm[datPerm$group == "NGP-Female" | datPerm$group == "SGP-Female", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
  #Permute
  datPermuted <- datPerm
  datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
  #calculate
  siber.perm <- createSiberObject(datPermuted)
  ellipse1 <- "SGP-Female.SGP-Female"
  ellipse2 <- "NGP-Female.NGP-Female"
  sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.perm, 
                               p.interval = NULL, n = 100)
  prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                   sea.overlap[1] -
                                   sea.overlap[3])
  # prop.over <- sea.overlap[3] #use for non-proportional overlap
  perm.overlaps[i] <- list(prop.over)
}

overlaps.F <- unlist(perm.overlaps)

mean(overlaps.F) #mean of permutations
sd(overlaps.F) #sd of permuations
prop.over.F #real value

#Plot
hist(overlaps.F, breaks = seq(0, 1, 0.01))
abline(v = prop.over.F, col = "red")

f.p <- sum(prop.over.F >= overlaps.F)
zapsmall(binconf(f.p, nperm, method = 'exact'))

#---------------------------------

#4. Males
nperm <- 1000

datPerm <- datSiber
datPerm <- datPerm[datPerm$group == "NGP-Male" | datPerm$group == "SGP-Male", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
  #Permute
  datPermuted <- datPerm
  datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
  #calculate
  siber.perm <- createSiberObject(datPermuted)
  ellipse1 <- "SGP-Male.SGP-Male"
  ellipse2 <- "NGP-Male.NGP-Male"
  sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.perm, 
                               p.interval = NULL, n = 100)
  prop.over <- sea.overlap[3] / (sea.overlap[2] + 
                                   sea.overlap[1] -
                                   sea.overlap[3])
  # prop.over <- sea.overlap[3] #use for non-proportional overlap
  perm.overlaps[i] <- list(prop.over)
}

overlaps.M <- unlist(perm.overlaps)

mean(overlaps.M) #mean of permutations
sd(overlaps.M) #sd of permuations
prop.over.M #real value

# Plot
hist(overlaps.M, breaks = seq(0, 1, 0.01))
abline(v = prop.over.M, col = "red")

m.p <- sum(prop.over.M >= overlaps.M)
zapsmall(binconf(m.p, nperm, method = 'exact'))

#---------------------------------
#---------------------------------

# This done in earlier scripts:

# Test residual values
# library(MASS)
# library(ggbeeswarm)
# 
# dat.r <- dat2 # Use all data
# dat.r <- dat.comp # Or use only data with tracking data
# 
# dat.r$longname <- NA
# dat.r[dat.r$ss == "FN", "longname"] <- "NGP female"
# dat.r[dat.r$ss == "MN", "longname"] <- "NGP male"
# dat.r[dat.r$ss == "FS", "longname"] <- "SGP female"
# dat.r[dat.r$ss == "MS", "longname"] <- "SGP male"
# 
# plot(dat.r$d13C, dat.r$d15N)
# 
# # d13C response
# m1 <- lm(formula = d13C ~ d15N, data = dat.r)
# dat.r$d13C.residuals <- m1$residuals
# dat.r$d13C.studentized <- studres(m1)
# 
# ggplot(dat.r, aes(x = longname, y = d13C.studentized, colour = longname)) +
#   # geom_boxplot() +
#   geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
#   geom_quasirandom(width = 0.15) +
#   scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
#   labs(x = "", y = "d13C Studentised residuals")
# 
# 
# # d15N response
# m2 <- lm(formula = d15N ~ d13C, data = dat.r)
# dat.r$d15N.residuals <- m2$residuals
# dat.r$d15N.studentized <- studres(m2)
# 
# ggplot(dat.r, aes(x = d13C, y = d15N.studentized, colour = ss)) +
#   geom_point()
# 
# ggplot(dat.r, aes(x = ss, y = d15N.studentized, colour = ss)) +
#   geom_boxplot()
# 
# 
# # d15N response to latitude
# # m3 <- lm(formula = d15N ~ Lat, data = dat.comp)
# # dat.r$d15N.residuals <- m3$residuals
# # dat.r$d15N.studentized <- studres(m3)
# 
# ggplot(dat.r, aes(x = d13C, y = d15N.studentized, colour = ss)) +
#   geom_point()
# 
# ggplot(dat.r, aes(x = longname, y = d15N.studentized, colour = longname)) +
#   # geom_boxplot() +
#   geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
#   geom_quasirandom(width = 0.15) +
#   scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
#   labs(x = "", y = "d15N Studentised residuals")

#---------------------------------
#---------------------------------

#----------------------------------------------
## Biplot only for LONG TRIPS

## Format the dataframe for SIBER
datSiber <- d[!is.na(d$long.trip) & d$long.trip == "Y", c("d13C_DP", "d15N_Plas", "sex", "Species")]
names(datSiber) <- c("iso1", "iso2", "group", "community")

# Create a new community field
datSiber$community <- paste0(datSiber$community, "-", datSiber$group)

# Copy the community field back to group
datSiber$group <- datSiber$community

# Create SIBER object
siber.example <- createSiberObject(datSiber)

#----------------------------------------------
## Plot

# Create lists of plotting arguments to be passed onwards to the  
# plotting functions. With p.interval = NULL, these are SEA. NB not SEAc, though
# which is what we will base our overlap calculations on. This implementation 
# needs to be added in a future update. For now, the best way to plot SEAc is to
# add the ellipses manually following the vignette on this topic.
group.ellipses.args  <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2)

## plot in GGplot

first.plot <- ggplot(data = datSiber, aes(iso1, iso2)) +
  geom_point(aes(color = factor(group)), size = 2) +
  scale_x_continuous(limits = c(-24, -17.5)) +
  scale_y_continuous(limits = c(11, 16)) +
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a")) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme_bw() +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank())
print(first.plot)

# Add frontal locations from Jaeger et al

second.plot <- first.plot + geom_vline(xintercept = -22.9, colour = "gray") + #PF
  geom_vline(xintercept = -20.1, colour = "gray") #STF

print(second.plot)

#---------
# Loop over data to create ellipses
grps <- siber.example$all.groups
list.ell <- list()

for (i in 1:length(grps)) {
  c.id <- grps[i] # specify the community ID
  g.id <- grps[i] # specify the group ID within the community
  
  coords <- addEllipse(siber.example$ML.mu[[c.id]][ , , g.id],
                       siber.example$ML.cov[[c.id]][ , , g.id],
                       n = 500,
                       m = siber.example$sample.sizes[c.id, c.id],
                       p.interval = 0.4,
                       ci.mean = FALSE,
                       small.sample = T)
  coords <- as.data.frame(coords)
  names(coords) <- c("x", "y")
  list.ell[i] <- list(coords)
}

#---------

second.plot <- second.plot +
  geom_path(data = list.ell[[1]], aes(x = x, y = y), colour = "#377eb8", size = 1) +
  geom_path(data = list.ell[[2]], aes(x = x, y = y), colour = "#e41a1c", size = 1) +
  geom_path(data = list.ell[[3]], aes(x = x, y = y), colour = "#984ea3", size = 1) +
  geom_path(data = list.ell[[4]], aes(x = x, y = y), colour = "#4daf4a", size = 1)

print(second.plot)

third.plot <- second.plot +
  geom_polygon(data = list.ell[[1]], aes(x = x, y = y), fill = "#377eb8", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[2]], aes(x = x, y = y), fill = "#e41a1c", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[3]], aes(x = x, y = y), fill = "#984ea3", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[4]], aes(x = x, y = y), fill = "#4daf4a", alpha = 0.3, size = 1)

print(third.plot)

pdf("./Plots/SIAbiplotLong.pdf",
    width = single.col/fig.scale + 1, height = single.col/fig.scale,
    useDingbats = FALSE)
print(third.plot)
dev.off()

rm(list.ell)