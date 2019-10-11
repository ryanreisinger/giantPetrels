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

# Set wd
setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

# Set a palette
palette(viridis::viridis(4))

#----------------------------------------------
## Data

# Load plasma data
dat <- read.csv("./Data/GP_norm_plas.csv",
                stringsAsFactors = F)

dat[dat$species == "N", "species"] <- "NGP"
dat[dat$species == "S", "species"] <- "SGP"

# SGP20 should be female, not male
dat[dat$id == "SGP20_KD_SEP_2015", "sex"] <- "F"

#----------------------------------------------
# Compare plasma values to total blood values
tb <- read.csv("./Data/GP_TB.csv",
               stringsAsFactors = F)
tb <- tb[ , 1:3]

dat.comp <- merge(dat, tb)

plot(dat.comp$DL_d13C, dat.comp$d13C)
abline(a = 0, b = 1)
cor(dat.comp$d13C, dat.comp$DL_d13C)

plot(dat.comp$DL_d15N, dat.comp$d15N)
abline(a = 0, b = 1)
cor(dat.comp$d15N, dat.comp$DL_d15N)

#----------------------------------------------
# Compare plasma values to red blood cell values

rbc <- read.csv("./Data/GP_RBC.csv",
               stringsAsFactors = F)
rbc <- rbc[ , 1:3]
names(rbc) <- c("id", "d13C_RBC", "d15N_RBC")

dat.comp <- merge(dat.comp, rbc)

# d15N
plot(dat.comp$DL_d15N, dat.comp$d15N_RBC)
abline(a = 0, b = 1)

# d13C
plot(dat.comp$DL_d13C, dat.comp$d13C_RBC)
abline(a = 0, b = 1)

#----------------------------------------------
## Summary stats
s.ave <- aggregate(dat[ , c("d13C", "d15N")], by = list(dat$species, dat$sex), FUN = mean)
s.sd <- aggregate(dat[ , c("d13C", "d15N")], by = list(dat$species, dat$sex), FUN = sd)

s.ave$d13Csd <- s.sd$d13C
s.ave$d15Nsd <- s.sd$d15N

names(s.ave) <- c("species", "sex", "ave_d13C", "ave_d15N", "sd_d13C", "sd_d15N")

s.ave <- s.ave[order(s.ave$ave_d13C), ]

#----------------------------------------------
# Merge SIA with summary table for paper

sumt <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)
sumt <- merge(sumt, dat, by.x = "ID", by.y = "id", all = T)


for (i in 1:nrow(sumt)) {
  if (is.na(sumt$sex.x[i])) {
    sumt$sex.x[i] <- sumt$sex.y[i]
  }
}

sumt[sumt$sex.x == "Male", "sex.x"] <- "M"
sumt[sumt$sex.x == "Female", "sex.x"] <- "F"

for (i in 1:nrow(sumt)) {
  if (is.na(sumt$Species[i])) {
    sumt$Species[i] <- sumt$species[i]
  }
}

sumt[sumt$Species == "NGP", "Species"] <- "Northern giant petrel"
sumt[sumt$Species == "SGP", "Species"] <- "Soutern giant petrel"

sumt$tracked <- rep("Y", nrow(sumt))
sumt[is.na(sumt$duration), "tracked"] <- "N"

sumt$sampled <- rep("Y", nrow(sumt))
sumt[is.na(sumt$d13C), "sampled"] <- "N"


write.csv(sumt, "./Output/SummaryTableSIA.csv", row.names = F, na = "-")

#----------------------------------------------
# Correlate values against tracking locations

track <- readRDS("./Output/tr_03-2.RDS")

ids <- unique(track$ID)

hold <- data.frame("id" = ids, "Lat" = numeric(length = 33))

for (i in 1:length(ids)) {
  sub <- track[track$ID == ids[i], ]
# thresh <- 0.25*max(sub$Distance, na.rm = T)
# sub <- sub[sub$Distance > thresh, ]
  hold[i, "Lat"] <- mean(sub$Latitude, na.rm = T)
}

dat.comp <- merge(dat.comp, hold, by = "id")

# Add latitude to summary table
foo <- dat.comp[, c("id", "Lat")]
names(foo) <- c("ID", "Lat")
sumt <- merge(sumt, foo)

# Correlation between Latitude and d13C
cor(dat.comp$d13C, dat.comp$Lat)
ggplot(data = dat.comp, aes(x = d13C, y = Lat, colour = ss)) +
  geom_point()

# But not as strong for d15N, although the shape of the relationship is similar
cor(dat.comp$d15N, dat.comp$Lat)
ggplot(data = dat.comp, aes(x = d15N, y = Lat, colour = ss)) +
  geom_point()

# Relationship is much weaker in plasma
cor(dat.comp$DL_d13C, dat.comp$Lat)
ggplot(data = dat.comp, aes(x = DL_d13C, y = Lat, colour = ss)) +
  geom_point()

cor(dat.comp$DL_d15N, dat.comp$Lat)
ggplot(data = dat.comp, aes(x = DL_d15N, y = Lat, colour = ss)) +
  geom_point()

#----------------------------------------------
# Correlate d13C v d15N
ggplot(dat.comp, aes(x = d13C, y = d15N)) +
  geom_point() +
  geom_smooth(method = "lm")
cor(dat.comp$d13C, dat.comp$d15N)

summary(lm(formula = d15N ~ d13C + Lat, data = dat.comp))
summary(lm(formula = d15N ~ d13C, data = dat.comp))

#----------------------------------------------
# 3D plot

library(plot3D)

dat.comp$group <- NA
dat.comp[dat.comp$ss == "MN", "group"] <- 1
dat.comp[dat.comp$ss == "FN", "group"] <- 2
dat.comp[dat.comp$ss == "MS", "group"] <- 3
dat.comp[dat.comp$ss == "FS", "group"] <- 4

# Plotting symbols
dat$symbol <- NA
dat.comp[dat.comp$group == 1, "symbol"] <- 15
dat.comp[dat.comp$group == 2, "symbol"] <- 16
dat.comp[dat.comp$group == 3, "symbol"] <- 17
dat.comp[dat.comp$group == 4, "symbol"] <- 18

scatter3D(x = dat.comp$d13C, y = dat.comp$Lat, z = dat.comp$d15N, colvar = dat.comp$group,
          col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"),
          bty = "b",
          ticktype = "detailed",
          pch = dat.comp$symbol,
          phi = 17, theta = 35,
          xlab = "d13C", ylab = "Latitude", zlab = "d15N",
          colkey = list(at = c(1, 2, 3, 4),labels = c("NGP Male", "NGP Female", "SGP Male", "SGP Female")))

dat$symbol <- NULL

#----------------------------------------------
## Test for multivariate normality
mvn(data = dat[dat$sex == "M" & dat$species == "NGP", c("d13C", "d15N")], mvnTest = "mardia") #fine
mvn(data = dat[dat$sex == "F" & dat$species == "NGP", c("d13C", "d15N")], mvnTest = "mardia") #fine
mvn(data = dat[dat$sex == "M" & dat$species == "SGP", c("d13C", "d15N")], mvnTest = "mardia") #fine
mvn(data = dat[dat$sex == "F" & dat$species == "SGP", c("d13C", "d15N")], mvnTest = "mardia") #fine

#----------------------------------------------
## PERMANOVA/MANOVA

# Manova
dat2 <- dat
dat2$group <- paste0(dat2$species, "-", dat2$sex)
man <- manova(data = dat2, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")
summary.aov(man)

# Compare males
dat3 <- dat2[dat2$sex == "M", ]
man <- manova(data = dat3, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")

# Compare females
dat3 <- dat2[dat2$sex == "F", ]
man <- manova(data = dat3, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")

# Compare NGPs
dat3 <- dat2[dat2$species == "NGP", ]
man <- manova(data = dat3, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")

# Compare SGPs
dat3 <- dat2[dat2$species == "SGP", ]
man <- manova(data = dat3, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")

# Compare males which remained near the island
shortmales <- unique(track[track$Maxdist < 100, "ID"])
dat3 <- dat2[dat2$id %in% shortmales, ]
man <- manova(data = dat3, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")

# Summary stats for them
mean(dat3[dat3$species == "SGP", "d13C"])
sd(dat3[dat3$species == "SGP", "d13C"])
mean(dat3[dat3$species == "SGP", "d15N"])
sd(dat3[dat3$species == "SGP", "d15N"])

mean(dat3[dat3$species == "NGP", "d13C"])
sd(dat3[dat3$species == "NGP", "d13C"])
mean(dat3[dat3$species == "NGP", "d15N"])
sd(dat3[dat3$species == "NGP", "d15N"])

# What about NGP males that stayed near the island v. those that made distant trips?
shortmales <- shortmales[grep("NGP", shortmales)]
d3 <- dat2[dat2$species == "NGP" & dat2$sex == "M", ]
d3[d3$id %in% shortmales, "group"] <- "NGP-M-short"

# Difference?
man <- manova(data = d3, cbind(d13C, d15N) ~ group)
summary(man, test = "Wilks")

# Summary stats
# Stayed near isdland
mean(d3[d3$group == "NGP-M-short", "d13C"])
sd(dat3[d3$group == "NGP-M-short", "d13C"])
mean(d3[d3$group == "NGP-M-short", "d15N"])
sd(d3[d3$group == "NGP-M-short", "d15N"])

# For the other males
mean(d3[d3$group == "NGP-M", "d13C"])
sd(d3[d3$group == "NGP-M", "d13C"])
mean(d3[d3$group == "NGP-M", "d15N"])
sd(d3[d3$group == "NGP-M", "d15N"])

# Permanova
adonis(dat2[ , c("d13C", "d15N")] ~ group, data = dat2, method = "euclidean")

#----------------------------------------------
## Format the dataframe for SIBER
datSiber <- dat[ , 2:(ncol(dat) - 1)]
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

#Add frontal locations from Jaeger et al
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

third.plot <- second.plot +
  geom_path(data = list.ell[[1]], aes(x = x, y = y), colour = "#377eb8", size = 1) +
  geom_path(data = list.ell[[2]], aes(x = x, y = y), colour = "#e41a1c", size = 1) +
  geom_path(data = list.ell[[3]], aes(x = x, y = y), colour = "#984ea3", size = 1) +
  geom_path(data = list.ell[[4]], aes(x = x, y = y), colour = "#4daf4a", size = 1)

print(third.plot)

third.plot <- second.plot +
  geom_polygon(data = list.ell[[1]], aes(x = x, y = y), fill = "#377eb8", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[2]], aes(x = x, y = y), fill = "#e41a1c", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[3]], aes(x = x, y = y), fill = "#984ea3", alpha = 0.3, size = 1) +
  geom_polygon(data = list.ell[[4]], aes(x = x, y = y), fill = "#4daf4a", alpha = 0.3, size = 1)

print(third.plot)

#----------------------------------------------
## With respect to other species from Reisinger et al. 2016
prey <- read.csv("./Data/Prey_SIA.csv", stringsAsFactors = F)

prey.plot <- third.plot +
  geom_point(data = prey, aes(x = d13C, y = d15N, shape = Species), inherit.aes = F) +
  geom_errorbar(data = prey, aes(x = d13C, ymin = d15N - d15N.SD, ymax = d15N + d15N.SD),
                inherit.aes = F) + 
  geom_errorbarh(data = prey, aes(y = d15N, xmin = d13C - d13C.SD, xmax = d13C + d13C.SD),
                 inherit.aes = F)

print(prey.plot)

#----------------------------------------------
## SEAc sizes
groupMetricsML(siber.example)


#----------------------------------------------
## Calculate overlap

# NGPs
ellipse1 <- "NGP-F.NGP-F"
ellipse2 <- "NGP-M.NGP-M"

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
ellipse1 <- "SGP-F.SGP-F"
ellipse2 <- "SGP-M.SGP-M"

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
ellipse1 <- "NGP-M.NGP-M"
ellipse2 <- "SGP-M.SGP-M"

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
ellipse1 <- "NGP-F.NGP-F"
ellipse2 <- "SGP-F.SGP-F"

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
datPerm <- datPerm[datPerm$group == "NGP-M" | datPerm$group == "NGP-F", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
#Permute
datPermuted <- datPerm
datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
#calculate
siber.perm <- createSiberObject(datPermuted)
ellipse1 <- "NGP-F.NGP-F"
ellipse2 <- "NGP-M.NGP-M"
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
datPerm <- datPerm[datPerm$group == "SGP-M" | datPerm$group == "SGP-F", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
  #Permute
  datPermuted <- datPerm
  datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
  #calculate
  siber.perm <- createSiberObject(datPermuted)
  ellipse1 <- "SGP-F.SGP-F"
  ellipse2 <- "SGP-M.SGP-M"
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
datPerm <- datPerm[datPerm$group == "NGP-F" | datPerm$group == "SGP-F", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
  #Permute
  datPermuted <- datPerm
  datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
  #calculate
  siber.perm <- createSiberObject(datPermuted)
  ellipse1 <- "SGP-F.SGP-F"
  ellipse2 <- "NGP-F.NGP-F"
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
datPerm <- datPerm[datPerm$group == "NGP-M" | datPerm$group == "SGP-M", ] #select the group of interest

perm.overlaps <- list()

for (i in 1:nperm) {
  #Permute
  datPermuted <- datPerm
  datPermuted[ , c("group", "community")] <- datPerm[sample(nrow(datPerm)), c("group", "community")]
  #calculate
  siber.perm <- createSiberObject(datPermuted)
  ellipse1 <- "SGP-M.SGP-M"
  ellipse2 <- "NGP-M.NGP-M"
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
# Test residual values
library(MASS)
library(ggbeeswarm)

dat.r <- dat2 # Use all data
dat.r <- dat.comp # Or use only data with tracking data

dat.r$longname <- NA
dat.r[dat.r$ss == "FN", "longname"] <- "NGP female"
dat.r[dat.r$ss == "MN", "longname"] <- "NGP male"
dat.r[dat.r$ss == "FS", "longname"] <- "SGP female"
dat.r[dat.r$ss == "MS", "longname"] <- "SGP male"

plot(dat.r$d13C, dat.r$d15N)

# d13C response
m1 <- lm(formula = d13C ~ d15N, data = dat.r)
dat.r$d13C.residuals <- m1$residuals
dat.r$d13C.studentized <- studres(m1)

ggplot(dat.r, aes(x = longname, y = d13C.studentized, colour = longname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "", y = "d13C Studentised residuals")


# d15N response
m2 <- lm(formula = d15N ~ d13C, data = dat.r)
dat.r$d15N.residuals <- m2$residuals
dat.r$d15N.studentized <- studres(m2)

ggplot(dat.r, aes(x = d13C, y = d15N.studentized, colour = ss)) +
  geom_point()

ggplot(dat.r, aes(x = ss, y = d15N.studentized, colour = ss)) +
  geom_boxplot()


# d15N response to latitude
# m3 <- lm(formula = d15N ~ Lat, data = dat.comp)
# dat.r$d15N.residuals <- m3$residuals
# dat.r$d15N.studentized <- studres(m3)

ggplot(dat.r, aes(x = d13C, y = d15N.studentized, colour = ss)) +
  geom_point()

ggplot(dat.r, aes(x = longname, y = d15N.studentized, colour = longname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "", y = "d15N Studentised residuals")
