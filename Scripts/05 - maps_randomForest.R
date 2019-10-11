setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

library(raster)
library(ggplot2)
library(pals)
library(marmap)
library(ncdf4)
library(rgeos)
library(gridExtra)
library(grid)
library(reshape2)
library(plyr)
library(rgdal)
library(maptools)

#-----------------------------------
# Environmental data have been extracted
# using raadtools run on R server
# (script 04 -)

# Read in these data
tr <- readRDS("./Output/tracks_trips_envar.RDS")

#-----------------------------------------------------
# Maps
#-----------------------------------------------------

# Define extent
minx <- min(tr$lon)
maxx <- max(tr$lon)
miny <- min(tr$lat)
maxy <- max(tr$lat)

# For presentation
minx <- 5
maxx <- 65
miny <- -65
maxy <- -30

## Plotting stuff

## Figure widths in mm
single.col <- 84*0.0393701
double.col <- 174*0.0393701
double.col.sup <- 150*0.0393701

## Scaling for font size
fig.scale <- 8/11

## Custom theme for plots
source("./Scripts/theme_rr.R")


# Get bathymetry
# b <- getNOAA.bathy(minx,maxx, miny, maxy, resolution = 1)
# b2 <- as.raster(b)
# writeRaster(b2, filename = "b2.grd", format = "raster", overwrite = T)
# b2 <- raster("b2.grd")

# A smaller version for large plot
# b3 <- getNOAA.bathy(minx,maxx, miny, maxy, resolution = 4)
# b4 <- as.raster(b3)
# b4[b4>0] <- NA #set land
# writeRaster(b4, filename = "b4.grd", format = "raster")

# Or read in existing bathymetry
b4 <- raster("./Data/Mapping/bathymetry.grd")

# Crop
b4 <- crop(b4, extent(minx, maxx, miny, maxy))

# Convert to points
b4.p <- rasterToPoints(b4)
b4.p <- data.frame(b4.p)
colnames(b4.p) <- c("lon", "lat", "val")

# Read in saved SST
sst <- raster("./Data/Mapping/ghrsst.grd")
sst.large <- aggregate(sst, 10)

# Convert to points
sst.p <- rasterToPoints(sst.large)
sst.p <- data.frame(sst.p)
colnames(sst.p) <- c("lon", "lat", "val")

# Get Southern Ocean fronts from Park & Durand 2019
# https://doi.org/10.17882/59800

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

# Crop the fronts
frnts <- frnts[frnts$lat >= miny & frnts$lat <= maxy, ]
frnts <- frnts[frnts$lon >= minx & frnts$lon <= maxx, ]

# Create a dummy variable for plotting
tr$Species2 <- paste(tr$sp_code, tr$sex, sep = " - ")

# World map
world <- borders("world", colour = "gray50", fill = "gray50", xlim = c(minx, maxx), ylim = c(miny, maxy))


# # Bathymetry
# # Coloured
# p1 <- ggplot(data = b4.p, aes(x = lon, y = lat))
# p1 <- p1 + geom_tile(aes(fill = val))
# p1 <- p1 + scale_x_continuous(expand = c(0,0))
# p1 <- p1 + scale_y_continuous(expand = c(0,0))
# p1 <- p1 + scale_fill_gradient2(low = "dodgerblue4", mid = "gainsboro", high = "burlywood4", midpoint = 0)
# p1 <- p1 + geom_point(data = dat, aes(x = dat$Longitude, y = dat$Latitude, colour = Species)) + facet_wrap(Species ~ sex) +
#   scale_colour_manual(values = c("#4daf4a", "#984ea3"))
# p1 <- p1 + coord_quickmap()
# p1 <- p1 + theme_bw() + theme(panel.grid.minor = element_blank(),
#                               panel.grid.major = element_blank())
# p1

# # Greyscale
# p2 <- ggplot(data = b4.p, aes(x = lon, y = lat))
# p2 <- p2 + geom_tile(aes(fill = val))
# # p2 <- p2 + scale_fill_gradient(low = "gray60", high = "gray99", name = "Depth (m)", limits = c(-7000, 0)) #if land is NA
# # p2 <- p2 + geom_path(data = myOrsi, aes(x = long, y = lat, group = group), linetype = 2)
# p2 <- p2 + geom_point(data = dat, aes(x = dat$Longitude, y = dat$Latitude, colour = Species), size = 0.6) + facet_wrap(Species ~ sex) +
#   scale_colour_manual(values = c("#4daf4a", "#984ea3"), guide = FALSE)
# p2 <- p2 + coord_quickmap()
# p2 <- p2 + scale_x_continuous(expand = c(0,0), limits = c(minx, maxx))
# p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(miny, maxy))
# p2 <- p2 + theme_bw() + theme(panel.grid.minor = element_blank(),
#                               panel.grid.major = element_blank(),
#                               strip.background = element_blank(),
#                               panel.border = element_rect(colour = "black"),
#                               axis.title.x = element_blank(),
#                               axis.title.y = element_blank())
# world <- borders("world", colour = "black", fill = "grey", xlim = c(min(dat$Longitude), max(dat$Longitude)), ylim = c(min(dat$Latitude), max(dat$Latitude)))
# p2 <- p2 + world #+ xlim(minx, maxx) + ylim(miny, maxy)
# p2

# With SST
p2 <- ggplot(data = sst.p, aes(x = lon, y = lat))
p2 <- p2 + geom_tile(aes(fill = val))
p2 <- p2 + scale_fill_gradientn(colours=parula(100),
                                guide = FALSE,
                                name = "SST",
                                limits = c(-2, 24))
# p2 <- p2 + geom_path(data = myOrsi, aes(x = long, y = lat, group = group), colour = "white")
p2 <- p2 + geom_point(data = tr, aes(x = tr$lon, y = tr$lat, colour = sp_code), size = 0.6) +
  facet_wrap(~ Species2, ncol = 2) +
  # scale_colour_manual(values = c("#4daf4a", "#984ea3"), guide = FALSE) # coloured points
  scale_colour_manual(values = c("black", "black"), guide = FALSE)
p2 <- p2 + coord_quickmap()
p2 <- p2 + labs(title = "a")
p2 <- p2 + scale_x_continuous(expand = c(0,0), limits = c(minx, maxx))
p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(miny, maxy))
p2 <- p2 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              axis.text = element_text(colour = "black"),
                              panel.border = element_rect(colour = "black"),
                              plot.title = element_text(size = rel(1)),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank())
world <- borders("world", colour = "black", fill = "grey", xlim = c(minx, maxx), ylim = c(miny, maxy))
p2 <- p2 + world #+ xlim(minx, maxx) + ylim(miny, maxy)
p2 <- p2 + annotate("point", x = 37.74, y = -46.91, size = 3, colour = "red")
# p2

#-----------------------------------

# Close-up map of island

# Get island shapefile
island = readOGR(dsn = "./Data/Mapping", layer = "Islands_Polygonizer")
island@data$id = rownames(island@data)
island.points = fortify(island, region = "id")
island.df = join(island.points, island@data, by = "id")

# Crop to Marion only
island.df <- island.df[island.df$lat < -46.7, ]

# Quick check
ggplot(island.df) + aes(long, lat, group = group) + geom_polygon(fill = "grey") + geom_path(color = "black")

# Crop raster
# Bathymetry
b5 <- getNOAA.bathy(37.45, 38.2, -47.2, -46.8, resolution = 1)
b5 <- marmap::as.raster(b5)
b5[b5 > 0] <- NA # Set land
b5.p <- rasterToPoints(b5)
b5.p <- data.frame(b5.p)
colnames(b5.p) <- c("lon", "lat", "val")

# SST
sst.small <- crop(sst, extent(37.45, 38.2, -47.2, -46.8))
sst.small.p <- rasterToPoints(sst.small)
sst.small.p <- data.frame(sst.small.p)
colnames(sst.small.p) <- c("lon", "lat", "val")

# Only near trips
dat.short <- tr[tr$Maxdist < 150, ]

# Bathymetry
# p3 <- ggplot(data = b5.p, aes(x = lon, y = lat))
# p3 <- p3 + geom_tile(aes(fill = val))
# p3 <- p3 + scale_fill_gradient(low = "gray60", high = "gray99", name = "Depth (m)", limits = c(-7000, 0)) #if land is NA
# p3 <- p3 + geom_polygon(data = island.df, aes(x = long, y = lat, group = group), fill = "grey")
# p3 <- p3 + geom_path(data = island.df, aes(x = long, y = lat, group = group), colour = "black")
# p3 <- p3 + geom_point(data = dat.short, aes(x = dat.short$Longitude, y = dat.short$Latitude, colour = Species), size = 0.6) + facet_wrap(~Species, nrow = 2) +
#   scale_colour_manual(values = c("#4daf4a", "#984ea3"))
# p3 <- p3 + coord_quickmap()
# p3 <- p3 + scale_x_continuous(expand = c(0,0), limits = c(37.45, 38.2))
# p3 <- p3 + scale_y_continuous(expand = c(0,0), limits = c(-47.2, -46.8))
# p3 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),
#                               panel.grid.major = element_blank(),
#                               strip.background = element_blank(),
#                               panel.border = element_rect(colour = "black"),
#                               axis.title.x = element_blank(),
#                               axis.title.y = element_blank())
# p3

# SST
p3 <- ggplot(data = sst.small.p, aes(x = lon, y = lat))
p3 <- p3 + geom_tile(aes(fill = val))
p3 <- p3 + scale_fill_gradientn(colours=parula(100),
                                guide = "colourbar",
                                name = "SST",
                                limits = c(-2, 24))
p3 <- p3 + geom_polygon(data = island.df, aes(x = long, y = lat, group = group), fill = "grey")
p3 <- p3 + geom_path(data = island.df, aes(x = long, y = lat, group = group), colour = "black")
p3 <- p3 + geom_point(data = dat.short, aes(x = dat.short$lon, y = dat.short$lat, colour = sp_code), size = 0.6) +
  facet_wrap(~sp_code, nrow = 2) +
  # scale_colour_manual(values = c("#4daf4a", "#984ea3"), guide = FALSE)
  scale_colour_manual(values = c("black", "black"), guide = FALSE)
p3 <- p3 + coord_quickmap()
p3 <- p3 + labs(title = "b")
p3 <- p3 + scale_x_continuous(expand = c(0,0), limits = c(37.45, 38.2))
p3 <- p3 + scale_y_continuous(expand = c(0,0), limits = c(-47.2, -46.805))
p3 <- p3 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              axis.text = element_text(colour = "black"),
                              panel.border = element_rect(colour = "black"),
                              plot.title = element_text(size = rel(1)),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank())
# p3

# grid.arrange(p2, p3, ncol = 2)

# Plot to file

# Distant trips
pdf("./Plots/mapDistant.pdf",
    width = double.col/fig.scale, height = (double.col/1.6)/fig.scale,
    useDingbats = FALSE)
print(p2)
dev.off()

# Nearby trips
pdf("./Plots/mapNearby.pdf",
    width = double.col/fig.scale, height = (double.col/1.6)/fig.scale,
    useDingbats = FALSE)
print(p3)
dev.off()

# Combined
pdf("./Plots/mapCombined.pdf",
    width = double.col/fig.scale, height = double.col/fig.scale,
    useDingbats = FALSE)
cowplot::plot_grid(p2, p3, axis = "t", ncol = 2, rel_widths = c(2, 1.5))
dev.off()

tiff("./Plots/mapCombined.tiff",
     width = double.col/fig.scale, height = double.col/fig.scale,
     res = 600,
     units = "in")
cowplot::plot_grid(p2, p3, axis = "t", ncol = 2, rel_widths = c(2, 1.5))
dev.off()

# Remove dummy variable
tr$Species2 <- NULL

#-----------------------------------------------------
#Visualization of environmental covariates
#-----------------------------------------------------

# Lat & Lon
library(beanplot)

# Latitude
beanplot(lat ~ sp_code, data = tr[tr$sex == "Male", ],
         main = NULL, ylab = "Latitude", side = "both",
         ylim = c(min(tr[tr$sex == "Male", "lat"], na.rm = T), max(tr[tr$sex == "Male", "lat"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

beanplot(lat ~ sp_code, data = tr[tr$sex == "Female", ],
         main = NULL, ylab = "Latitude", side = "both",
         ylim = c(min(tr[tr$sex == "Female", "lat"]), max(tr[tr$sex == "Female", "lat"])),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

# SST
beanplot(SST ~ Species, data = tr3[tr3$sex == "Male", ],
         main = NULL, ylab = "SST", side = "both",
         ylim = c(min(tr3[tr3$sex == "Male", "SST"], na.rm = T), max(tr3[tr3$sex == "Male", "SST"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

beanplot(SST ~ Species, data = tr3[tr3$sex == "Female", ],
         main = NULL, ylab = "SST", side = "both",
         ylim = c(min(tr3[tr3$sex == "Female", "SST"], na.rm = T), max(tr3[tr3$sex == "Female", "SST"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

# Depth
beanplot(DEPTH ~ Species, data = tr3[tr3$sex == "Male", ],
         main = NULL, ylab = "DEPTH", side = "both",
         ylim = c(min(tr3[tr3$sex == "Male", "DEPTH"], na.rm = T), max(tr3[tr3$sex == "Male", "DEPTH"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

beanplot(DEPTH ~ Species, data = tr3[tr3$sex == "Female", ],
         main = NULL, ylab = "DEPTH", side = "both",
         ylim = c(min(tr3[tr3$sex == "Female", "DEPTH"], na.rm = T), max(tr3[tr3$sex == "Female", "DEPTH"], na.rm = T)),
         bw = "nrd",
         overallline = "median",
         what = c(1,1,1,0), frame.plot = F,
         border = NA, col = list("#4daf4a", "#984ea3"))

#-----------------------------------
# Try split violins with GGPlot
source("./Scripts/function_splitViolin.R")
library(tidyr)

datV <- tr
datV$CHL <- log10(datV$CHL)

datV <- gather(data = datV,
               key = "Covariate",
               value = "Value",
               "lat",
               "SST",
               "CHL",
               "DEPTH2")

datV[datV$Covariate == "DEPTH2", "Covariate"] <- "DEPTH"

# Get in the correct order
datV$Covariate <- factor(datV$Covariate, levels = c("lat", "SST", "CHL", "DEPTH"))

pdf("./Plots/envarDensity.pdf",
    width = (single.col*1.6/fig.scale),
    height = (single.col*1.6/fig.scale),
    useDingbats = FALSE)
beans <- ggplot(data = datV, aes(x = sex, y = Value, fill = sp_code)) +
  geom_split_violin() +
  scale_fill_manual(values = c("#4daf4a", "#984ea3"), name = "Species") +
  facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Sex", y = "Covariate value") +
  theme_rr()
print(beans)
dev.off()


#-----------------------------------
## Random forests

library(randomForest)

dat.f <- tr[ , c("individual_id", "sp_code", "sex", "Distance", "DEPTH2", "SST",
                      "SSHA", "EKE", "CHL", "WINDU", "WINDV")]
# Select only complete cases
dat.f <- dat.f[complete.cases(dat.f), ]

# Log Chlorophyll
dat.f$CHL <- log10(dat.f$CHL)

# Only at-sea
dat.f <- dat.f[dat.f$DEPTH2 <= 0, ]

# First check for colinearity
cor(x = dat.f[ , c(4:11)], method = "spearman")

library(corrplot)
cr <- cor(x = dat.f[ , c(4:16)], method = "pearson")
corrplot(cr, method = "number", type = "lower")

library(GGally)
ggpairs(dat.f[ , c(4:16)])

# Leave out Distance, DIST.200, PROD, EKE

# Add species abreviations

dat.f$Group <- rep(NA, nrow(dat.f))
for (i in 1:nrow(dat.f)) {
  dat.f$Group[i] <- paste0(dat.f$sp_code[i], " ", dat.f$sex[i])
}

# Build RF
rf <- randomForest(data = dat.f, as.factor(Group) ~ DEPTH2 + SST +
                     SSHA + EKE + CHL + WINDU + WINDV, ntree = 1000,
                   na.action = na.omit, proximity = T, importance = F)
rf

varImpPlot(rf, pch = 16)

# Plots
pdf("./Plots/variableImportance.pdf",
    width = single.col/fig.scale,
    height = single.col/fig.scale,
    useDingbats = FALSE)
# par(ps = 9) #set font size
varImpPlot(rf, pch = 16)
dev.off()


plot(dat.f$SST, dat.f$CHL, col = as.factor(dat.f$Group)) # Plot two most important vars



pdf("./Plots/envars2D.pdf",
    width = single.col/fig.scale,
    height = single.col/fig.scale,
    useDingbats = FALSE)
p <- ggplot(dat.f, aes(x = SST, y = CHL, col = Group)) +
  geom_point() +
  scale_colour_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a")) +
  theme_rr() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
print(p)
dev.off()

MDSplot(rf, fac = as.factor(dat.f$Group), k = 5)

# MDS with rfPermute
library(rfPermute)
proximityPlot(rf, legend.loc = "right", circle.size = NULL)


# Make own plot - proximity as MDS
rf.mds <- cmdscale(1 - rf$proximity, k = 2)
mdsframe <- cbind.data.frame(dat.f$Group, rf.mds[ , 1], rf.mds[ ,2])
names(mdsframe) <- c("Group", "Dim1", "Dim2")

p.mds <- ggplot(mdsframe, aes(x = Dim1, y = Dim2, colour = Group))
p.mds <- p.mds + geom_point(size = 1.5, alpha = 1) + scale_colour_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"))
p.mds <- p.mds + theme_bw() + theme(panel.grid.minor = element_blank(),
                                    panel.grid.major = element_blank())
p.mds <- p.mds + labs(x = "Dimension 2", y = "Dimension 1")
p.mds

# 3D plots
library(plot3D)
Groups <- unique(dat.f$Group)
dat.f$g <- rep(NA, nrow(dat.f))
dat.f[dat.f$Group == Groups[1], "g"] <- 1
dat.f[dat.f$Group == Groups[2], "g"] <- 2
dat.f[dat.f$Group == Groups[3], "g"] <- 3
dat.f[dat.f$Group == Groups[4], "g"] <- 4

scatter3D(x = dat.f$DEPTH2, y = dat.f$SST, z = dat.f$CHL, colvar = dat.f$g,
          col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"),
          bty = "f",
          pch = 20,
          phi = 30, theta = 50,
          ticktype = "detailed",
          xlab = "DEPTH", ylab = "SST", zlab = "log(CHL)",
          colkey = list(at = c(1,2,3,4),labels = c("NGP Male", "NGP Female", "SGP Male", "SGP Female")))