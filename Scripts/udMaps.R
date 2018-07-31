#### Maps of giant petrel utilization distributions

## Ryan Reisinger

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

library(raster)
library(viridis)
library(ggplot2)
library(orsifronts)

# UDs

vud.males.ngp <- raster("./Output/vud.males.raster.ngp.grd")
vud.males.sgp <- raster("./Output/vud.males.raster.sgp.grd")
vud.females.ngp <- raster("./Output/vud.females.raster.ngp.grd")
vud.females.sgp <- raster("./Output/vud.females.raster.sgp.grd")

# NGP Males
ngpm <- rasterToPoints(vud.males.ngp)
ngpm <- data.frame(ngpm)
colnames(ngpm) <- c("lon", "lat", "val")
ngpm$sp <- "Northern giant petrel"
ngpm$sex <- "Male"

# NGP Females
ngpf <- rasterToPoints(vud.females.ngp)
ngpf <- data.frame(ngpf)
colnames(ngpf) <- c("lon", "lat", "val")
ngpf$sp <- "Northern giant petrel"
ngpf$sex <- "Female"

# SGP Males
sgpm <- rasterToPoints(vud.males.sgp)
sgpm <- data.frame(sgpm)
colnames(sgpm) <- c("lon", "lat", "val")
sgpm$sp <- "Southern giant petrel"
sgpm$sex <- "Male"

# SGP Females
sgpf <- rasterToPoints(vud.females.sgp)
sgpf <- data.frame(sgpf)
colnames(sgpf) <- c("lon", "lat", "val")
sgpf$sp <- "Southern giant petrel"
sgpf$sex <- "Female"

# Bind
kern <- rbind(ngpf, ngpm, sgpf, sgpm)
kern <- kern[kern$val <= 95, ]

# Add the Orsifronts
myOrsi <- fortify(orsifronts)
myOrsi <- myOrsi[myOrsi$id %in% c("stf", "saf", "pf"), ]

# Map
minx <- 8
maxx <- 60
miny <- -64
maxy <- -30

# Plot
p2 <- ggplot(data = kern, aes(x = lon, y = lat))
p2 <- p2 + geom_tile(aes(fill = val)) + facet_wrap(sp ~ sex)
p2 <- p2 + scale_fill_viridis(direction = 1, option = "plasma")
p2 <- p2 + geom_path(data = myOrsi, aes(x = long, y = lat, group = group), linetype = 2)
p2 <- p2 + coord_quickmap()
p2 <- p2 + scale_x_continuous(expand = c(0,0), limits = c(minx, maxx))
p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(miny, maxy))
p2 <- p2 + theme_bw() + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank())
world <- borders("world", colour = "black", fill = "grey", xlim = c(minx, maxx), ylim = c(miny, maxy))
p2 <- p2 + world
p2