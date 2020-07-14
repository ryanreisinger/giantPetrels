## Create individual summary maps

library(ggplot2)
library(gridExtra)
library(rJava)
library(OpenStreetMap)

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

tracks <- readRDS("./Output/tracks_trips.RDS")
d <- readRDS("./Output/Summary_track_sia.RDS")

ids <- unique(tracks$individual_id)

pdf("./Output/trackmaps.pdf",
    paper = "a4r")

for (i in 1:length(ids)) {
  this.id <- ids[i]
  if(!is.na(this.id)) {
  print(this.id)
  this.track <- tracks[tracks$individual_id == this.id, ]
  this.met <- d[d$Individual_id == this.id, ]
  
  # Set up limits
  yrange <- diff(range(this.track$decimal_latitude, na.rm = T))*0.1
  xrange <- diff(range(this.track$decimal_longitude, na.rm = T))*0.1
  
  topy <- max(this.track$decimal_latitude, na.rm = T) + yrange
  topx <- min(this.track$decimal_longitude, na.rm = T) - xrange
  bottomy <- min(this.track$decimal_latitude, na.rm = T) - yrange
  bottomx <- max(this.track$decimal_longitude, na.rm = T) + xrange
  
  # Create a table for display
  this.table <- this.met[, c("d13C_Plas", "d13C_DP", "d13C_RBC",
                             "d15N_Plas", "d15N_DP", "d15N_RBC",
                             "Culmen_length", "Maxdist")]
  this.table$Maxdist = round(this.table$Maxdist, 1)
  
  # Get map data
  map <- openmap(c(topy,
                   topx),
                 c(bottomy,
                   bottomx),
                 type="bing")
  
  # Project
  map <- openproj(map, projection = "+proj=longlat +datum=WGS84 +no_defs")
  
  # Plot
  p1 <- autoplot(map) + 
    geom_path(data = this.track, aes(x = decimal_longitude, y = decimal_latitude),
              colour = "grey") +
    geom_point(data = this.track, aes(x = decimal_longitude, y = decimal_latitude),
               colour = "red") +
    labs(title = this.met$Individual_id,
         subtitle = this.met$sex,
         x = "Longitude",
         y = "Latitude")
  p2 <- tableGrob(this.table,
                  theme = ttheme_minimal(base_size = 9))
  p3 <- grid.arrange(p1, p2, nrow = 2, heights = c(4, 1))
  
  print(p3)
  }
}

dev.off()
