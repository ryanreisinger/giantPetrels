# Stable isotopes

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

library(ggplot2)
library(dplyr)
library(data.table)

met <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)
dat <- read.csv("./Data/GP_sia_2019-07-04.csv", stringsAsFactors = F)

## Get sex
dat$sex <- "Unknown"
dat[!is.na(dat$Culmen_Length) & dat$Culmen_Length < 97, "sex"] <- "Female"
dat[!is.na(dat$Culmen_Length) & dat$Culmen_Length > 97, "sex"] <- "Male"

## Initial look
ggplot(data = dat, aes(x = d13C, y = d15N, colour = Species)) +
  geom_point() +
  facet_wrap(~Sample_Type, ncol = 1)

# ----------------------------------
## Wide table for easier plotting
d.temp <- dcast(data = setDT(dat), Individual_id + Species + sex ~ Sample_Type, value.var = c("d13C", "d15N"))

## Match tracking info
d <- merge(x = d.temp, y = met, by.x = "Individual_id", by.y = "individual_id",
           all.x = T, all.y = T)
rm(d.temp)

## Update sex
d[is.na(d$sex.y), "sex.y"] <- d[is.na(d$sex.y), "sex.x"]
d$sex <- d$sex.y
d$sex.x <- NULL
d$sex.y <- NULL

# Update 'Species' field
d[d$scientific_name == "Northern Giant Petrel", "scientific_name"] <- "NGP"
d[d$scientific_name == "Southern Giant Petrel", "scientific_name"] <- "SGP"
d[is.na(d$sex), "sex"] <- d[is.na(d$sex), "scientific_name"]

# Write summary
write.csv(d, "./Output/SummaryTable_track_sia.csv", row.names = F)
saveRDS(d, "./Output/Summary_track_sia.RDS")

# ----------------------------------
## Summary stats
# ----------------------------------

d.group <- group_by(d, Species, sex)

d.summary <- summarise(d.group,
                       d13C_RBC_ave = mean(d13C_RBC, na.rm = T),
                       d13C_Plas_ave = mean(d13C_Plas, na.rm = T),
                       d13C_DP_ave = mean(d13C_DP, na.rm = T),
                       d15N_RBC_ave = mean(d15N_RBC, na.rm = T),
                       d15N_Plas_ave = mean(d15N_Plas, na.rm = T),
                       d15N_DP_ave = mean(d15N_DP, na.rm = T),
                       d13C_RBC_sd = sd(d13C_RBC, na.rm = T),
                       d13C_Plas_sd = sd(d13C_Plas, na.rm = T),
                       d13C_DP_sd = sd(d13C_DP, na.rm = T),
                       d15N_RBC_sd = sd(d15N_RBC, na.rm = T),
                       d15N_Plas_sd = sd(d15N_Plas, na.rm = T),
                       d15N_DP_sd = sd(d15N_DP, na.rm = T)
)
# ----------------------------------

# ----------------------------------
## Some exploratory plots

# ------------------
## Plasma v. RBC
# ------------------

# ------------------
## d13C
# Correlation
tmp <- cor(d$d13C_DP[!is.na(d$d13C_DP) & !is.na(d$d13C_RBC)],
    d$d13C_RBC[!is.na(d$d13C_DP) & !is.na(d$d13C_RBC)])

ggplot(data = d, aes(x = d13C_DP, y = d13C_RBC, colour = Species, shape = sex)) +
  geom_point() +
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, lty = 1) +
  # geom_abline(slope = 1, intercept = -1, lty = 2) +
  # geom_abline(slope = 1, intercept = +1, lty = 2) +
  geom_smooth(method = "lm", data = d, aes(x = d13C_DP, y = d13C_RBC), inherit.aes = F,
              se = FALSE) +
  labs(title = "d13C: delipated plasma v. red blood cells",
       subtitle = paste0("r = ", round(tmp, 2)))

# ------------------
## d15N
# Correlation
tmp <- cor(d$d15N_Plas[!is.na(d$d15N_Plas) & !is.na(d$d15N_RBC)],
    d$d15N_RBC[!is.na(d$d15N_Plas) & !is.na(d$d15N_RBC)])

ggplot(data = d, aes(x = d15N_Plas, y = d15N_RBC, colour = Species, shape = sex)) +
  geom_point() +
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, lty = 1) +
  # geom_abline(slope = 1, intercept = -1, lty = 2) +
  # geom_abline(slope = 1, intercept = +1, lty = 2) +
  geom_smooth(method = "lm", data = d, aes(x = d15N_Plas, y = d15N_RBC), inherit.aes = F,
              se = FALSE)  +
  labs(title = "d15N: plasma v. red blood cells",
       subtitle = paste0("r = ", round(tmp, 2)))

# ------------------
## Isospace
# ------------------

## Plasma
ggplot(data = d, aes(x = d13C_DP, y = d15N_Plas, colour = Species, shape = sex)) +
  geom_point() +
  facet_wrap(~long.trip, ncol = 1) +
  coord_equal()

## RBC
ggplot(data = d, aes(x = d13C_RBC, y = d15N_RBC, colour = Species, shape = sex)) +
  geom_point() +
  facet_wrap(~long.trip, ncol = 1) +
  coord_equal()

# ------------------
## Latitude v. SIA values
# ------------------

## Only for long-trippers

# ------------------
## Plasma
# Correlation
temp <- cor(d$d13C_DP[!is.na(d$d13C_DP) & !is.na(d$mean.lat)],
    d$mean.lat[!is.na(d$d13C_DP) & !is.na(d$mean.lat)])

ggplot(data = d[d$long.trip == "Y", ], aes(x = d13C_DP, y = mean.lat, colour = Species, shape = sex)) +
  geom_point() +
  geom_smooth(method = "lm", data = d, aes(x = d13C_DP, y = mean.lat),
              inherit.aes = F,
              se = FALSE) +
  labs(title = "Latitude v. d13C delipidated plasma",
       subtitle = paste0("r = ", round(temp, 2)))

# ------------------
## RBC
# Correlation
temp <- cor(d$d13C_RBC[!is.na(d$d13C_RBC) & !is.na(d$mean.lat)],
    d$mean.lat[!is.na(d$d13C_RBC) & !is.na(d$mean.lat)])

ggplot(data = d[d$long.trip == "Y", ], aes(x = d13C_RBC, y = mean.lat, colour = Species, shape = sex)) +
  geom_point() +
  geom_smooth(method = "lm", data = d, aes(x = d13C_RBC, y = mean.lat),
              inherit.aes = F,
              se = FALSE) +
  labs(title = "Latitude v. d13C red blood cells",
       subtitle = paste0("r = ", round(temp, 2)))


# ------------------
## d13C v. d15N
# ------------------

# ------------------
# RBC
temp <- cor(d$d13C_RBC[!is.na(d$d13C_RBC) & !is.na(d$d15N_RBC)],
            d$d15N_RBC[!is.na(d$d13C_RBC) & !is.na(d$d15N_RBC)])

ggplot(data = d, aes(x = d13C_RBC, y = d15N_RBC, colour = Species, shape = sex)) +
  geom_point() +
  geom_smooth(method = "lm", data = d, aes(x = d13C_RBC, y = d15N_RBC),
              inherit.aes = F,
              se = FALSE) +
  labs(title = "d13C v. d15N red blood cells",
       subtitle = paste0("r = ", round(temp, 2)))

# ------------------
# Plasma
temp <- cor(d$d13C_DP[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas)],
            d$d15N_Plas[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas)])

ggplot(data = d, aes(x = d13C_DP, y = d15N_Plas, colour = Species, shape = sex)) +
  geom_point() +
  geom_smooth(method = "lm", data = d, aes(x = d13C_DP, y = d15N_Plas),
              inherit.aes = F,
              se = FALSE) +
  labs(title = "d13C v. d15N plasma",
       subtitle = paste0("r = ", round(temp, 2)))

# ------------------
# Look at the residuals

#---------------------------------
#---------------------------------
# Test residual values
library(MASS)
library(ggbeeswarm)

d$groupname <- paste0(d$Species, " - ", d$sex)

# -----------------------------
# d13C response - RBC
d1 <- d[!is.na(d$d13C_RBC) & !is.na(d15N_RBC), ]
m1 <- lm(formula = d13C_RBC ~ d15N_RBC, data = d1)

d1$d13C_RBC_residuals <- m1$residuals
d1$d13C_RBC_studentized <- studres(m1)

ggplot(d1, aes(x = groupname, y = d13C_RBC_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "", y = "d13C Studentised residuals")


# -----------------------------
# d15N response - RBC
# d2 <- d[!is.na(d$d13C_RBC) & !is.na(d15N_RBC) & d$long.trip == "N", ]
d2 <- d[!is.na(d$d13C_RBC) & !is.na(d15N_RBC), ]
m2 <- lm(formula = d15N_RBC ~ d13C_RBC, data = d2)

d2$d15N_RBC_residuals <- m2$residuals
d2$d15N_RBC_studentized <- studres(m2)

ggplot(d2, aes(x = groupname, y = d15N_RBC_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "", y = "d15N Studentised residuals")

# -----------------------------
# d15N response - plasma
# d2 <- d[!is.na(d$d13C_DP) & !is.na(d15N_Plas) & d$long.trip == "N", ]
d2 <- d[!is.na(d$d13C_DP) & !is.na(d15N_Plas), ]
m2 <- lm(formula = d15N_Plas ~ d13C_DP, data = d2)

d2$d15N_Plas_residuals <- m2$residuals
d2$d15N_Plas_studentized <- studres(m2)

ggplot(d2, aes(x = groupname, y = d15N_Plas_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "", y = "d15N Studentised residuals")
