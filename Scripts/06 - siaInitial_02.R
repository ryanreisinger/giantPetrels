# Stable isotopes

setwd("D:/PEI_Toppredators/Giant Petrels/Working/giantPetrels")

library(ggplot2)
library(dplyr)
library(data.table)

met <- read.csv("./Output/SummaryTable.csv", stringsAsFactors = F)
dat <- read.csv("./Data/GP_sia_2019-09-20.csv", stringsAsFactors = F)

# --------------------------------
## Plotting stuff

## Figure widths in mm
single.col <- 70*0.0393701
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
           all.x = T, all.y = T)

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
           y = met[, c("individual_id", "Culmen_length", "date.start", "duration", "Maxdist", "mean.lat", "long.trip", "short.trip")],
           by.x = "Individual_id",
           by.y = "individual_id",
           all.x = TRUE,
           all.y = FALSE)

# -----------------------------------------------
## Lipid normalization of d13C where there are no delipidated plasma samples

# Which are these? (13 values)
d[is.na(d$d13C_DP) & !is.na(d$d13C_Plas), ]

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
d[is.na(d$d13C_DP) & !is.na(d$d13C_Plas), "Normalized"] <- "Y"
d[is.na(d$d13C_Plas) & is.na(d$d15N_Plas), "Normalized"] <- NA

# Replace missing values with predictions
d[is.na(d$d13C_DP) & !is.na(d$d13C_Plas), "d13C_DP"] <- d[is.na(d$d13C_DP) & !is.na(d$d13C_Plas), "Predicted_d13C_DP"]

# Tidy up
rm(d13C_delta, CN_delta, m1, preds, foo)

# -----------------------------------------------
## How many birds with plasma samples?
nrow(d[!is.na(d$d13C_Plas), ])

# Remove birds where we still don't have sex
d <- filter(d, !is.na(sex))

## Write summary
write.csv(d, "./Output/SummaryTable_track_sia.csv", row.names = F)
saveRDS(d, "./Output/Summary_track_sia.RDS")

## Now, continue with only the plasma sampled birds
d <- filter(d, !is.na(d13C_Plas) & !is.na(d15N_Plas))

## The linear model is not great, so proceed without the normalized samples
d <- filter(d, Normalized == "N")

## How many birds for SIA analyses?
nrow(d)

## How many of these birds were tracked
nrow(d[!is.na(d$Maxdist), ])

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

# Write to file
write.csv(d.summary, "./Output/SIAgroupSummary.csv", row.names = F)
# ----------------------------------

#----------------------------------------------
## TESTS
#----------------------------------------------

#----------------------------------------------
## Test for multivariate normality
library(MVN)

mvn(data = d[d$sex == "Male" & d$Species == "NGP", c("d13C_DP", "d15N_Plas")], mvnTest = "mardia") #fine
mvn(data = d[d$sex == "Female" & d$Species == "NGP", c("d13C_DP", "d15N_Plas")], mvnTest = "mardia") #fine
mvn(data = d[d$sex == "Male" & d$Species == "SGP", c("d13C_DP", "d15N_Plas")], mvnTest = "mardia") #fine
mvn(data = d[d$sex == "Female" & d$Species == "SGP", c("d13C_DP", "d15N_Plas")], mvnTest = "mardia") #fine

#----------------------------------------------
## PERMANOVA/MANOVA

# # Permanova
# adonis(dat2[ , c("d13C", "d15N")] ~ group, data = dat2, method = "euclidean")

# Manova
dat2 <- d
dat2$group <- paste0(dat2$Species, "-", dat2$sex)
man <- manova(data = dat2, cbind(d13C_DP, d15N_Plas) ~ group)
summary(man, test = "Wilks")
summary.aov(man)

# Compare males
dat3 <- dat2[dat2$sex == "Male", ]
man <- manova(data = dat3, cbind(d13C_DP, d15N_Plas) ~ Species)
summary(man, test = "Wilks")

# Compare females
dat3 <- dat2[dat2$sex == "Female", ]
man <- manova(data = dat3, cbind(d13C_DP, d15N_Plas) ~ Species)
summary(man, test = "Wilks")

# Compare NGPs
dat3 <- dat2[dat2$Species == "NGP", ]
man <- manova(data = dat3, cbind(d13C_DP, d15N_Plas) ~ sex)
summary(man, test = "Wilks")

# Compare SGPs
dat3 <- dat2[dat2$Species == "SGP", ]
man <- manova(data = dat3, cbind(d13C_DP, d15N_Plas) ~ sex)
summary(man, test = "Wilks")

# Only inidivuals making short trips
dat3 <- dat2[!is.na(dat2$long.trip) & dat2$long.trip == "N", ]
man <- manova(data = dat3, cbind(d13C_DP, d15N_Plas) ~ Species)
summary(man, test = "Wilks")

# What are these values?
temp <- group_by(dat3, Species)
temp <- summarise(temp,
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

as.data.frame(temp)
rm(temp)

# Only individuals that made long trips
dat3 <- dat2[!is.na(dat2$long.trip) & dat2$long.trip == "Y", ]
man <- manova(data = dat3, cbind(d13C_DP, d15N_Plas) ~ Species)
summary(man, test = "Wilks")

# What are these values?
# But take care, with d13C/latitude influence on d15N here.
temp <- group_by(dat3, Species)
temp <- summarise(temp,
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

as.data.frame(temp)
rm(temp)


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

p <- ggplot(data = d, aes(x = d13C_DP, y = d13C_RBC, colour = interaction(Species, sex))) +
  geom_point() +
  coord_equal() +
  scale_colour_manual(values = c("#377eb8", "#984ea3", "#e41a1c", "#4daf4a")) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey") +
  geom_smooth(method = "lm", data = d, aes(x = d13C_DP, y = d13C_RBC), inherit.aes = F,
              se = FALSE,
              colour = "black") +
  labs(title = "a. d13C: delipated plasma v. red blood cells",
       subtitle = paste0("r = ", round(tmp, 2)),
       x = "d13C delipidated plasma",
       y = "d13C red blood cells") +
  theme_rr()

p

pdf("./Plots/plasmaRBC_d13C.pdf",
    width = single.col/fig.scale + 2, height = single.col/fig.scale,
    useDingbats = FALSE)
print(p)
dev.off()

# ------------------
## d15N
# Correlation
tmp <- cor(d$d15N_Plas[!is.na(d$d15N_Plas) & !is.na(d$d15N_RBC)],
    d$d15N_RBC[!is.na(d$d15N_Plas) & !is.na(d$d15N_RBC)])

p <- ggplot(data = d, aes(x = d15N_Plas, y = d15N_RBC, colour = interaction(Species, sex))) +
  geom_point() +
  coord_equal() +
  scale_colour_manual(values = c("#377eb8", "#984ea3", "#e41a1c", "#4daf4a")) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey") +
  geom_smooth(method = "lm", data = d, aes(x = d15N_Plas, y = d15N_RBC), inherit.aes = F,
              se = FALSE,
              colour = "black") +
  labs(title = "b. d15N: plasma v. red blood cells",
       subtitle = paste0("r = ", round(tmp, 2)),
       x = "d15N plasma",
       y = "d15N red blood cells") +
  theme_rr()

p

pdf("./Plots/plasmaRBC_d15N.pdf",
    width = single.col/fig.scale + 2, height = single.col/fig.scale,
    useDingbats = FALSE)
print(p)
dev.off()

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
# Correlation in d13C
temp <- cor(d$d13C_DP[!is.na(d$d13C_DP) & !is.na(d$mean.lat) & d$Normalized == "N"],
    d$mean.lat[!is.na(d$d13C_DP) & !is.na(d$mean.lat) & d$Normalized == "N"],
    method = "pearson")

ggplot(data = d[d$long.trip == "Y", ], aes(x = d13C_DP, y = mean.lat, colour = Species, shape = sex)) +
  geom_point() +
  geom_smooth(method = "lm", data = d, aes(x = d13C_DP, y = mean.lat),
              inherit.aes = F,
              se = FALSE) +
  labs(title = "Latitude v. d13C delipidated plasma",
       subtitle = paste0("r = ", round(temp, 2)))

# Correlation in d15N
temp <- cor(d$d15N_Plas[!is.na(d$d15N_Plas) & !is.na(d$mean.lat) & d$Normalized == "N"],
            d$mean.lat[!is.na(d$d15N_Plas) & !is.na(d$mean.lat) & d$Normalized == "N"],
            method = "pearson")

ggplot(data = d[d$long.trip == "Y", ], aes(x = d15N_Plas, y = mean.lat, colour = Species, shape = sex)) +
  geom_point() +
  geom_smooth(method = "lm", data = d, aes(x = d15N_Plas, y = mean.lat),
              inherit.aes = F,
              se = FALSE) +
  labs(title = "Latitude v. d15N plasma",
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
  # geom_smooth(method = "lm", data = d, aes(x = d13C_RBC, y = d15N_RBC),
  #             inherit.aes = F,
  #             se = FALSE) +
  labs(title = "d13C v. d15N red blood cells",
       subtitle = paste0("r = ", round(temp, 2)))

## Only long trips
ggplot(data = d[d$long.trip == "Y", ], aes(x = d13C_RBC, y = d15N_RBC, colour = Species, shape = sex)) +
  geom_point() +
  labs(title = "d13C v. d15N red blood cells")

# ------------------
# Plasma
temp <- cor(d$d13C_DP[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas)],
            d$d15N_Plas[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas)])

ggplot(data = d, aes(x = d13C_DP, y = d15N_Plas, colour = Species, shape = sex)) +
  geom_point() +
  # geom_smooth(method = "lm", data = d, aes(x = d13C_DP, y = d15N_Plas),
  #             inherit.aes = F,
  #             se = FALSE) +
  labs(title = "d13C v. d15N plasma",
       subtitle = paste0("r = ", round(temp, 2)))

# Only long trips
ggplot(data = d[d$long.trip == "Y", ], aes(x = d13C_DP, y = d15N_Plas, colour = Species, shape = sex)) +
  geom_point() +
  labs(title = "d13C v. d15N plasma")

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

#Long trips only?
d1 <- d[!is.na(d$long.trip), ]
d1 <- d1[!is.na(d1$d13C_RBC) & !is.na(d1$d15N_RBC) & d1$long.trip == "Y", ]

# d1 <- d[!is.na(d$d13C_RBC) & !is.na(d$d15N_RBC), ]
m1 <- lm(formula = d13C_RBC ~ d15N_RBC, data = d1)

d1$d13C_RBC_residuals <- m1$residuals
d1$d13C_RBC_studentized <- studres(m1)

ggplot(d1, aes(x = groupname, y = d13C_RBC_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "d13C Residuals in red blood cells",
       x = "", y = "d13C Studentised residuals")

# -----------------------------
# d13C response - Plasma

#Long trips only?
d1 <- d[!is.na(d$long.trip), ]
d1 <- d1[!is.na(d1$d13C_DP) & !is.na(d1$d15N_Plas) & d1$long.trip == "N", ]

# d1 <- d[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas), ]
m1 <- lm(formula = d13C_DP ~ d15N_Plas, data = d1)

d1$d13C_DP_residuals <- m1$residuals
d1$d13C_DP_studentized <- studres(m1)

ggplot(d1, aes(x = groupname, y = d13C_DP_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "d13C Residuals in plasma",
       x = "", y = "d13C Studentised residuals")

# -----------------------------
# d15N response - RBC

#Long trips only?
d2 <- d[!is.na(d$long.trip), ]
d2 <- d2[!is.na(d2$d13C_RBC) & !is.na(d2$d15N_RBC) & d2$long.trip == "Y", ]

# d2 <- d[!is.na(d$d13C_RBC) & !is.na(d$d15N_RBC), ]
m2 <- lm(formula = d15N_RBC ~ d13C_RBC, data = d2)

d2$d15N_RBC_residuals <- m2$residuals
d2$d15N_RBC_studentized <- studres(m2)

ggplot(d2, aes(x = groupname, y = d15N_RBC_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "d15N Residuals in red blood cells",
       x = "", y = "d15N Studentised residuals")

# -----------------------------
# d15N response - plasma

#All
d2 <- d[!is.na(d$long.trip), ]
d2 <- d2[!is.na(d2$d13C_DP) & !is.na(d2$d15N_Plas), ]

# d2 <- d[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas), ]
m2 <- lm(formula = d15N_Plas ~ d13C_DP, data = d2)

d2$d15N_Plas_residuals <- m2$residuals
d2$d15N_Plas_studentized <- studres(m2)

p <- ggplot(d2, aes(x = groupname, y = d15N_Plas_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  scale_color_manual(guide = F,
                     values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "",
       y = "d15N Studentised residuals") +
  theme_rr()

pdf("./Plots/d15NresidualsAllTrips.pdf",
    width = single.col/fig.scale, height = single.col/fig.scale,
    useDingbats = FALSE)
print(p)
dev.off()


#Long trips only?
d2 <- d[!is.na(d$long.trip), ]
d2 <- d2[!is.na(d2$d13C_DP) & !is.na(d2$d15N_Plas) & d2$long.trip == "Y", ]

# d2 <- d[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas), ]
m2 <- lm(formula = d15N_Plas ~ d13C_DP, data = d2)

d2$d15N_Plas_residuals <- m2$residuals
d2$d15N_Plas_studentized <- studres(m2)

p <- ggplot(d2, aes(x = groupname, y = d15N_Plas_studentized, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  scale_color_manual(guide = F,
                     values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(x = "",
       y = "d15N Studentised residuals") +
  theme_rr()

pdf("./Plots/d15NresidualsLongTrips.pdf",
    width = single.col/fig.scale, height = single.col/fig.scale,
    useDingbats = FALSE)
print(p)
dev.off()

#Short trips only?
d2 <- d[!is.na(d$long.trip), ]
d2 <- d2[!is.na(d2$d13C_DP) & !is.na(d2$d15N_Plas) & d2$long.trip == "N", ]

# d2 <- d[!is.na(d$d13C_DP) & !is.na(d$d15N_Plas), ]
m2 <- lm(formula = d15N_Plas ~ d13C_DP, data = d2)

d2$d15N_Plas_residuals <- m2$residuals
d2$d15N_Plas_studentized <- studres(m2)

p <- ggplot(d2, aes(x = groupname, y = d13C_DP, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  scale_color_manual(guide = F,
                     values = c("#e41a1c", "#4daf4a"), name = "") +
  labs(x = "",
       y = "d15N Studentised residuals") +
  theme_rr()

pdf("./Plots/d15NresidualsShortTrips.pdf",
    width = single.col/fig.scale, height = single.col/fig.scale,
    useDingbats = FALSE)
print(p)
dev.off()

# -----------------------------
# Short trip males
d3 <- d[!is.na(d$long.trip), ]
d3 <- d3[d3$long.trip == "N" & d3$sex == "Male", ]

# d15N Plasma
ggplot(d3, aes(x = groupname, y = d15N_Plas, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "plasma d15N in short trip males",
       x = "", y = "plasma d15N")

# d15N Red blood cells
ggplot(d3, aes(x = groupname, y = d15N_RBC, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "red blood cell d15N in short trip males",
       x = "", y = "red blood cell d15N")

# d13C Plasma
ggplot(d3, aes(x = groupname, y = d13C_DP, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "plasma d13C in short trip males",
       x = "", y = "delpidated plasma d13C")

# d13C Red blood cells
ggplot(d3, aes(x = groupname, y = d13C_RBC, colour = groupname)) +
  # geom_boxplot() +
  geom_boxplot(fill = grey(0.8), colour = grey(0.8), outlier.shape = NA, width = 0.25) +
  geom_quasirandom(width = 0.15) +
  # scale_color_manual(values = c("#377eb8", "#e41a1c", "#984ea3", "#4daf4a"), name = "") +
  labs(title = "red blood cell d13C in short trip males",
       x = "", y = "red blood cell d13C")
