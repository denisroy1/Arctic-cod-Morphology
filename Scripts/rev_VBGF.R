#### rev_VBGF.R ####
## Script reads in raw Arctic cod data and uses the von Bertalanffy growth function (VBGF)
## to estimate age from length (in mm). We have data from Arctic cod collected in the 
## Beaufort Sea in 2018-2019 by our collaborators at DFO, and want to assess their 
## traditional and geometric morphometrics to determine whether there are differences
## among habitat aggregations (as determined by Majewski et al. 2017) and/or age classes.

## As the fish have not been aged we need to sort them into age classes before we can proceed with the age class 
## analyses. To do this we need to use the VBGF parameter for Arctic cod and these were recently published 
## by Forster et al. 2020 Deep-Sea Research Part II doi 10.1016/j.dsr2.2020.104779.
## Here we use the length determined for our fish and reverse the VBGF to estimate age. 

## Written by jm and dr May 2022, updated June 2023.

#### Set up ####

## clearing all instances (data and variables)
rm(list = ls())

## loading the libraries needed 
{
  library(geomorph)
  library(ggplot2)
  library(colorspace)
  library(mosaic)
  library(emmeans)
}

## Setting the working directory to where data and scripts are located
setwd("/Users/denis/Library/CloudStorage/Dropbox/McGill/students/Honours/Juliano Mazalia/Data/manuscript/ArcticScience/Reviews/resubmission/2ndsubmission/2nd review/working replies/2023-03-31/3rd submission/4th submission/Data and Scripts/")

## Reading in the data from the file. Looking for "Arctic_cod_raw.csv" file.
acdata <- read.csv("Arctic_cod_raw.csv", header=TRUE, stringsAsFactors = T)

## View the incoming data
View(acdata)

#### Arctic cod VBGF parameters ####
## Use the Arctic cod parameters from Forster et al. 2020 and set these as
## variables.
Linf_ac <- 290.7462
K_ac <- 0.2175
t0_ac <- -0.8603

## Reversed Function for age estimates for given lengths in mm.
VBage_est <- ((log(1 - (as.numeric(acdata$length_mm) / Linf_ac))) / -K_ac) + t0_ac

## Round the VBage_est to the nearest integer to create year classes.
VBage <- ceiling(VBage_est)

## Place VBage data at the end of the acdata dataframe.
acdata$VBage <- as.factor(VBage)

## Run a basic linear model to assess the relationship between length and VBage:
agelen <- lm(length_mm ~ as.factor(VBage), data = acdata)

## Use results of lm above to generate the means, SEs, and 95% CIs for the means 
## of each age class.
agelenmat <- as.data.frame(emmeans(agelen, "VBage",type = "response"))

## Potting the length ~ age relationship from the VBGF for Arctic cod using our own 
## data.
ggplot(agelenmat, aes(x = emmean, y = as.factor(VBage))) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.5, linewidth = 2) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(acdata, mapping = aes(x = length_mm, y = VBage, color = VBage), 
              position = position_jitter(0.3), size = 5, stroke = 1, alpha = 0.5) +
  labs(x = "Length (mm)", y = "VB age", color = "VBage") +
  scale_color_manual(values = c("firebrick4", "#E69F00", "#56B4E9", "#009E73","#F0E442", "burlywood4")) +
  scale_x_continuous(limits = c(0, 210), breaks = seq(0, 210, by = 50)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 35),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")

## Filter out the individuals useful for GM and TM analyses according 
## to picture (=True), and picture quality (=g, for good).
acdata <- acdata[which(acdata$GM == "TRUE" & acdata$picqual == "g"),]

## Reorganise the acdata dataframe to something that can be usable in the other scripts.
acdata <- cbind.data.frame(year = acdata$year, station = acdata$station, fishID = acdata$fishID, 
                           depth = acdata$depth_m, assem_depth = acdata$assemblage_depth, 
                           assem_shore = acdata$assemblage_shore, dist_to_shore = acdata$shoredistance_km, 
                           length_mm = acdata$length_mm, VBage = acdata$VBage)

## Options to save Arctic cod data with age classes assessed
write.csv(acdata,"ACmeta1.csv", row.names = F)

#### END ####
