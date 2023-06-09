#### VBGF test curve.R ####
## Generation of the age estimate based on VBGF curve in relation to Malizia et al.
## This is a script that we use to respond to a reviewer comment that would like to
## know how the rounding of lengths and/or VBGF generated age estimate could impact 
## our data.

## Written by dr May 29 2023

## We aim to demonstrate the shape of the VBGF as it pertains to Arctic cod in the 
## Beaufort Sea - as derived using the parameters listed in Forster et al. 2020
## Deep-Sea Research Part II doi 10.1016/j.dsr2.2020.104779 - is:
## a - adequate
## b - close to a gradual relationship and so not likely to have a big influence 
## on our age classes as derived in our paper.

#### set up ####

## clearing all instances (data and variables)
rm(list = ls())

## loading the libraries needed 
{
  library(geomorph)
  library(ggplot2)
  library(colorspace)
}

#### Importing and Cleaning the data ####

## Setting the working directory to where data and scripts are located
setwd("/pathtofile/")

## Reading in the data from the file. Looking for "cod_photo.csv" file.
ARCD_GM_data <- read.csv("Arctic_cod_raw.csv", header=TRUE)

## To make age estimates, we used the VBGF with parameters as determined from the 
## the Forster et al. 2020 paper (see reference above). 

#### Mean and sdev of data ####
## First calculate the mean and standard deviation of the length data 
## without including the NAs (can't use these). 
ml<-mean(na.omit(ARCD_GM_data$length_mm))
sdl<-sd(na.omit(ARCD_GM_data$length_mm))

#### Length - Freq hist of data ####
## Second, plot the length frequency distribution of the overall data
## while fitting a normal distribution.
ggplot(ARCD_GM_data, aes(x=length_mm)) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 8, fill = "firebrick4", color = "black") +
  stat_function(fun = dnorm, args = list(mean = ml, sd = sdl), linewidth = 2, color = "steelblue3", linetype = 2) +
  labs(x="Length in mm", y="Frequency") +
  theme_classic()+
  theme(axis.title.x = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 15))

## While the data may not be fully normal, we will assume that it is close and that the 
## mean and standard deviation of the distribution do an approximately good job at 
## characterising the central tendencies of the data. 

## We will then use the mean and SD to generate a random distribution of lengths
## that we can then use to visualise the VBGF.

#### Generating data ####
## Generating a random normal distribution of length based on the mean and 
## SD from our actual data.
L <- rnorm(1000,ml, sdl)
hist(L, breaks = 10, col = "steelblue4")

## The idea here is that if the VBGF curve is steep then small changes in length (especially at 
## larger sizes can lead to large changes in estimated ages (and possibly introduce bias). 
## On the other hand, a more gradual curve means a small change in length is also likely 
## to translate to more consistent estimates of age.

## To demonstrate this, we compare the VBGF for Myoxocephalus scorpius (an Arctic sculpin)
## which has a rather steep VBGF and whose parameters are also available in Forster et al. 2020

#### Myoxocephalus scorpius ####
## Forster derived parameters
Linf_ms <- 164
K_ms <- 0.69
t0_ms <- -0.9

## Function for Age estimates and estimated age values for given lengths 
Age <- ((log(1 - (L / Linf_ms))) / -K_ms) + t0_ms
ms125 <- ((log(1 - (125 / Linf_ms))) / -K_ms) + t0_ms
ms126 <- ((log(1 - (126 / Linf_ms))) / -K_ms) + t0_ms
ms152 <- ((log(1 - (152 / Linf_ms))) / -K_ms) + t0_ms
ms153 <- ((log(1 - (153 / Linf_ms))) / -K_ms) + t0_ms

## Combining general length and age estimated from VBGF
msdata<-cbind.data.frame(length = L, age = Age)

## Potting the VBGF for Myoxocephalus scorpius. It has a very steep VBGF.
## While small changes at low values generate age estimates that are relatively
## easy to differentiate, small differences at larger lengths can translate 
## into very different age estimates.
ggplot(msdata, aes(x=length, y=age)) +
  geom_point(size=2.5, col="steelblue4") +
  ggtitle("Shorthorn sculpin") +
  geom_segment(aes(x = 125, y = 0, xend = 125, yend = ms125), 
               col = "springgreen4",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ms125, col = "springgreen4", lty = 2, linewidth = 0.25) +
  geom_segment(aes(x = 126, y = 0, xend = 126, yend = ms126), 
               col = "springgreen4",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ms126, col = "springgreen4", lty = 2, linewidth = 0.25) +
  geom_segment(aes(x = 152, y = 0, xend = 152, yend = ms152), 
               col = "red",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ms152, col = "red", lty = 2, linewidth = 0.25) +
  geom_segment(aes(x = 153, y = 0, xend = 153, yend = ms153), 
               col = "red",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ms153, col = "red", lty = 2, linewidth = 0.25) +
  xlab("Length (in mm)") +
  ylab("VBGF est. Age") +
  theme_classic() +
  xlim(0, 200) +
  ylim(0,5) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(color = "black", size = 35, face = "bold.italic"))

## Using the parameters from Forster for Arctic cod shows a much more gradual VBGF 
## than for Myoxocephalus scorpius. Thus small changes, even at longer lengths are 
## still relatively similar to one another thus likely minimising (though not 
## completely eliminating) age estimated biases or errors. 

#### Arctic cod ####
## Forster derived parameters
Linf_ac <- 290.7462
K_ac <- 0.2175
t0_ac <- -0.8603

## Function for Age estimates and estimated age values for given lengths
Age <- ((log(1 - (L / Linf_ac))) / -K_ac) + t0_ac
ac125 <- ((log(1 - (125 / Linf_ac))) / -K_ac) + t0_ac
ac126 <- ((log(1 - (126 / Linf_ac))) / -K_ac) + t0_ac
ac152 <- ((log(1 - (152 / Linf_ac))) / -K_ac) + t0_ac
ac153 <- ((log(1 - (153 / Linf_ac))) / -K_ac) + t0_ac

## Combining general length and age estimated from VBGF
acdata<-cbind.data.frame(length = L, age = Age)

## Potting the VBGF for Arctic cod.
ggplot(acdata, aes(x=length, y=age)) +
  geom_point(size=2.5, col="steelblue4") +
  ggtitle("Arctic cod") +
  geom_segment(aes(x = 125, y = 0, xend = 125, yend = ac125), 
               col = "springgreen4",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ac125, col = "springgreen4", lty = 2, linewidth = 0.25) +
  geom_segment(aes(x = 126, y = 0, xend = 126, yend = ac126), 
               col = "springgreen4",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ac126, col = "springgreen4", lty = 2, linewidth = 0.25) +
  geom_segment(aes(x = 152, y = 0, xend = 152, yend = ac152), 
               col = "red",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ac152, col = "red", lty = 2, linewidth = 0.25) +
  geom_segment(aes(x = 153, y = 0, xend = 153, yend = ac153), 
               col = "red",lty = 2, linewidth = 0.25) +
  geom_hline(yintercept = ac153, col = "red", lty = 2, linewidth = 0.25) +
  xlab("Length (in mm)") +
  ylab("VBGF est. Age") +
  theme_classic() +
  xlim(0, 200) +
  ylim(0,5) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(color = "black", size = 35, face = "bold.italic"))

#### END ####
