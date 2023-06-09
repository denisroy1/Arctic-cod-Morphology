## ARCD_GM_trad.R
## Script that takes the data collected from Arctic cod standardised pictures, and performs the 
## Geometric morphometric analyses and the traditional morphometric analyses. Both the GM and the 
## traditional data are taken from the standardised photos.

## The fish pictures were generate by the DFO scientist on the CBSMEA 2018-2019
## cruise to the Beaufort Sea. Photos were digitised using the tpsUtil32, and 
## 22 landmarks placed on each using tpsDIGw64. Landmarks were placed as 
## described in Malizia et al. Generated landmark configuration for each
## specimen were saved in a .TPS 

## In this script we use the TPS data (ARC_22_rmv.TPS) and it's assocaited metadata
## (ACmeta.csv)
## Written by dr and jm May 29 2023

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#### Set up ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## clearing and resetting all instances (data and variables)
rm(list = ls())

## Loading needed libraries
{
  library(geomorph)
  library(ggplot2)
  library(ggfortify)
  library(ellipse)
  library(car)
  library(dplyr)
  library(effectsize)
  library(visreg)
  library(readr)
  library(vegan)
}

## Setting working directory (path to where the data files located)
setwd("/pathtofile/")

## Reading in the metadata
ARCD_GM_data <- read.csv("ACmeta.csv", header = T, na.string = "", 
                         stringsAsFactors = T, row.names = NULL)

## Making sure the VBage estimate is considered a factor and not a numeric value
ARCD_GM_data$VBage <- as.factor(ARCD_GM_data$VBage) 

#Import TPS file
ARCD_22 = readland.tps("ARC_22_rmv.tps", specID = "ID")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#### Geometric Morphometrics ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Perform Generalised Procrustes Analysis using gpagen from geomorph (Adams et al. 2021)
ARCD_22_GPA = gpagen(ARCD_22,curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                   max.iter = NULL, ProcD = TRUE, Proj = TRUE,
                   print.progress = TRUE)

## Use the Consensus configuration as reference for shape deformations. This is 
## only important much later when making the habitat and age specific deformation
## plots
ref <- ARCD_22_GPA$consensus

## Create a dataframe with geometric morphometric data and metadata
## Here, we include age (VBage) and habitat class (maj_class)
gdf = geomorph.data.frame(ARCD_22_GPA,
                          VBage = ARCD_GM_data$VBage, 
                          maj_class = ARCD_GM_data$maj_class,
                          flength = ARCD_GM_data$length_mm)

## Need to check the data to see if there is an allometric relationship that 
## should be accounted for. To do this we follow the details outlined in 
## Adams et al. 2021, in the vignette associated with allometric analyses.
## Performing size correction of shape data using the procD.lm cmd.
## This will create an object we can then use to test for a significant relationship.
sc_gdf<-procD.lm(coords~log(flength),iter=999, data=gdf)
sc_gdf

## Testing whether the allometric relationship is significant
anova(sc_gdf)

## The relationship is significant and there is likely a shape size based 
## shape effect that should be considered. We can use the residuals of the 
## sc_gdf to account for possible shape difference based on size.

## However, it would be important to check whether the allometric relationship 
## is different for the different habitat types and/or VB age groups. For this 
## we plot the allometry and look for likely large deviations in the slopes 
## of the differnt groups. So, first by habitat:
plotAllometry(sc_gdf, size = gdf$Csize, logsz = TRUE, method = "PredLine",
              pch = 19, col = as.numeric(gdf$maj_class))

## Second by age: 
plotAllometry(sc_gdf, size = gdf$Csize, logsz = TRUE, method = "PredLine",
              pch = 19, col = as.numeric(gdf$VBage))

## Plots of the allometric relationship demonstrates the best fit line would likely apply 
## to all habitat and age classes equally well and so one common allotmetric 
## model is likely to work to size correct the shape data. So, we take the residuals 
## from the procD.lm with size and consider these as size corrected shape variables in 
## the GM remaining analyses.

## Perform a PCA analysis on the size corrected shape variables, using the default 
## parameters fo gm.
gmpc <- gm.prcomp(sc_gdf$residuals,phy=NULL,align.to.phy=NULL,GLS=FALSE,
                       transform=FALSE)

## Compute total variance explained by each PC using the sdev of each 
## as calculated in the gm.prcomp.
gmpcvar <- gmpc$sdev^2 / sum(gmpc$sdev^2) * 100

## Generate a barplot <- screeplot to demonstrate the % variance explained 
## by each shape PC.
barplot(gmpcvar, ylab = "% Variance explained", cex.names = 0.8,
        col = "steelblue4", ylim = c(0,30), las = 1)

## Construct a dataframe of the pc scores and add the major habitat 
## and age classes for each individual. Also going to attach log length 
## in mm.
gmpc_data <- as.data.frame(gmpc$x)
gmpc_data$maj_class <- gdf$maj_class
gmpc_data$VBage <- gdf$VBage
gmpc_data$loglen <- log(gdf$flength)

## This is a helper function that can transform the raw Eigenvectors to % contributions to of each landmark 
## to each PC using the sdev of each PC in the data.
scalecont <- function(rotation, sdev) {rotation*sdev}

## Assigning parts of the gmpc (PCA object) profile to objects for easier manipulation
gmrot <- gmpc$rotation
gmsd <- gmpc$sdev

## Applying the helper function above transforming each rotation by the sdev of that PC.
gmscont <- t(apply(gmrot, 1, scalecont, gmsd))
head(gmscont)

# Squaring the gmcont1 object to get positive values
gmscont2 <- gmscont^2
head(gmscont2)

## Sum the contributions for all landmarks
gmscontT <- apply(gmscont2, 2, sum)

## Creating a function that converts the total contribution to a percentage  
claccont = function(x, y){x*100/y}

## Using the apply function to apply the function across all PCs
gmlmcont = t(apply(gmscont2, 1, claccont, gmscontT))
print(gmlmcont)

## Option to save the gmpc_data in the working  directory
write_csv(gmpc_data,"GMPCAscores.csv")

#### Plotting the GM PCA By Habitat ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Plotting the GM PCA data along PC1 and PC2 and using habitat as a factor by which to 
## colour the data using ggplot2.
ggplot(data = gmpc_data,aes(x = Comp1, y = Comp2, color = maj_class,
                           shape = maj_class))+
  geom_hline(yintercept = 0, linewidth = 1.4)+
  geom_vline(xintercept = 0, linewidth = 1.4)+
  geom_point(mapping = aes(Comp1,Comp2), size = 6, alpha = 0.8, stroke = 2.5)+
  stat_ellipse(aes(color = maj_class), linewidth = 2, level = 0.95)+
  ## Note these next percent variance explained were determined above using cmd summary gmpc
  labs(x ="PC1: 24.9%", y = "PC2: 11.6%", color = "Habitat aggregation")+
  scale_color_manual(name = "", labels = c("LS", "NS", "OS", "US"), values = c("#999999", "#D55E00","#CC79A7","#6A5ACD")) +
  scale_shape_manual(name = "", labels = c("LS", "NS", "OS", "US"), values = c(15, 19, 17, 18)) +
  scale_x_continuous(limits = c(-0.155, 0.155),
                     breaks = c(-0.10, 0, 0.10)) +
  scale_y_continuous(limits = c(-0.09, 0.07),
                     breaks = c(-0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075)) +
  theme_classic() +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(color = "grey"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))

#### Plotting the GM PCA By Age ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Plotting the GM PCA data along PC1 and PC2 and using age as a factor by which 
## to colour the data using ggplot2.
ggplot(data = gmpc_data,aes(x = Comp1, y = Comp2, color = VBage,
                            shape = VBage))+
  geom_hline(yintercept = 0, linewidth = 1.4)+
  geom_vline(xintercept = 0, linewidth = 1.4)+
  geom_point(mapping = aes(Comp1,Comp2), size = 6, alpha = 0.8, stroke = 2.5)+
  stat_ellipse(aes(color = VBage), linewidth = 2, level = 0.95)+
  ## Note these next percent variance explained were determined above using cmd summary gmpc
  labs(x ="PC1: 24.9%", y = "PC2: 11.6%", color = "Age classes")+
  scale_color_manual(name = "", labels = c("1", "2", "3", "4"), values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  scale_shape_manual(name = "", labels = c("1", "2", "3", "4"), values = c(19, 17, 18, 15)) +
  scale_x_continuous(limits = c(-0.155, 0.155),
                     breaks = c(-0.10, 0, 0.10)) +
  scale_y_continuous(limits = c(-0.09, 0.07),
                     breaks = c(-0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075)) +
  theme_classic() +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(color = "grey"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))

#### Shape based MANCOVA ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## MANOVA for shape data to see if there is a significant difference in the size corrected 
## shape variables fo the sampled Arctic cod.
shape_mancova <- manova(cbind(Comp1, Comp2, Comp3, Comp4) ~ maj_class * VBage, data = gmpc_data)

## Summarizing the results in a table
summary(shape_mancova)

## Calculating the partial variance explained by the variables using eta_square
eta_squared(shape_mancova)

#### Habitat and age class specific shapes ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## the next set of cmds create mean shape outline (consensus configuration) and deviation of 
## from this for each habitat aggregation and age class

## First create a cc and link the appropriate landmarks as outlined in Malizia et al.
arcd_links <- define.links(ref, ptsize = 2)
#ref <- mshape(ARCD_22_GPA$coords)

maj_ns <- mshape(ARCD_22_GPA$coords[,,which(gdf$maj_class=="NS")])
plotRefToTarget(ref, maj_ns, method = "points", mag = 2, links = arcd_links, 
                gridPars=gridPar(tar.pt.bg = "#D55E00", tar.pt.size = 2,
                                 tar.link.col = "#D55E00", tar.link.lwd = 2))

maj_os <- mshape(ARCD_22_GPA$coords[,,which(gdf$maj_class=="OS")])
plotRefToTarget(ref, maj_os, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "#CC79A7", tar.pt.size = 2,
                                 tar.link.col = "#CC79A7", tar.link.lwd = 2))

maj_us = mshape(ARCD_22_GPA$coords[,,which(gdf$maj_class=="US")])
plotRefToTarget(ref, maj_us, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "#6A5ACD", tar.pt.size = 2,
                                 tar.link.col = "#6A5ACD", tar.link.lwd = 2))

maj_ls = mshape(ARCD_22_GPA$coords[,,which(gdf$maj_class=="LS")])
plotRefToTarget(ref, maj_ls, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "lightblue", tar.pt.size = 2,
                                 tar.link.col = "lightblue", tar.link.lwd = 2))

## Doing the same for age VBGF derived age classes 
age1 <- mshape(ARCD_22_GPA$coords[,,which(gdf$VBage=="1")])
plotRefToTarget(ref, age1, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "#E69F00", tar.pt.size = 2,
                                 tar.link.col = "#E69F00", tar.link.lwd = 2))

age2 <- mshape(ARCD_22_GPA$coords[,,which(gdf$VBage=="2")])
plotRefToTarget(ref, age2, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "#56B4E9", tar.pt.size = 2,
                                 tar.link.col = "#56B4E9", tar.link.lwd = 2))

age3 <- mshape(ARCD_22_GPA$coords[,,which(gdf$VBage=="3")])
plotRefToTarget(ref, age3, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "#009E73", tar.pt.size = 2,
                                 tar.link.col = "#009E73", tar.link.lwd = 2))

age4 <- mshape(ARCD_22_GPA$coords[,,which(gdf$VBage=="4")])
plotRefToTarget(ref, age4, method = "points", mag = 2, links = arcd_links,
                gridPars=gridPar(tar.pt.bg = "#F0E442", tar.pt.size = 2,
                                 tar.link.col = "#F0E442", tar.link.lwd = 2))


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#### Linear Morphometrics ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Input traditional linear measurements using their corresponding landmarks
trad_lmks <- data.frame(CP = c(16, 18), DF1 = c(10, 11), DF2 = c(12, 13),
                       DF3 = c(14, 15), VF1 = c(19, 20), VF2 = c(21, 22),
                       ML = c(1, 2), HD = c(4, 7), PF = c(8, 9), 
                       ABD = c(10, 22), PBD = c(17, 22), 
                       row.names = c("start", "end"))

## Calculate interlandmark distance (i.e. length of linear measurements)
## and log transform the data
trad_dist <- as.data.frame(interlmkdist(ARCD_22_GPA$coords, trad_lmks))
trad_dist <- log(trad_dist)

## Add age (VBage) and habitat (maj_class) to the trad_dist data
trad_dist$VBage <- ARCD_GM_data$VBage
trad_dist$maj_class <- ARCD_GM_data$maj_class

## Also add the log transformed length (in mm).
trad_dist$fork_length <- log(ARCD_GM_data$length_mm)

## As with the gm data, here we will size correct the trad_dist data to account for possible 
## allometric relationships

## Create an empty vector for the linear model (lm) coefficients used to size 
## correct linear traits.
lm_coeff <- vector()
  
## Use a loop to size correct the linear morphological traits using the formulation
## as provided by Skoglund et al. 2015. *** Note - we know there are 11 traits to be 
## size corrected.
for (i in 1:11) { 
  ## generate the lm for the ith trait
  modi <- lm(trad_dist[,i] ~ fork_length, data = trad_dist)
  
  ## get the slope coefficient for the ith trait
  lm_coeff[i] <- coefficients(modi)[2]
  
  ## use the slope to size correct the ith trait
  trad_dist[,i] <- trad_dist[,i] + lm_coeff[i]*(log10(mean(trad_dist$fork_length)) - trad_dist$fork_length)
}

## option to remove outliers
## trad_dist <- trad_dist[!(row.names(trad_dist) %in% c("60095","60109","46610")),]

## Perform PCA for the size corrected linear morphometric traits 
tmpc = prcomp(trad_dist[c(1:11)], scale = TRUE)

## As with the GM data - compute total variance explained by each PC using the sdev of each.
tmpcvar <- tmpc$sdev^2 / sum(tmpc$sdev^2) * 100

## Generate a barplot equivalent to a screeplot to demonstrate the % variance explained 
## by each shape PC.
barplot(tmpcvar, ylab = "% Variance explained", cex.names = 0.8, xlab = "Principal components",
        names.arg = colnames(tmpc$rotation), col = "steelblue4", ylim = c(0,30), las = 1)

## Create a dataframe with the PC scores 
tmpc_data = as.data.frame(tmpc$x)

## Add the maj_class and VBage to the tmpc_data dataframe
tmpc_data$maj_class <- as.factor(trad_dist$maj_class)
tmpc_data$VBage <- as.factor(trad_dist$VBage)

## Obtain the contributions of each linear trait to the different PCs
## Getting the rotations and sdev from the PCA object (here tmpc)
tmrot <- tmpc$rotation
tmsd <- tmpc$sdev

## Use the helper function "scalecont" above to transform the raw contributions to 
## those scaled by sdev in the data. Using the apply function to multiply rotation by 
## the sdev of each PC.
tmscont <- t(apply(tmrot, 1, scalecont, tmsd))
head(tmscont)

## Square the tmscont object to get positive values
tmscont2 <- tmscont^2
head(tmscont2)

## Sum the contributions for all landmarks
tmscontT <- apply(tmscont2, 2, sum)

## And find the percent contribution of each trait to each PC 
## claccont <- function(x, y){x*100/y}
tmlmcont  <- t(apply(tmscont2, 1, claccont, tmscontT))
print(tmlmcont)

## Option to save the tmpc_data in the working  directory
write_csv(tmpc_data,"TMPCAscores.csv")

#### Plotting the Linear PCA By Habitat ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Plotting the linear trait PCA data along PC1 and PC2 and using habitat as a factor by 
## which to colour the data using ggplot2.
ggplot(data = tmpc_data,aes(x = PC1, y = PC2, color = maj_class,
                            shape = maj_class))+
  geom_hline(yintercept = 0, linewidth = 1.4)+
  geom_vline(xintercept = 0, linewidth = 1.4)+
  geom_point(mapping = aes(PC1,PC2), size = 6, alpha = 0.8, stroke = 2.5)+
  stat_ellipse(aes(color = maj_class), linewidth = 2, level = 0.95)+
  ## Note these next percent variance explained were determined above using cmd summary gmpc
  labs(x ="PC1: 17.9%", y = "PC2: 16.8%", color = "Habitat aggregation")+
  scale_color_manual(name = "", labels = c("LS", "NS", "OS", "US"), values = c("#999999", "#D55E00","#CC79A7","#6A5ACD")) +
  scale_shape_manual(name = "", labels = c("LS", "NS", "OS", "US"), values = c(15, 19, 17, 18)) +
  scale_x_continuous(limits = c(-6, 6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(-6, 6),
                     breaks = c(-6, -4,-2, 0, 2, 4 ,6)) +
  theme_classic() +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(color = "grey"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))

#### Plotting the Linear PCA By age ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Plotting the linear trait PCA data along PC1 and PC2 and using age as a factor by 
## which to colour the data using ggplot2.
ggplot(data = tmpc_data,aes(x = PC1, y = PC2, color = VBage,
                            shape = VBage))+
  geom_hline(yintercept = 0, linewidth = 1.4)+
  geom_vline(xintercept = 0, linewidth = 1.4)+
  geom_point(mapping = aes(PC1,PC2), size = 6, alpha = 0.8, stroke = 2.5)+
  stat_ellipse(aes(color = VBage), linewidth = 2, level = 0.95)+
  ## Note these next percent variance explained were determined above using cmd summary gmpc
  labs(x ="PC1: 17.9%", y = "PC2: 16.8%", color = "Age classes")+
  scale_color_manual(name = "", labels = c("1", "2", "3", "4"), values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  scale_shape_manual(name = "", labels = c("1", "2", "3", "4"), values = c(19, 17, 18, 15)) +
  scale_x_continuous(limits = c(-6, 6),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(-6, 6),
                     breaks = c(-6, -4,-2, 0, 2, 4 ,6)) +
  theme_classic() +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(color = "grey"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))


#### Linear trait based MANCOVA ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## MANOVA for shape data to see if there is a significant difference in the size corrected 
## shape variables fo the sampled Arctic cod.
trad_mancova = manova(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11) ~ maj_class * VBage, data = tmpc_data)

## Summarizing the results in a table
summary(trad_mancova)

## Calculating the partial variance explained by the variables using eta_square
eta_squared(trad_mancova)

#### END ####
