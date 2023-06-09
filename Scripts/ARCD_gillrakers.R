#### ARCD_GR-DR.R ####
## Script designed for the analyses of Arctic cod Gillrakers. Gillrakers on the
## first gill arch on the right side of 37 individual fish were extracted. The 
## length of the gill arch was recorded and the number of gillrakers on the 
## gill arch was also recorded. The gillraker density was calculated by dividing the 
## number of gillrakers by the length of the gill arch.

## All three variables were tested for adherence to normality and assessed for 
## significant differences using linear models in R. 

## written by dr and jm May 2022, updated June 2023.

#### Set up ####

## Resetting all instances
rm(list=ls())

## Loading required libraries
{
  library(ggplot2)
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
  library(mosaic)
  library(tidyr)
  library(lsr)
  library(car)
  library(emmeans)
}

## Set up working directory
setwd("/pathtofile/")

## Reading in the data file "ARCD_GR_data.csv"
gillrakers <- read.csv("ARCD_GR_data.csv", header = T, stringsAsFactors = T, na.strings = "")

## View the data 
View(gillrakers)

## Eliminate any individuals with missing values
gillrakers <- gillrakers[complete.cases(gillrakers), ]

## Set the age (VBage) variable as a factor 
gillrakers$VBage <- as.factor(gillrakers$VBage)

## Make variables available without the need for subsetting 
attach(gillrakers)


#### Histograms ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Making histograms of the raw data to look at their distribution
## Histogram function for looking at the distribution of the variables 
histoloop<-function (x,...) {
  par(mfrow=c(2,2))
  for (a in 1:ncol(x)) {
    hist(x[,a], main=colnames(x)[a], xlab=colnames(x)[a],
         xlim = , ylim = NULL,
         cex.lab=1.5, cex.main=1.7, col="steelblue4", cex.axis=1.4, lwd=2)
  }
}  

## Re-assembling the gillrakers data to make it easier to plot. 
newdat<-cbind.data.frame(length=length_mm, lenarch=len_arch, numrakers=num_rakers,
                         rakerdens=num_len_ratio)

## applying the histoloop function to see the histograms
histoloop(newdat)

#### Analysis of Gillraker numbers by habitat and age ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Resetting the order of the habitat classes to print out as outlined in
## Malizia et al. 
gillrakers$maj_class = factor(gillrakers$maj_class,
                              levels = c("NS", "OS", "US", "LS"))

## Use Levene's test to assess departures from normality in the gillraker numbers.
## Results that are not sign (< 0.05), are not considered to be deviating from 
## a normal distribution.
levgrnhab <- leveneTest(num_rakers ~ maj_class)
levgrnhab

## Run a linear model to test for significant differences between gill raker numbers 
## among habitats.
grnum_mcL <- lm(num_rakers ~ maj_class, data = gillrakers)

## Running a Shapiro-Wilk's test to assess the normality of residuals from 
## the linear model. Results that are not sign. (< 0.05), are not considered to 
## be deviating from a normal distribution.
shapiro_test(residuals(grnum_mcL))

## Quickly assess major deviations in the qqplot.
ggqqplot(residuals(grnum_mcL))

## Summarising the linear model to get the coefficients
summary(grnum_mcL)

## Testing the hypothesis of no sign. differences among habitats
anova(grnum_mcL)
etaSquared(grnum_mcL)

## Generating a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
grnum_mcL_res <- as.data.frame(emmeans(grnum_mcL, "maj_class",type = "response"))

## Applying a post-hoc Tukey's Honestly Significant test to determine which habitats 
## are significantly different form which using the emmeans function pairs.
pairs(emmeans(grnum_mcL, "maj_class"))

## The tabulated transformed coeff and 95% CI can be visualised: 
grnum_mcL_res

#### Supposed to have been done ####
## this would have been a more appropriate test for gillraker numbers.
#grnum_mc <- glm(num_rakers ~ maj_class, family = poisson(link="log"), data = gillrakers)
#summary(grnum_mc)
#anova(grnum_mc, test = "Chisq")
#grnum_mc_res <- as.data.frame(emmeans(grnum_mc, "maj_class",type = "response"))
#marg_grnum_mcres <- emmeans(grnum_mc, ~ maj_class)
#pairs(marg_grnum_mcres)
## However, given that this would not change the overall results of the analyses,
## we default to the current accepted analyses.

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the gillraker counts by habitat
ggplot(grnum_mcL_res, aes(x = maj_class, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(gillrakers, mapping = aes(x = maj_class, y = num_rakers, 
              color = maj_class), position = position_jitter(0.3), size = 5, stroke = 2.5) +
  labs(x = "", y = "Gillraker numbers", color = "maj_class") +
  scale_color_manual(values = c("#D55E00", "#CC79A7", "#6A5ACD","#999999")) +
  scale_y_continuous(limits = c(20, 50), breaks = seq(20, 50, by = 5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 35),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")


#@@@@ Do same thing as above but for the age variable @@@@#

## Use Levene's test to assess departures from normality in the gillraker numbers.
## Results that are not sign (< 0.05), are not considered to be deviating from 
## a normal distribution.
levgrnage <- leveneTest(num_rakers ~ VBage)
levgrnage

## Run a linear model to test for significant differences between gill raker numbers 
## among age classes.
grnum_aL <- lm(num_rakers ~ VBage, data = gillrakers)

## Running a Shapiro-Wilk's test to assess the normality of residuals from 
## the linear model. Results that are not sign. (< 0.05), are not considered to 
## be deviating from a normal distribution.
shapiro_test(residuals(grnum_aL))

## Quickly assess major deviations in the qqplot.
ggqqplot(residuals(grnum_aL))

## Summarising the linear model to get the coefficients
summary(grnum_aL)

## Testing the hypothesis of no sign. differences among age classes
anova(grnum_aL)
etaSquared(grnum_aL)

## Generating a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
grnum_aL_res <- as.data.frame(emmeans(grnum_aL, "VBage",type = "response"))

## Applying a post-hoc Tukey's Honestly Significant test to determine which habitats 
## are significantly different form which using the emmeans function pairs.
pairs(emmeans(grnum_aL, "VBage"))

## The tabulated transformed coeff and 95% CI can be visualised: 
grnum_aL_res

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the gillraker counts by age class
ggplot(grnum_aL_res, aes(x = VBage, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(gillrakers, mapping = aes(x = VBage, y = num_rakers, 
                                        color = VBage), position = position_jitter(0.3), size = 5, stroke = 2.5) +
  labs(x = "", y = "Gillraker numbers", color = "maj_class") +
  scale_color_manual(name = "Age class", values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_y_continuous(limits = c(20, 50), breaks = seq(20, 50, by = 5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 35),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")

#### Analysis of Gill arch length habitat and age ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Starting again with Levene's test assessing departures from normality in the gill arch lengths.
## Results that are not sign (< 0.05), are not considered to be deviating from 
## a normal distribution.
levgamc <- leveneTest(len_arch ~ maj_class)
levgamc

## Run a linear model to test for significant differences between gill arch length 
## among habitats.
gra_mcL <- lm(len_arch ~ maj_class, data = gillrakers)

## Running a Shapiro-Wilk's test to assess the normality of residuals from 
## the linear model. Results that are not sign. (< 0.05), are not considered to 
## be deviating from a normal distribution.
shapiro_test(residuals(gra_mcL))

## Quickly assess major deviations in the qqplot.
ggqqplot(residuals(gra_mcL))

## Summarising the linear model to get the coefficients
summary(gra_mcL)

## Testing the hypothesis of no sign. differences among habitats
anova(gra_mcL)
etaSquared(gra_mcL)

## Generating a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
gra_mcL_res <- as.data.frame(emmeans(gra_mcL, "maj_class",type = "response"))

## Applying a post-hoc Tukey's Honestly Significant test to determine which habitats 
## are significantly different form which using the emmeans function pairs.
pairs(emmeans(gra_mcL, "maj_class"))

## The tabulated transformed coeff and 95% CI can be visualised: 
gra_mcL_res

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the gill arch by habitat
ggplot(gra_mcL_res, aes(x = maj_class, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(gillrakers, mapping = aes(x = maj_class, y = len_arch, 
                                        color = maj_class), position = position_jitter(0.3), size = 5, stroke = 2.5) +
  labs(x = "", y = "Length of gill arch (mm)", color = "maj_class") +
  scale_color_manual(values = c("#D55E00", "#CC79A7", "#6A5ACD","#999999")) +
  scale_y_continuous(limits = c(10, 40), breaks = seq(10, 40, by = 5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 35),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")


#@@@@ Do same thing as above but for the age variable @@@@#

## Use Levene's test to assess departures from normality in the gill arch.
## Results that are not sign (< 0.05), are not considered to be deviating from 
## a normal distribution.
levgaa <- leveneTest(len_arch ~ VBage)
levgaa

## Run a linear model to test for significant differences between gill arch 
## among age classes.
gra_aL <- lm(len_arch ~ VBage, data = gillrakers)

## Running a Shapiro-Wilk's test to assess the normality of residuals from 
## the linear model. Results that are not sign. (< 0.05), are not considered to 
## be deviating from a normal distribution.
shapiro_test(residuals(gra_aL))

## Quickly assess major deviations in the qqplot.
ggqqplot(residuals(gra_aL))

## Summarising the linear model to get the coefficients
summary(gra_aL)

## Testing the hypothesis of no sign. differences among age classes
anova(gra_aL)
etaSquared(gra_aL)

## Generating a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
gra_aL_res <- as.data.frame(emmeans(gra_aL, "VBage",type = "response"))

## Applying a post-hoc Tukey's Honestly Significant test to determine which habitats 
## are significantly different form which using the emmeans function pairs.
pairs(emmeans(gra_aL, "VBage"))

## The tabulated transformed coeff and 95% CI can be visualised: 
gra_aL_res

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the gill arch by age class
ggplot(gra_aL_res, aes(x = VBage, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(gillrakers, mapping = aes(x = VBage, y = len_arch, 
                                        color = VBage), position = position_jitter(0.3), size = 5, stroke = 2.5) +
  labs(x = "", y = "Length of gill arch (mm)", color = "maj_class") +
  scale_color_manual(name = "Age class", values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_y_continuous(limits = c(10, 40), breaks = seq(10, 40, by = 5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 35),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")

#### Analysis of Gillraker density by habitat and age ####
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Starting with Levene's test assessing departures from normality in the gillraker density.
## Results that are not sign (< 0.05), are not considered to be deviating from 
## a normal distribution.
levgdmc <- leveneTest(num_len_ratio ~ maj_class)
levgdmc

## Run a linear model to test for significant differences between gillraker density 
## among habitats.
grd_mcL <- lm(num_len_ratio ~ maj_class, data = gillrakers)

## Running a Shapiro-Wilk's test to assess the normality of residuals from 
## the linear model. Results that are not sign. (< 0.05), are not considered to 
## be deviating from a normal distribution.
shapiro_test(residuals(grd_mcL))

## Quickly assess major deviations in the qqplot.
ggqqplot(residuals(grd_mcL))

## Summarising the linear model to get the coefficients
summary(grd_mcL)

## Testing the hypothesis of no sign. differences among habitats
anova(grd_mcL)
etaSquared(grd_mcL)

## Generating a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
grd_mcL_res <- as.data.frame(emmeans(grd_mcL, "maj_class",type = "response"))

## Applying a post-hoc Tukey's Honestly Significant test to determine which habitats 
## are significantly different form which using the emmeans function pairs.
pairs(emmeans(grd_mcL, "maj_class"))

## The tabulated transformed coeff and 95% CI can be visualised: 
grd_mcL_res

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the gillraker density by habitat
ggplot(grd_mcL_res, aes(x = maj_class, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(gillrakers, mapping = aes(x = maj_class, y = num_len_ratio, 
                                        color = maj_class), position = position_jitter(0.3), size = 5, stroke = 2.5) +
  labs(x = "", y = "Gillraker density", color = "maj_class") +
  scale_color_manual(values = c("#D55E00", "#CC79A7", "#6A5ACD","#999999")) +
  scale_y_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, by = 0.5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 35),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")


#@@@@ Do same thing as above but for the age variable @@@@#

## Use Levene's test to assess departures from normality in the gillraker density.
## Results that are not sign (< 0.05), are not considered to be deviating from 
## a normal distribution.
levgda <- leveneTest(num_len_ratio ~ VBage)
levgda

## Run a linear model to test for significant differences between gillraker density
## among age classes.
grd_aL <- lm(num_len_ratio ~ VBage, data = gillrakers)

## Running a Shapiro-Wilk's test to assess the normality of residuals from 
## the linear model. Results that are not sign. (< 0.05), are not considered to 
## be deviating from a normal distribution.
shapiro_test(residuals(grd_aL))

## Quickly assess major deviations in the qqplot.
ggqqplot(residuals(grd_aL))

## Summarising the linear model to get the coefficients
summary(grd_aL)

## Testing the hypothesis of no sign. differences among age classes
anova(grd_aL)
etaSquared(grd_aL)

## Generating a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
grd_aL_res <- as.data.frame(emmeans(grd_aL, "VBage",type = "response"))

## Applying a post-hoc Tukey's Honestly Significant test to determine which habitats 
## are significantly different form which using the emmeans function pairs.
pairs(emmeans(grd_aL, "VBage"))

## The tabulated transformed coeff and 95% CI can be visualised: 
grd_aL_res

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the gillraker density by age class
ggplot(grd_aL_res, aes(x = VBage, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  geom_point(shape = 15, size = 9) +
  geom_jitter(gillrakers, mapping = aes(x = VBage, y = num_len_ratio, 
                                        color = VBage), position = position_jitter(0.3), size = 5, stroke = 2.5) +
  labs(x = "", y = "Gillraker density", color = "maj_class") +
  scale_color_manual(name = "Age class", values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_y_continuous(limits = c(0.5, 3), breaks = seq(0.5, 3, by = 0.5)) +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 40),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 35),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.4,"cm"),
        legend.position = "none")

#### END ####
