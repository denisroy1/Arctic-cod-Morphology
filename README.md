# Arctic-cod-Morphology
Script repository for Arctic cod morphological assessments using fish from the Beaufort Sea. This repository comprises all the relevant scripts and data files for the manuscript titled: 
### Morphology of Arctic cod (Boreogadus saida) assessed according to habitat preference and age in the Beaufort Sea 

Authored by: Juliano Malizia, Marie Launay, Ingrid Marie Bruvold, María Quintela, Torild Johansen, James D. Reist, Andrew R. Majewski, and Denis Roy 

#
### Scripts:

To see the main scripts used for this research, go to 'scripts' folder. Script definitions:

* rev_VBGF.R: This script takes in the raw Arctic cod meta data from the file “Arctic_cod_raw.csv” and uses the parameters derived by Forster et al. 2020 Deep-Sea Research Part II: [doi 10.1016/j.dsr2.2020.104779](https://www.sciencedirect.com/science/article/pii/S0967064519302668) to estimate the age of individual collected fish from their lengths (in mm). It plots out the data into year classes and creates a revised dataset “ACmeta.csv” that can be used for the geometric and traditional morphometric analyses.

* ARCD_GM_trad.R: This script takes in the “ACmeta.csv” data, along with the digitised landmark data “ARC_22_rmv.TPS” generated from the tpsUTIL32 and tpsDIGw64 softwares available at: https://www.sbmorphometrics.org/ 
Using these data, the script performs a Generalised Procrustes Analysis (GPA) on the landmark data and converts it to a geometric morphometrics dataframe used in the remaining analyses. It tests for allometry and corrects the GM data for size-based relationship to shape. It will then perform a PCA on the shape variables, calculate the individual PC scores, and assess the percent contribution of each landmark to the PCs. It will the plot the PCA first by habitat aggregation and then by age class. It also performs the MANCOVA for the shape data and will generate habitat and age specific deformation figures. It will do the same analyses for the traditional morphology as well.

* ARCD_gillrakers.R: This script use the “ARCD_GR_data.csv” file which is a recapitulation of the metadata but with  the gillraker data for the selected individuals for which the data was available. Gillraker data are tested for normality and then tested for significant differences using the lm function in base R. This generates the coefficients for the model and the model is tested using the anova() function. Post-hoc Tukey’s HSD tests are run to determine which habitats and which age classes are significantly different from one another. Results are visualised using stripcharts with the mean and 95%CI superimposed.

* VBGF test curve.R: A script that is used to demonstrate how the rounding of lengths and/or VBGF generated age estimate could impact our age classes. Script demonstrates how the steepness of the VBGF curve can introduce bias in the age estimates especially with fish at larger lengths. It compares the age estimates for toy datasets generated using the shorthorn sculpin and the Arctic cod parameters as listed in Forster et al. 2020.

#
### Data:

To see the datasets used for this research go to the 'data' folder or to the dryad repository at the following link:XXXXXXXXXXXXXX

* Arctic_cod_raw.csv: Contains the raw data forming the base of the analyses. Individual fish are listed, where they were collected from, their lengths (in mm), whether a picture was taken, and if that picture was good enough for GM analyses.

* ACmeta.csv: A filtered version of the above (with the same FishID) that lists the relevant data for the main GM and TM analyses, including age (generated from the rev_VBGF.R above) and the habitat aggregations as assessed using Majewski et al. 2017; Deep Sea Res. Part I: Oceanogr. Res. Pap. 121: 169-182.

* ARC_22_rmv.TPS: Landmark configurations for the individual Arctic cod used in our analyses. Data generated from standardised fish pictures using the ‘tps’ softwares listed above.

* ARCD_GR_data.csv: The subset of individuals (indexed with the same FIshID) for which the gillraker data was generated. 
