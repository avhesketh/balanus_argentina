# Ecological and environmental context shape the differential effects of a facilitator in its native and invaded ranges

## Abstract

Invasive species are often studied for the disproportionately strong negative effects they have in their introduced range compared to their native range and the mechanisms causing these effects. Less attention has been paid to whether invasive species, particularly those that are important as facilitators in their native range, have persistent positive effects in their invaded range. In this study, we manipulated the density of the high intertidal acorn barnacle Balanus glandula in its native (Bluestone Point, British Columbia, Canada) and invaded range (Punta Ameghino, Chubut Province, Argentina) in combination with herbivore density to determine how this facilitator differentially affects associated species at these two sites. Given that high intertidal species at Punta Ameghino (PA) are evolutionarily naïve to barnacles, we expected the positive effects of B. glandula at Bluestone Point (BP) to be weakened or absent at PA. Balanus glandula had a positive association with perennial algal cover at BP, but a negative association with perennial algal cover at PA. Ephemeral algal cover at BP was negatively associated with barnacle cover, but the same association was positive at PA. Herbivore abundance, meanwhile, was positively correlated with barnacle cover in both systems, and the strength of this interaction was similar in both the native and invaded range. These results suggest that shared evolutionary history is not a prerequisite for species to benefit from a novel facilitator, and further reinforce the importance of the traits of associated species and environmental stress in governing the strength of facilitative interactions.


## Description of raw data files

Within the raw data folder, there are a number of data frames used for eventual data visualization and analysis.

1) community_PA_200602.csv: monitoring data collected at Punta Ameghino including _Balanus glandula_ recruitment, and algal cover changes through time.

2) community_PA_200512.csv: community data from Punta Ameghino for the first timepoint (when treatments were established) including adult limpet and starting algal cover.

3) siphonaria_PA.csv: _Siphonaria lessonii_ abundances within plots for timepoints after 2005-12.

4) humidity_PA.csv: relative humidity and other environmental data from Punta Ameghino collected January 2007 to May 2020 at a CENPAT meteorological station

5) temperature_PA_200603.csv: intertidal temperature data collected at Punta Ameghino from December 2005-March 2006

6) community_BP_200707.csv: monitoring data collected at Bluestone Point at all timepoints

7) temperature_BP_200609.csv: intertidal temperature data collected at Bluestone Point from June 2006-September 2006

8) humidity_BP.csv: relative humidity collected in Tofino, British Columbia, at a DFO weather station with very similar conditions to the Bluestone Point study site from January 2015 to September 2020.

9) fucus_BP_200707.csv: Cover of _Fucus distichus_ perennial algae at Bluestone Point over the course of the experiment extracted from either data notebooks or images using ImageJ.

10) siphonaria_dwsl_Tablado_SciMar_2001.csv: dry weight and shell length data for _Siphonaria lessonii_ at different shores in Argentina extracted from Tablado and Gappa, 2001.

11) grazer_size_BP_PA.csv: the size (shell length for _Lottia digitalis_ and _Siphonaria lessonii_ and shell height for _Littorina scutulata_) of grazers estimated from photographs taken at each site at three different timepoints for use in approximating biomass in plots from abundance data.

12) littorina_dwsl_North_1954.png: figure of the relationship between shell height and tissue weight of _Littorina scutulata_ ,taken from North 1954 for extracting data with metaDigitise to establish a relationship between log(shell height) and log(tissue weight).

13) lottia_dwsl_Frank_1965.png: figure of the relationship between shell length and volume (which is easy related to dry weight) for _Lottia digitalis_, taken from Frank 1965 for extracting data with metaDigitise to establish a relationship between log(shell length) and log(tissue weight).

14) littorina_lottia_dwsl_North_Frank.csv: extracted data from figures 12 and 13.

Note that the 'caldat' folder contains objects created with the metaDigitise package that hold extracted data from the images (items 12 and 13). These data have been converted to a csv file (#14).

## Description of scripts and outputs

1) 1_biotic_tidying.R: tidying of raw data and conversion into a useable format (clean data outputs = bio_responses.csv, fucus_clean.csv)
2) 2_abiotic.r: analysis and plotting of temperature and humidity data (clean data outputs = humidity_clean.csv, temperature_clean.csv, Figures 1a-1c)
3) 3_herbivores.R: conversion of herbivore abundance to biomass estimates, analysis and plotting of these data (clean data outputs = herbivore_abundance.csv, herbivore_biomass.csv, Figures 2a-2c, components of S3, S4)
4) 4_algae.R: analysis and plotting of ephemeral and perennial algal cover (Figures 3-4, S5, S6)
5) 5_barnacles.R: analysis and plotting of barnacle recruitment data (Figures 5, S7)
6) 6_maps.R: code for generating maps (Figure S1)

## Description of variables

Metadata, including a description of variables and their units of measurement, can be found under the 'metadata' folder in the file 'metadata_variables.csv'. Further descriptions of factors and their levels can be foundin 'metadata_factors.csv'. A description of files can be found (in addition to within this README) in 'metadata_files.csv'.
