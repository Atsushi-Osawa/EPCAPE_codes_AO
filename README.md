The author (s) of each of the MATLAB programs below are listed within the .m file.

EPCAPE_SEMS_APS_merge_0716.m: This MATLAB program reads in EPCAPE SEMS and APS data measured at Mt. Soledad for EPCAPE. 
The SEMS and APS data used for this study can be found in the UC San Diego library's digital collections (https://doi.org/10.6075/J0Q81D8N).
This code merges the measurement data to form a 5-minute averaged merged size distribution following Khlystov et al (2004), Modini et al (2015), and Dedrick et al (2024).
Outputs .mat files of merged size distribution.


EPCAPE_SMPS_APS_merge_0716.m: This MATLAB program reads in EPCAPE SMPS and APS data measured at Mt. Soledad for EPCAPE. 
The SMPS and APS data used for this study were obtained from DOE ARM Data Discovery. APS (http://dx.doi.org/10.5439/1407135) and 
SMPS (http://dx.doi.org/10.5439/1476898). This code merges the measurement data to form a 5-minute averaged merged size distribution following 
Khlystov et al (2004), Modini et al (2015), and Dedrick et al (2024). Outputs .mat files of merged size distribution.

EPCAPE_FittingModel_soledad_v20240716.m: This MATLAB program fits three lognormal modes to 5-minute averaged merged measured size distribution
from APS and SMPS during EPCAPE at Pier. Saves lognormal modal fits to a .mat file. A section dedicated to graphing  2 by 2 subplots of the merged
size distribution and their modal fit is available.

EPCAPE_FittingModel_pier_v20240716.m: This MATLAB program fits three lognormal modes to 5-minute averaged merged measured size distribution
from APS and SMPS during EPCAPE at Pier. Saves lognormal modal fits to a .mat file. A section dedicated to graphing  2 by 2 subplots of the merged
size distribution and their modal fit is available.

create_PSD.m: This MATLAB program creates number and mass size distribution based on lognormal fitting parameters (number, mean size, and width). 
Outputs particle number size distribution (PNSD) and particle mass size distribution(PMSD).

Hoppel_min_soledad_v20240905.m: This MATLAB program uses the fitted aerosol modes to retrieve the Hoppel Minimum diameter at Mt. Soledad. Program for retrieving
Hoppel Minimum diameter using the intersection of fitted modes and the minimum from measured size distribution is also included (Dedrick et al., 2024).
Section to retrieve the Hoppel Minimum only for the distribution where the category for distribution is bimodal is also included. Concatenates Hoppel Minimum
identified from merged size distribution to those identified from SEMS distribution when APS was not available.

Hoppel_min_pier_v2024_0905.m:  This MATLAB program uses the fitted aerosol modes to retrieve the Hoppel Minimum diameter at Pier. Program for retrieving
Hoppel Minimum diameter using the intersection of fitted modes and the minimum from measured size distribution is also included (Dedrick et al., 2024).
Section to retrieve the Hoppel Minimum only for the distribution where the category for distribution is bimodal is also included. Concatenates Hoppel Minimum
identified from merged size distribution to those identified from SMPS distribution when APS was not available.

EPCAPE_SMPS_mode_fitting_v2024_0905.m: This MATLAB program fits two lognormal modes of accumulation and Aitken to 5-minute averaged merged measured size distribution
of SMPS during EPCAPE at Pier. Saves lognormal modal fits to a .mat file. A modified version of EPCAPE_FittingModel_pier_v20240716.m to allow for the fitting of accumulation and
Aitken mode when APS were not available.

EPCAPE_SEMS_mode_fitting_v2024_0906.m: This MATLAB program fits two lognormal modes of accumulation and Aitken to 5-minute averaged merged measured size distribution
of SMPS during EPCAPE at Pier. Saves lognormal modal fits to a .mat file. A modified version of EPCAPE_FittingModel_soledad_v20240716.m to allow for the fitting of accumulation and
Aitken mode when APS were not available.

