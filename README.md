The author (s) of each of the MATLAB programs below are listed within the .m file.

EPCAPE_SEMS_APS_merge_0716.m: This MATLAB program reads in EPCAPE SEMS and APS data measured at Mt. Soledad. 
The SEMS and APS data used for this study can be found in the UC San Diego library's digital collections (https://doi.org/10.6075/J0Q81D8N).
This code merges the APS and SEMS measurement data to form a 5-minute averaged merged size distribution following Khlystov et al (2004), Modini et al (2015), and Dedrick et al (2024).
Outputs are .mat files of merged-size distribution.


EPCAPE_SMPS_APS_merge_0716.m: This MATLAB program reads in EPCAPE SMPS and APS data measured at Pier.
The SMPS and APS data used for this study were obtained from DOE ARM Data Discovery. APS (http://dx.doi.org/10.5439/1407135) and 
SMPS (http://dx.doi.org/10.5439/1476898). This code merges the measurement data to form a 5-minute averaged merged size distribution following 
Khlystov et al (2004), Modini et al (2015), and Dedrick et al (2024). Outputs a .mat file of merged size distribution.

EPCAPE_FittingModel_soledad_v20240716.m: This MATLAB program fits three lognormal modes to 5-minute averaged merged measured size distribution
from APS and SMPS during EPCAPE at Mt. Soledad. Saves lognormal modal fits to a .mat file. A section dedicated to graphing  2 by 2 subplots of the merged
size distribution and their modal fit is available.

EPCAPE_FittingModel_pier_v20240716.m: This MATLAB program fits three lognormal modes to 5-minute averaged merged measured size distribution
from APS and SMPS during EPCAPE at Pier. Saves lognormal modal fits to a .mat file. A section dedicated to graphing  2 by 2 subplots of the merged
size distribution and their modal fit is available.

create_PSD.m: This MATLAB program creates number and mass size distribution based on lognormal fitting parameters (number, mean size, and width). 
Outputs particle number size distribution (PNSD) and particle mass size distribution(PMSD).

Hoppel_min_soledad_v2024_0905.m: This MATLAB program uses the fitted aerosol modes to retrieve the Hoppel Minimum diameter at Mt. Soledad. Program for retrieving
Hoppel Minimum diameter using the intersection of fitted modes and the minimum from measured size distribution is also included (Dedrick et al., 2024).
Section to retrieve the Hoppel Minimum only for the distribution where the category for distribution is bimodal is also included. Concatenates Hoppel Minimum
identified from merged size distribution to those identified from SEMS distribution when APS was unavailable. Whether to include bimodal distribution conditions or not can be customized.

Hoppel_min_pier_v2024_0905.m:  This MATLAB program uses the fitted aerosol modes to retrieve the Hoppel Minimum diameter at Pier. Program for retrieving
Hoppel Minimum diameter using the intersection of fitted modes and the minimum from measured size distribution is also included (Dedrick et al., 2024).
Section to retrieve the Hoppel Minimum only for the distribution where the category for distribution is bimodal is also included. Concatenates Hoppel Minimum
identified from merged size distribution to those identified from SMPS distribution when APS was not available. Whether to include bimodal distribution conditions or not can be customized.

EPCAPE_SMPS_mode_fitting_v2024_0905.m: This MATLAB program fits two lognormal modes of accumulation and Aitken to 5-minute averaged SMPS measured size distribution
of SMPS during EPCAPE at Pier. Saves lognormal modal fits to a .mat file. A modified version of EPCAPE_FittingModel_pier_v20240716.m to allow for the fitting of accumulation and
Aitken mode when APS was not available.

EPCAPE_SEMS_mode_fitting_v2024_0906.m: This MATLAB program fits two lognormal modes of accumulation and Aitken to 5-minute averaged SEMS measured size distribution
of SEMS during EPCAPE at Mt. Soledad. Saves lognormal modal fits to a .mat file. A modified version of EPCAPE_FittingModel_soledad_v20240716.m to allow for the fitting of accumulation and
Aitken mode when APS was not available.

Pier_category_v20240819 (1).csv and Soledad_category_v20240807 (2).csv

Dates are in UTC time and categorization of distribution and fits were made as follows:
value of 0: PNSD data was missing for the measurement period (this value seems to be pretty high for Soledad(1/3 of all period) most likely due to SEMS not being run while Soledad was in-cloud)
value of 1: Observed bimodal distribution in the submicron range (<1μm) and fits of those two modes in the submicron range were correct
value of 2: Observed unimodal distribution in the submicron range
value of 3: All other cases. Such as when PNSD distribution was bimodal but not correctly fit, more than two modes in the submicron range, etc...

Link to CSV file for size distributions and modal fit: https://drive.google.com/file/d/1_mUsjD_1pXbC4ro78WDKl0JG4Ydm3euH/view?usp=sharing
File was too big to be uploaded to GitHub


EPCAPE_SMPS_APS_merged_5min_20240722.csv for Pier and EPCAPE_SEMS_APS_merged_5min_20240722.csv for Mt. Soledad:

These merged size distributions only include the time when both APS and SEMS/SMPS time were available from measurement; Either as NaN or with value. 

Notably, Pier does not have merged distribution 9/18~10/23

Variables:
time_5min: UTC of recorded measurements; They are five minutes apart as these were 5-minute averaged distributions.
PNSD_merged_5min_1, PNSD_merged_5min_2,.....PNSD_merged_5min_100:Aerosol particle size distribution corresponding to D_merged. In units of cm^-3
D_merged (first 100 rows): Aerosol particle diameter corresponding to the 5-minute averaged size distribution data. These diameters are logarithmically spaced.  (0.01-10 µm) 
e.g. PNSD_merged_5min_1 corresponds to the 1st row diameter of D_merged and so on.
Nt_merge_5min:5 minute averaged total integrated number concentration of aerosols cm^-3 based on the particle size distribution
Rho_merge_5min: Effective density used for merging APS and SEMS/ APS and SMPS


EPCAPE_pier_mode_fit_out_0722.csv for Pier and EPCAPE_soledad_mode_fit_out_0722 for Mt. Soledad

This distribution fits include those when there were only submicron SEMS/SMPS available but not APS. They were fit under the assumption that sea spray mode was in APS, giving two modal fits in the submicron range. All of the values below were constrained by the range developed by Dr.Derdrick in previous LASIC study

Variables:
Time_5min: UTC at which either merged size distribution was used to fit the modes (or SEMS/SMPS if APS was not available)
sea_spray_N, accumulation_N,  aitken_N: Fitted number concentration of each of the modes cm^-3
sea_spray_diam, accumlation_diam, aitken_diam: Center/peak of each of the modes that is fitted based on the fitting model; In units of micrometer
sea_spray_gsd, accumulation_gsd, aitken_gsd: Residual of the fitted model for the objective function that was used; 
fit_flag: A value of 1 here indicates that only SMPS/SEMS were used to fit that particular size distribution. Only accumulation and Aitken modal parameters would be available.



