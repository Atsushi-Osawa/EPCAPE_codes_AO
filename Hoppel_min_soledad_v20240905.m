%% hoppel_min_soledad_v20240905.m
%
% PURPOSE: Uses the fitted aerosol modes to retrieve the Hoppel Minimum
% diameter. Also retrieves Hoppel minimum diameter using intersection of
% fitted modes and the minimum of the measured size distribution. 
%
% INPUTS: Modal fit from the fitting model
% Requires the function create_PSD to be in the directory.
%
% OUTPUTS: Saves lognormal modal fits to a .mat and .txt file.
%
% AUTHOR: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         July 16, 2023
%
%         edited August 8, 2023: added comments. 
%
%
%AUTHOR: Atsushi Osawa
%        SIO
%        August 7th, 2024
%        Edited so it is compatible for Soledad and the code can be run
%        targettting the bimodal distributions.
%
%

clear all; close all; clc



%% load data
%bimodal category data of the 3 modes avaialble
load ('soledad_086.mat')
%% 

% measured size distributions
load('EPCAPE_SEMS_APS_merged_5min_20240721.mat')

% fit parameters of size distributons
load('no_CPCfilter_mode_fit_out_all_n5_SEMS_Fit_0.70_0716.mat')
%% 

% define time as the EPCAPE 5 min time
time = time_5min;

% measured size distribution (For us, PNSD of merged distirbution here)
PNSD_meas = PNSD_merged_5min;

% measured size distribution diameter (for merged distribution)
D_meas = D_merged;



%% Make Aerosol Size Distributions

% create a diameter grid over 0.1 - 10 um
D = logspace(log10(0.0104), log10(11), 250)'; 

% compute dlogDp
for i = 1:length(D)-1

    dlogDp_r(i) = log10(D(i+1)) - log10(D(i));

end
dlogDp = mean(dlogDp_r);


% Number and Mass Size distributions (reads modal fits)

[aitken_PNSD, aitken_PMSD] = create_PSD(D, ...
                        mode3_N, ...
                        mode3_diam, ...
                        mode3_gsd);
                                    
[accum_PNSD, accum_PMSD] = create_PSD(D, ...
                        mode2_N, ...
                        mode2_diam, ...
                        mode2_gsd);

[sea_spray_PNSD, sea_spray_PMSD] = create_PSD(D, ...
                        mode1_N, ...
                        mode1_diam, ...
                        mode1_gsd);

                                 


%% Hoppel Minimum Retrieval

for i = 1:length(time_5min)
%%%%%%%%%%%%%%%%%%%%%%% Get Hoppel Minimum:%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check for missing mode peaks
   check_mode_peaks = sum(isnan([mode3_diam(i), mode2_diam(i)]));
   % if no available mode data
   if check_mode_peaks > 0
       hoppel_min_idx(i)           = NaN;
       hoppel_min_meas_idx(i)      = NaN;
       hoppel_min_intersect_idx(i) = NaN;
       hoppel_min(i)               = NaN;
       hoppel_min_meas(i)          = NaN;
       hoppel_min_intersect(i)     = NaN;
       dHM_ait(i)                  = NaN;
       dHM_acc(i)                  = NaN;
       dHM_ait_all(i)              = NaN;
       num_hoppel_min(i)                  = NaN;
       num_aitken_hoppel_min(i)           = NaN;
       num_accum_hoppel_min(i)            = NaN;
      
   
      
   else
       % summed PNSD
       PNSD_sum(:,i) = nansum([aitken_PNSD(:,i), ...
                          accum_PNSD(:,i), ...
                          sea_spray_PNSD(:,i)],2);
       % summed NSD
       NSD_sum(:,i) = PNSD_sum(:,i) .* dlogDp;
       % find minimum between modes mean diameters (fitted)
       low_d_idx = find(D >= mode3_diam(i), 1, 'first'); % aitken mode peak
       hi_d_idx  = find(D <= mode2_diam(i), 1, 'last');  % accumulation mode peak   
       D_look_idx = low_d_idx:hi_d_idx; % indices to look over find indicies between max of mode 2 and 3
       % find minimum between modes mean diameters (measured)
       low_d_meas_idx  = find(D_meas >= mode3_diam(i), 1, 'first'); % aitken mode peak
       hi_d_meas_idx   = find(D_meas <= mode2_diam(i), 1, 'last');  % accumulation mode peak   
       D_look_meas_idx = low_d_meas_idx:hi_d_meas_idx; % indices to look over
       if isempty(D_look_meas_idx)
           hoppel_min_idx(i)           = NaN;
           hoppel_min_meas_idx(i)      = NaN;
           hoppel_min_intersect_idx(i) = NaN;
           hoppel_min(i)               = NaN;
           hoppel_min_meas(i)          = NaN;
           hoppel_min_intersect(i)     = NaN;
           dHM_ait(i)                  = NaN;
           dHM_acc(i)                  = NaN;
           dHM_ait_all(i)              = NaN;
  
           num_hoppel_min(i)                  = NaN;
           num_aitken_hoppel_min(i)           = NaN;
           num_accum_hoppel_min(i)            = NaN;
       else
      
       % Hoppel minimum diameter index from fitted size distribution
       hoppel_min_idx(i) = find(PNSD_sum(:,i) ==  ...
                             nanmin(PNSD_sum(D_look_idx,i)), 1, 'first');
      
       % Hoppel minimum diameter index from measured size distribution
       % (smoothed distribution)
       n_avg                  = 1; % number of points to smooth over
       smooth_PNSD_meas       = movmean(PNSD_meas(D_look_meas_idx,i),n_avg);
       hoppel_min_meas_idx_find = D_look_meas_idx(find(smooth_PNSD_meas ==  ...
                                     nanmin(smooth_PNSD_meas),1,'first'));
       if length(hoppel_min_meas_idx_find) > 1
           hoppel_min_meas_idx_find(i) = NaN;
       else
           hoppel_min_meas_idx(i) = hoppel_min_meas_idx_find;
       end
       % Hoppel minimum diameter index from intersection of curves
       % different approach to find using fitted curve interesection point
       % between mode 2 and mode 3
       %hoppel_min_intersect_idx(i) = D_look_idx(find(aitken_PNSD(D_look_idx,i) <= ...
                                                  %accum_PNSD(D_look_idx,i), 1, 'first') - 0);
       % Hoppel minimum diameter (at index value)
       hoppel_min(i) = D(hoppel_min_idx(i));
       % Hoppel minimum diameter (at index value of measured size
       % distribution)
       hoppel_min_meas(i) = D_meas(hoppel_min_meas_idx(i));
       % Hoppel minimum diameter (ait index of intersected fitted modes)
       % hoppel_min_intersect(i) = D(hoppel_min_intersect_idx(i));
       % difference between Hoppel Minimum diameter and mean Aitken
       % diameter
       dHM_ait_all(i) = abs(mode3_diam(i) - hoppel_min(i));
       % difference between Hoppel Minimum diameter and mean Accumulation
       % diameter
       dHM_acc(i) = abs(mode2_diam(i) - hoppel_min(i));
       % if less than 10 nm difference, not a real HM (morphs into acc
       % mode at LASIC)
       if dHM_ait_all(i) < 0.00
  
       hoppel_min_idx(i)           = NaN;
       hoppel_min_meas_idx(i)      = NaN;
       hoppel_min_intersect_idx(i) = NaN;
       hoppel_min(i)               = NaN;
       hoppel_min_meas(i)          = NaN;
       hoppel_min_intersect(i)     = NaN;
       dHM_ait(i)                  = NaN;
       num_hoppel_min(i)                  = NaN;
       num_aitken_hoppel_min(i)           = NaN;
       num_accum_hoppel_min(i)            = NaN;
       else
       dHM_ait(i) = dHM_ait_all(i);
   
       % integrated size distribution from Hoppel minimum diameter
       % (effectively CCN #)
       num_hoppel_min(i)    = nansum(PNSD_sum(hoppel_min_idx(i):end,i) .* dlogDp);
       % integrated modal distributions from Hoppel miniumum diameter
       num_aitken_hoppel_min(i)    = sum(aitken_PNSD(hoppel_min_idx(i):end, i) .* dlogDp);
       num_accum_hoppel_min(i)     = sum(accum_PNSD(hoppel_min_idx(i):end, i) .* dlogDp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end
       end % if not a real Hoppel Minimum
   end % mode peak check for Hoppel Minimum
end % time loop
%% 

for i=1:length(category)
    if category(i)==1
        hoppel_min(i)=hoppel_min(i);
    else
        hoppel_min(i)=NaN;
    end
end

%% 
hoeppL_min_merged=hoppel_min;
% save('Soledad_merged_hoppel_min.mat','hoppel_min','time_merged')
%% Hoppel minimum for SEMS only portion of only 2 modes
close all; clear all; clc
%% load data
TT=readmatrix('Soledad_category_v20240807.csv')
category_SEMS_only=TT(:,2);
%% 

% measured size distributions
load('SEMS_PNSD_0807.mat')

% fit parameters of size distributons
load('EPCAPE_soledad_mode_fit_out_0722.mat')
% define time as the LASIC 2-hour time
time = time_5min_final;

% measured size distribution (submicron, SMPS)
PNSD_meas = SEMS_PNSD;

% measured size distribution diameter (SMPS)
D_meas = SEMS_d_mid;


%% Make Aerosol Size Distributions

% create a diameter grid over 0.1 - 10 um
D = logspace(log10(0.0104), log10(11), 250)'; 

% compute dlogDp
for i = 1:length(D)-1

    dlogDp_r(i) = log10(D(i+1)) - log10(D(i));

end
dlogDp = mean(dlogDp_r);


% Number and Mass Size distributions (reads modal fits)
[aitken_PNSD, aitken_PMSD] = create_PSD(D, ...
                        mode3_N_final, ...
                        mode3_diam_final, ...
                        mode3_gsd_final);
                                    
[accum_PNSD, accum_PMSD] = create_PSD(D, ...
                        mode2_N_final, ...
                        mode2_diam_final, ...
                        mode2_gsd_final);

% [sea_spray_PNSD, sea_spray_PMSD] = create_PSD(D, ...
%                         mode1_N, ...
%                         mode1_diam, ...
%                         mode1_gsd);

                               
%% Hoppel Minimum Retrieval
hoppel_min=nan(length(time_5min_final),1);
for i = 1:length(time_5min_final)
   if isnan(mode1_diam_final(i)) & ~isnan(mode2_diam_final(i))
%%%%%%%%%%%%%%%%%%%%%%% Get Hoppel Minimum: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check for missing mode peaks
   check_mode_peaks = sum(isnan([mode3_diam_final(i), mode2_diam_final(i)]));
   % if no available mode data
   if check_mode_peaks > 0
       hoppel_min_idx(i)           = NaN;
       hoppel_min_meas_idx(i)      = NaN;
       hoppel_min_intersect_idx(i) = NaN;
       hoppel_min(i)               = NaN;
       hoppel_min_meas(i)          = NaN;
       hoppel_min_intersect(i)     = NaN;
       dHM_ait(i)                  = NaN;
       dHM_acc(i)                  = NaN;
       dHM_ait_all(i)              = NaN;
       num_hoppel_min(i)                  = NaN;
       num_aitken_hoppel_min(i)           = NaN;
       num_accum_hoppel_min(i)            = NaN;
      
   
      
   else
       % summed PNSD
       PNSD_sum(:,i) = nansum([aitken_PNSD(:,i), ...
                          accum_PNSD(:,i)],2); %sea spray part omitted
       % summed NSD
       NSD_sum(:,i) = PNSD_sum(:,i) .* dlogDp;
       % find minimum between modes mean diameters (fitted)
       low_d_idx = find(D >= mode3_diam_final(i), 1, 'first'); % aitken mode peak
       hi_d_idx  = find(D <= mode2_diam_final(i), 1, 'last');  % accumulation mode peak   
       D_look_idx = low_d_idx:hi_d_idx; % indices to look over find indicies between max of mode 2 and 3
       % find minimum between modes mean diameters (measured)
       low_d_meas_idx  = find(D_meas >= mode3_diam_final(i), 1, 'first'); % aitken mode peak
       hi_d_meas_idx   = find(D_meas <= mode2_diam_final(i), 1, 'last');  % accumulation mode peak   
       D_look_meas_idx = low_d_meas_idx:hi_d_meas_idx; % indices to look over
       if isempty(D_look_meas_idx)
           hoppel_min_idx(i)           = NaN;
           hoppel_min_meas_idx(i)      = NaN;
           hoppel_min_intersect_idx(i) = NaN;
           hoppel_min(i)               = NaN;
           hoppel_min_meas(i)          = NaN;
           hoppel_min_intersect(i)     = NaN;
           dHM_ait(i)                  = NaN;
           dHM_acc(i)                  = NaN;
           dHM_ait_all(i)              = NaN;
  
           num_hoppel_min(i)                  = NaN;
           num_aitken_hoppel_min(i)           = NaN;
           num_accum_hoppel_min(i)            = NaN;
       else
      
       % Hoppel minimum diameter index from fitted size distribution
       hoppel_min_idx(i) = find(PNSD_sum(:,i) ==  ...
                             nanmin(PNSD_sum(D_look_idx,i)), 1, 'first');
      
       % Hoppel minimum diameter index from measured size distribution
       % (smoothed distribution)
       n_avg                  = 5; % number of points to smooth over
       smooth_PNSD_meas       = movmean(PNSD_meas(D_look_meas_idx,i),n_avg);
       hoppel_min_meas_idx_find = D_look_meas_idx(find(smooth_PNSD_meas ==  ...
                                     nanmin(smooth_PNSD_meas),1,'first'));
       if length(hoppel_min_meas_idx_find) > 1
           hoppel_min_meas_idx_find(i) = NaN;
       else
           hoppel_min_meas_idx(i) = hoppel_min_meas_idx_find;
       end
       
       % ignore this 2 line for intersect hoppel_min_intersect_idx(i) = D_look_idx(find(aitken_PNSD(D_look_idx,i) <= ...
                                                  % accum_PNSD(D_look_idx,i), 1, 'first') - 0);
       % Hoppel minimum diameter (at index value)
       hoppel_min(i) = D(hoppel_min_idx(i));
       % Hoppel minimum diameter (at index value of measured size
       % distribution)
       hoppel_min_meas(i) = D_meas(hoppel_min_meas_idx(i));
       % Hoppel minimum diameter (ait index of intersected fitted modes)
       % hoppel_min_intersect(i) = D(hoppel_min_intersect_idx(i));
       % difference between Hoppel Minimum diameter and mean Aitken
       % diameter
       dHM_ait_all(i) = abs(mode3_diam_final(i) - hoppel_min(i));
       % difference between Hoppel Minimum diameter and mean Accumulation
       % diameter
       dHM_acc(i) = abs(mode2_diam_final(i) - hoppel_min(i));
       % if less than 10 nm difference, not a real HM (morphs into acc
       % mode at LASIC)
       if dHM_ait_all(i) < 0.00
  
       hoppel_min_idx(i)           = NaN;
       hoppel_min_meas_idx(i)      = NaN;
       hoppel_min_intersect_idx(i) = NaN;
       hoppel_min(i)               = NaN;
       hoppel_min_meas(i)          = NaN;
       hoppel_min_intersect(i)     = NaN;
       dHM_ait(i)                  = NaN;
       num_hoppel_min(i)                  = NaN;
       num_aitken_hoppel_min(i)           = NaN;
       num_accum_hoppel_min(i)            = NaN;
       else
       dHM_ait(i) = dHM_ait_all(i);
      
       % integrated size distribution from Hoppel minimum diameter
       % (effectively CCN #)
       num_hoppel_min(i)    = nansum(PNSD_sum(hoppel_min_idx(i):end,i) .* dlogDp);
       % integrated modal distributions from Hoppel miniumum diameter
       num_aitken_hoppel_min(i)    = sum(aitken_PNSD(hoppel_min_idx(i):end, i) .* dlogDp);
       num_accum_hoppel_min(i)     = sum(accum_PNSD(hoppel_min_idx(i):end, i) .* dlogDp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end
       end % if not a real Hoppel Minimum
   end % mode peak check for Hoppel Minimum
   end
  
 end % time loop

%% 

for i=1:length(category_SEMS_only)
    if category_SEMS_only(i)==1
        hoppel_min(i)=hoppel_min(i);
    else
        hoppel_min(i)=NaN;
    end
end
%% 
hoppel_min_SEMS=hoppel_min;
%% 
% 
% save('SEMS_only_hoppel_min_0807.mat','hoppel_min_SEMS','time_5min_final');

%% Concatenate SEMS only and Merged hoppel min

clear all; close all; clc;

load('Soledad_merged_hoppel_min.mat')
load('SEMS_only_hoppel_min_0807.mat')
%% 
hoppel_min_fitted=NaN(length(time_5min_final),1);

[intersect_time, time_5min_final_idx, time_merged_idx]=intersect(time_5min_final,time_merged);

for i=1:length(time_merged)
    k=time_merged_idx(i);
    j=time_5min_final_idx(i);
    hoppel_min_fitted(j)=hoppel_min(k);
end
count=0
for i=1:length(time_5min_final);
    if hoppel_min_SEMS(i)>0;
       hoppel_min_fitted(i)=hoppel_min_SEMS(i);
    end
end
%% saving hooppel min
% save('hoppel_min_fitted.mat','hoppel_min_fitted','time_5min_final');

        