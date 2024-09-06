%% EPCAPE_SEMS_APS_merge_v20231011.m
%
% PURPOSE: Reads in EPCAPE SMPS and APS data, computes averages, and merges
% the size distributions following Khlystov (2004) and Modini et al.
% (2015). 
%
% INPUTS: 5 minute SMPS and APS data.
%
% OUTPUTS: .MAT files of merged size distributions and CVI size
% distributions. 
%
% AUTHOR: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         July 24, 2023
%         
% EDITED: Ian Marroquin 
%         Washington University in St. Louis
%         August 15, 2023
%         
% 
% Note 1: I adjusted this file such that the input is already averaged hourly
% size distributions. Here, only the diameter is used from the febmar23.nc
% files. Also, even though the variables correspond to hourly distributions
% they are named _15min etc. Finally, line 101 is hard coded (1008 hourly size
% distributions translate to 42 days), yet may be easily corrected/adjusted.

% Note 2: The sems_hourly.csv and aps_hourly.csv are essentially generated
% by combining the individual sems files (see under the averaging_function_gamma)
% into one larger .csv file (same procedure for aps).
%
%
% EDITED: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         August 16, 2023
%
% Note 1: Re-added averaging of size distributions. Set up code to search
% files within processing directory for SEMS and APS. 
%
%
% EDITED: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         October 11, 2023
%
% Note 1: added 15 min average mergers. Added 5-min averages (Oct. 14)
%
%
% EDITED: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         November 13, 2023
%
% Note 1: added lines to remove negative data points in merged size
% distribution. 
%
% EDITED: Atsushi
%         Scripps Institution of Oceanography
%         November 13, 2023
%
% Note 1: Modified code to merge ARM SMPS and APS size distributions at 5
% min average.
% EDITED: ATSUSHI 
%         SIO
%         July 16th, 2023 
% Adjusted diameter bins to match that of APS+SEMS code that exists




clear all; close all; clc

cd '/Users/atsushiosawa/Downloads/EPCAPEAPSSMPS'
WD = cd;


%% load data

% data folder
data    = WD;


% read data
load(strcat(data, '/SMPS_PNSD_5min.mat'));
load(strcat(data, '/APS_PNSD_5min.mat'));



%% Intersect Times


[intersect_time, APS_intersect_time_idx, SMPS_intersect_time_idx] = intersect(APS_datenum,SMPS_datenum);

APS_PNSD_5min=APS_dndlodp_5min(:,APS_intersect_time_idx);
SMPS_PNSD_5min=SMPS_dndlodp_5min(:,SMPS_intersect_time_idx);

flag=NaN(2,length(intersect_time)); % flag varaible

%% Merge 5-min SMPS and APS
clc

%%%%%%%%%%%%%%%%% Setting/Getting Parameters and Indices %%%%%%%%%%%%%%%%%%

% log-spaced diameter grid to merge over (10  nm - 10 micron)
D_merged = logspace(log10(0.011), log10(10), 100)';

% dlogDp of diameters
dlogDp_merged = log10(D_merged(50)) - log10(D_merged(49));

% indices of SMPS tail (visual)
SMPS_tail_idx  = 155:171; 

% diameters of SMPS tail
SMPS_d_tail    = SMPS_PNSD_diam(SMPS_tail_idx); 

% smoothing window (points)
SMPS_smooth_window = 3;
APS_smooth_window  = 3;

% indices for which to grab SEMS for merger with APS
SMPS_d_to_merge_idx = 66:171; %1:58; 

% indices of APS that overlap with the SEMS
overlap_idx = 7:12;
overlap_d   = APS_PNSD_diam(overlap_idx);

% SMPS extended diameter for overlap
SMPS_extended_d = [APS_PNSD_diam(9:14)];

% diameters of SMPS that merge with APS
% SMPS_d_merge_idx = 55:60;
% SMPS_d_merge = SMPS_PNSD_diam(SMPS_d_merge_idx);

% density values to test (g cm-3) (Modini et al., 2015)
rho_test = 0.9:0.01:3.0;

% shape factor
chi_shape = 1;

% reference density (g cm-3)
rho_0 = 1;

%%%%%%%%%%%%%%%%%% Merging the size distributoins %%%%%%%%%%%%%%%%%%%%%%%%%
%% 
S2min = NaN(1,length(intersect_time))

for i = 1:length(intersect_time)

    disp(strcat('Merging SMPS+APS, 5-min: ', num2str(i), '/', num2str(length(intersect_time))))

    % check if size distributions exist for the specified time period
    SMPS_check      = all(isnan(SMPS_PNSD_5min(:,i)));
   
    APS_check       = all(isnan(APS_PNSD_5min(:,i)));
  
    if APS_check>0 & SMPS_check>0;
        flag(2,i)=3;
    elseif APS_check>0;
        flag(2,i)=2;
    elseif SMPS_check>0;
        flag(2,i)=1;
    elseif APS_check==0 & SMPS_check==0
        flag(2,i)=0;
    end;
    size_dist_check = SMPS_check + APS_check;

   

    if size_dist_check > 0 % at least one size distribution missing

        PNSD_merged_5min(:,i) = NaN(length(D_merged),1);
        Nt_merged_5min(:,i)   = NaN(1,1);
        rho_merge_5min(:,i)   = NaN(1,1);
        flag(1,i)=1;

    else

        % fix weird bins between 0.5 and 1 micron by spline interpolation 
        APS_interp_sizes_r         = APS_PNSD_5min(:,i);
        % APS_interp_sizes_r = interp1(APS_d_mid, ...
        %                              APS_PNSD_5min(7:,i), ...
        %                              SMPS_PNSD_diam()
        % APS_interp_sizes(2:13,:) = NaN;
        % APS_interp_sizes(isnan(APS_interp_sizes)) = interp1(APS_d_mid(~isnan(APS_interp_sizes)), ...
        %                                                     APS_interp_sizes(~isnan(APS_interp_sizes)), ...
        %                                                     APS_d_mid(isnan(APS_interp_sizes)), 'spline', 'extrap');

        % SMPS size distribution in tail
        SMPS_PNSD_tail = SMPS_PNSD_5min(SMPS_tail_idx,i);

        % Fit power law to tail of SMPS
        [SMPS_fit, SMPS_fit_gof] = fit(SMPS_PNSD_diam(SMPS_tail_idx), ... 
                                       movmean(SMPS_PNSD_tail, SMPS_smooth_window), ... % tail of SEMS smoothed
                                       'power1');

        SMPS_fit_params = coeffvalues(SMPS_fit); % coefficients from fit
        SMPS_fit_r2     = SMPS_fit_gof.rsquare;  % fit R^2


        if SMPS_fit_r2 < 0.50 % SMPS too noisy to fit a power law tail (subjective R2)

            PNSD_merged_5min(:,i) = NaN(length(D_merged),1);
            Nt_merged_5min(:,i)   = NaN(1,1);
            rho_merge_5min(:,i)   = NaN(1,1);

            flag(1,i)=2;

        elseif SMPS_fit_params(1) == 0 % if there's no fit

            PNSD_merged_5min(:,i) = NaN(length(D_merged),1);
            Nt_merged_5min(:,i)   = NaN(1,1);
            rho_merge_5min(:,i)   = NaN(1,1);

            flag(1,i)=3;

        else

            S2 = [];
            % Test density  
            for k = 1:length(rho_test)

                % convert aerodynamic diameter to geometric physical diameter
                APS_d_phys = SMPS_extended_d ./ ...
                             sqrt(rho_test(k) / (chi_shape .* rho_0));

                % APS_d_phys = APS_d_mid(overlap_idx) ./ ...
                %              sqrt(rho_test(k) / (chi_shape .* rho_0));

                % compute objective function in overlap region (S^2, Khlystov 2004)
                S2(:,k) = sum((log10( SMPS_fit_params(1) .* APS_d_phys .^ SMPS_fit_params(2) ) - ...
                          log10(movmean(APS_PNSD_5min(overlap_idx,i), APS_smooth_window))).^2);% ./ ...
                          % (length(overlap_idx)-1);

                % objective function not a real value          
                if isinf(S2(:,k))

                    S2(:,k) = NaN;

                    flag(1,i)=4;

                else

                    S2(:,k) = S2(:,k);

                end  

            end

            % stop the merger if can't find a minimum of objective
            % function
            if all(isnan(S2))
    
                PNSD_merged_5min(:,i) = NaN(length(D_merged),1);
                Nt_merged_5min(:,i)   = NaN(1,1);
                rho_merge_5min(:,i)   = NaN(1,1);
                flag(1,i)=5;
            else
                flag(1,i)=0;
    
            % minimum value of objective function
            S2min_5min(i) = nanmin(S2);

            % density at minimum value of objective function
            rho_merge_5min(i) = rho_test(S2 == nanmin(S2));

            % shifted APS diameter using effective density
            APS_diam_new =  APS_PNSD_diam ./ ...
                            sqrt(rho_merge_5min(i) / (chi_shape .* rho_0));

            % APS diameters to merge APS with SEMS (find between
            % 0.7 to 11 micorn
            % APS_d_merge_idx = find(APS_diam_new >= 0.75 & ...
            %                        APS_diam_new <= 11);

            % APS diameters to merge APS with SEMS (find between
            % 0.7 to 11 micorn
            APS_d_merge_idx = find(APS_diam_new>=max(SMPS_PNSD_diam(SMPS_d_to_merge_idx)), 1, 'first'):length(APS_diam_new);


            % merged SEMS and shifted APS diameters
            d_merged_SMPS_APS = [SMPS_PNSD_diam(SMPS_d_to_merge_idx); APS_diam_new(APS_d_merge_idx)];

            % merged size distributions 
            PNSD_merged_r = [SMPS_PNSD_5min(SMPS_d_to_merge_idx,i); APS_PNSD_5min(APS_d_merge_idx,i)];

            % interpolating merged size distribution onto diameter grid
            PNSD_merged_5min(:,i) = interp1(d_merged_SMPS_APS, PNSD_merged_r, D_merged);

            % merged total number concentration
            Nt_merged_5min(:,i) = nansum(PNSD_merged_5min(:,i) .* dlogDp_merged,1);

            end

        end

    end

end


%% Remove Merged Distributions where there are negative values

% find negative value size distributions
[neg5_r, neg5_c]     = find(PNSD_merged_5min < 0);


% remove negative value size distributions
PNSD_merged_5min(:,neg5_c)   = NaN;
Nt_merged_5min(:,neg5_c)    = NaN;
rho_merge_5min(:,neg5_c)    = NaN


%% Check how the merged size distribution compares to measured

close all; clc

% index/indices for checking merged size distribution
idx = 235;%1257, 1810, 2350, 3105, 3260, 3264

loglog(SMPS_PNSD_diam, SMPS_PNSD_5min(:,idx), 'ok', 'MarkerFaceColor', 'k')
hold on
loglog(APS_PNSD_diam, APS_PNSD_5min(:,idx), 'or', 'MarkerFaceColor', 'r')
loglog(D_merged, PNSD_merged_5min(:,idx), '-og', 'MarkerFaceColor', 'g')

% xline(SEMS_d_mid(52), 'Color', 'k');
xline(SMPS_PNSD_diam(155), 'Color', 'k');
xline(SMPS_PNSD_diam(171), 'Color', 'k');
xline(SMPS_PNSD_diam(165), '--', 'Color', 'k');
% xline(APS_d_mid(1), 'Color', 'r');
xline(APS_PNSD_diam(9), 'Color', 'r');
xline(APS_PNSD_diam(13), 'Color', 'r');

hold off
ylim([1e-3 1e4])
xlim([1e-3 1e1])
legend('SMPS', 'APS', 'merged')

title({datestr(intersect_time(idx)) strcat('\rho = ', num2str(rho_merge_5min(idx)))})



%% Saving data
clc

% save_filename_5min = strcat(out_dir, 'EPCAPE_SMPS_APS_merged_5min_', datestr(now, 'yyymdd')8 '.mat');

save('EPCAPE_SMPS_APS_merged_5min_20240822.mat', ...
    'PNSD_merged_5min', ...
    'Nt_merged_5min', ...
    'D_merged', ...
    'rho_merge_5min', ...
    'intersect_time',...
    'flag','S2min_5min')

%should have saved flag

disp('Data saved')



