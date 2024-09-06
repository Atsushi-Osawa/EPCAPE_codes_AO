%% EPCAPE_SEMS_mode_fitting_v2024_0906.m
%
% PURPOSE: Fits two lognormal modes to 5 minute-averaged measured size
% distributions from SEMS during EPCAPE. 
%
% INPUTS: SEMS size distribution; 

% OUTPUTS: Saves lognormal modal fits to a .mat file. 
%
% AUTHOR: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         March 18, 2022
%
% Author: Atsushi O
%         SIO
%         July 19,2024
% Notes: Modified variable name for EPCAPE analysis of Mt. Soledad distribution
% and added section for grpahing a well as odified it to only fit SMPS
% region (excluding mode1)

%% %% Clearing the variables, figures, command window, etc..
clear all; close all; clc



%% %  
% data folder
data    = '/Users/atsushiosawa/Downloads/'; %file name for the user

% output folder
out_dir = '/Users/atsushiosawa/Downloads/';  %file name for the user

% identifies files 
SEMS_files = dir(fullfile(strcat(data,'/SEMS_data_out/'),'EPCAPE*.nc'));

% identifies files 
APS_files = dir(fullfile(strcat(data,'/APS_data_out'),'EPCAPE*.nc'));


%% Read in data

disp(strcat('Opening and reading SEMS & APS files'))


% empty matrices for concatenation
SEMS_time = [];
SEMS_PNSD = [];
SEMS_Ntot = [];
SEMS_qc   = [];
SEMS_cvi  = [];

APS_time = [];
APS_PNSD = [];
APS_Ntot = [];
APS_qc   = [];
%% 

for i = 1:length(SEMS_files)

    % creates file strings
    SEMS_file = strcat(SEMS_files(i).folder, '/', SEMS_files(i).name);
    APS_file  = strcat(APS_files(i).folder, '/', APS_files(i).name); 

    % read data
    SEMS_time_r = ncread(SEMS_file, 'time')';
    SEMS_PNSD_r = ncread(SEMS_file, 'PNSD');
    SEMS_Ntot_r = ncread(SEMS_file, 'Ntot')';
    SEMS_d_mid  = ncread(SEMS_file, 'diam_mid');
    SEMS_qc_r   = ncread(SEMS_file, 'qc')';
    SEMS_cvi_r  = ncread(SEMS_file, 'cvi_flag')';
    

    APS_time_r = ncread(APS_file, 'time')';
    APS_PNSD_r = ncread(APS_file, 'PNSD');
    APS_Ntot_r = ncread(APS_file, 'Ntot')';
    APS_d_mid  = ncread(APS_file, 'diam_mid');
    APS_qc_r   = ncread(APS_file, 'APS_qc')';

    % quality control (only want QC = 1)
    % SEMS_PNSD_r(:,SEMS_qc<1) = NaN;
    % SEMS_Ntot_r(SEMS_qc<1)   = NaN;
    % APS_PNSD_r(:,APS_qc<1)  = NaN;
    % APS_Ntot_r(APS_qc<1)    = NaN;
% 62431,27601,12266,101,,0,129
    % concatenate
    SEMS_time = [SEMS_time SEMS_time_r];
    SEMS_PNSD = [SEMS_PNSD SEMS_PNSD_r];
    SEMS_Ntot = [SEMS_Ntot SEMS_Ntot_r];
    SEMS_qc   = [SEMS_qc SEMS_qc_r];
    SEMS_cvi  = [SEMS_cvi, SEMS_cvi_r];

    APS_time = [APS_time APS_time_r];
    APS_PNSD = [APS_PNSD APS_PNSD_r];
    APS_Ntot = [APS_Ntot APS_Ntot_r];
    APS_qc   = [APS_qc APS_qc_r];

end

%% sort so that time increases
[SEMS_time_sort, SEMS_time_sort_idx] = sort(SEMS_time, 'ascend');
[APS_time_sort, APS_time_sort_idx]   = sort(APS_time, 'ascend');

SEMS_time = SEMS_time(SEMS_time_sort_idx);
SEMS_PNSD = SEMS_PNSD(:,SEMS_time_sort_idx);
SEMS_Ntot = SEMS_Ntot(SEMS_time_sort_idx);
SEMS_cvi  = SEMS_cvi(SEMS_time_sort_idx);
SEMS_qc   = SEMS_qc(SEMS_time_sort_idx);

APS_time = APS_time(APS_time_sort_idx);
APS_PNSD = APS_PNSD(:,APS_time_sort_idx);
APS_Ntot = APS_Ntot(APS_time_sort_idx);
APS_qc   = APS_qc(APS_time_sort_idx);



%% 


PNSD_merged_5min=SEMS_PNSD;
D_merged=SEMS_d_mid;
% 
% % input size distribution (dN/dlog10Dp)
PNSD = PNSD_merged_5min.';  % [ROWS: time, COLUMNS: diameter]
% 
% % diameters (units = micrometers) 
D = D_merged; % [ROWS: diameters, COLUMNS: 1]
% 
% % "time" as the number of observations
time = 1:size(PNSD,1);
% %% Initial conditions and fitting options
% 
% % least-square fitting options
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
% 

% least squared fit bounds for middle mode (N, diam, gsd)
mode2_lb      = [0 0.06 1.2];   % lower-bound
mode2_ub      = [8000 0.4 2.0]; % upper-bound
mode2_x0      = [500 0.1 1.5]; % inital guess

% least squared fit bounds for smallest mode (N, diam, gsd)
mode3_lb      = [0 0.01 1.2];    % lower-bound
mode3_ub      = [8000 0.06 1.6]; % upper-bound
mode3_x0      = [500 0.03 1.5];  % inital guess

% % indices of diameter bins
% d_start_mode1 = 55;  % start diameter index for fit region of mode 1
% d_end_mode1   = 95;  % end diameter index for fit region of mode 1
d_start_mode2 = 29;  % start diameter index for fit region of mode 2
d_end_mode2   = 53;  % end diameter index for fit region of mode 2
d_start_mode3 = 1;   % start diameter index for fit region of mode 3
d_end_mode3   = 28;  % end diameter index for fit region of mode 3
% 
% 
% % number of datapoints in moving average for smoothing before fitting
n_avg    = 5; 

% counters
k        = 1; 
counter  = 0;

% non-linear fit parameters (for unconstrained fitting)
opts              = statset('nlinfit');
opts.Robust       = 'on';
opts.RobustWgtFun = 'bisquare';
opts.MaxIter      = 1000;
opts.Display      = 'off';


%% 
%% Storing fit parameters for modes

mode2_N    = NaN(1,length(time));
mode2_diam = NaN(1,length(time));
mode2_gsd  = NaN(1,length(time));

mode3_N    = NaN(1,length(time));
mode3_diam = NaN(1,length(time));
mode3_gsd  = NaN(1,length(time));

chi2_mode2      = NaN(1,length(time));
chi2_mode3        = NaN(1,length(time));


%% Fitting model 

% Eq. (1) Modini et al. (2015)           
modelfun = @(b,d)(b(1) ./ ...
                    (sqrt(2 .* pi) .* log10(b(3))) .* ...
                    exp(-((log10(d)-log10(b(2))).^2 ) ./ (2 .* log10(b(3)).^2))); 



%% %% Running of the fitting model
for ii = 1%start_idx:end
    disp(ii)
    %we could modify above to run the fitting model to be of certain range    disp(strcat('Fitting Size Distribution of smps:', num2str(ii), '/', num2str(length(time))))

    try
            counter = counter + 1;    


            % accum mode using constraints (largest mode)
            [xmode2{ii}, resnorm_mode2{ii}, residual_lsq_mode2{ii}, ...
             ~, ~, ~, J_lsq_mode2{ii}] = lsqcurvefit(modelfun, ...
                                               mode2_x0, ...
                                               D(d_start_mode2:d_end_mode2)', ...
                                               movmean(PNSD(ii,d_start_mode2:d_end_mode2),n_avg), ...
                                               mode2_lb,mode2_ub, ...
                                               options);

            % removing accum mode from size distribution
            PNSD_mode3_remainder(ii,:) = PNSD(ii,:) - modelfun(xmode2{ii},D)'; % this simply removes the distirbution associated with mode 1


            % ait mode using constraints (third largest mode)
            [xmode3{ii}, resnorm_mode3{ii}, residual_lsq_mode3{ii}, ...
             ~, ~, ~, J_lsq_mode3{ii}] = lsqcurvefit(modelfun, ...
                                       mode3_x0, ...
                                       D(d_start_mode3:d_end_mode3)', ...
                                       movmean(PNSD_mode3_remainder(ii,d_start_mode3:d_end_mode3),n_avg), ...
                                       mode3_lb,mode3_ub, ...
                                       options);



            % saving modal fit data
        
            mode2_N(ii)    = xmode2{ii}(1);
            mode2_diam(ii) = xmode2{ii}(2);
            mode2_gsd(ii)  = xmode2{ii}(3);

            mode3_N(ii)    = xmode3{ii}(1);
            mode3_diam(ii) = xmode3{ii}(2);
            mode3_gsd(ii)  = xmode3{ii}(3);

            % measured size distribution
            meas          = PNSD(ii,:)';
            meas(meas==0) = NaN;
            % get the fit of each mode using model
         
            mode2_fit_model       = modelfun([mode2_N(ii), mode2_diam(ii), mode2_gsd(ii)], D);
            mode3_fit_model       = modelfun([mode3_N(ii), mode3_diam(ii), mode3_gsd(ii)], D);
            % squared residual
           
            SQR_mode2      = (mode2_fit_model(d_start_mode2:d_end_mode2) - meas(d_start_mode2:d_end_mode2)).^2;
            SQR_mode3       = (mode3_fit_model(d_start_mode3 :d_end_mode3) - meas(d_start_mode3:d_end_mode3)).^2;
            % % chi square
           
            chi2_mode2(ii) = nansum( SQR_mode2 ./ meas(d_start_mode2:d_end_mode2) );
            chi2_mode3(ii) = nansum( SQR_mode3 ./ meas(d_start_mode3:d_end_mode3) );

            
    catch
    end
end

  
%% 


save('file_name',...
    'mode2_N','mode2_diam','mode2_gsd',...
    'mode3_N','mode3_diam','mode3_gsd',...
    'chi2_mode2','chi2_mode3',...
    'SEMS_time','SEMS_d_mid','SEMS_PNSD','SEMS_Ntot','SEMS_cvi')
