%% EPCAPE_SMPS_mode_fitting_v2024_0905.m
%
% PURPOSE: Fits two lognormal modes to 5 minute measured size
% distributions from SMPS during EPCAPE. 
%
% INPUTS: modal fit output from fitting model; requires the function
% create_psd.m in library

% OUTPUTS: Saves lognormal modal fits to a .mat file. 
%
% AUTHOR: Jeramy Dedrick
%         Scripps Institution of Oceanography
%         March 18, 2022
%
% Author: Atsushi O
%         SIO
%         July 19,2024
% Notes: Modified variable name for EPCAPE analysis of Pier distribution
% and added section for grpahing a well as odified it to only fit SMPS
% region (excluding mode1)
%% %% Clearing the variables, figures, command window, etc..
clear all; close all; clc

%% %  
% data folder
data    = '/Users/atsushiosawa/Downloads/'; %file name for the user

% output folder
out_dir = '/Users/atsushiosaa/Downloads/';  %file name for the user

cd '/Users/atsushiosawa/Downloads/EPCAPEAPSSMPS'
WD = cd;


%% load data

% data folder
data    = WD;
%% 


load(strcat(data, '/SMPS_PNSD_5min.mat'));

%% 


PNSD_merged_5min=SMPS_dndlodp_5min;

% % input size distribution (dN/dlog10Dp)
PNSD = PNSD_merged_5min.';  % [ROWS: time, COLUMNS: diameter]
% 
% % diameters (units = micrometers) 
D_merged=SMPS_PNSD_diam;
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
d_start_mode2 = 126;  % start diameter index for fit region of mode 2
d_end_mode2   = 159;  % end diameter index for fit region of mode 2
d_start_mode3 = 75;   % start diameter index for fit region of mode 3
d_end_mode3   = 125;  % end diameter index for fit region of mode 3
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
for ii = start:end
 

    disp(strcat('Fitting Size Distribution of smps:', num2str(ii), '/', num2str(length(time))))

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
save('pier_mode_fit_out_missing_sepandoct_nofilter_0722.mat',...
    'mode2_N','mode2_diam','mode2_gsd',...
    'mode3_N','mode3_diam','mode3_gsd',...
    'chi2_mode2','chi2_mode3',...
    'SMPS_datenum')
%% %% We can start with this section if we have already run the code above
%don't need to run above if all the line above was run
clear all; close all; clc
data    = '/Users/atsushiosawa/Downloads/'; %file name for the user

% output folder
out_dir = '/Users/atsushiosaa/Downloads/';  %file name for the user


cd '/Users/atsushiosawa/Downloads/EPCAPEAPSSMPS'
WD = cd;


%% load data

% data folder
data    = WD;


% read data
load(strcat(data, '/SMPS_PNSD_5min.mat'));
time_5min=SMPS_datenum;
PNSD_merged_5min=SMPS_dndlodp_5min;
D_merged=SMPS_PNSD_diam;
% input size distribution (dN/dlog10Dp)
PNSD = PNSD_merged_5min.';  % [ROWS: time, COLUMNS: diameter]

% diameters (units = micrometers) 
D = D_merged; % [ROWS: diameters, COLUMNS: 1]

% "time" as the number of observations
time = 1:size(PNSD,1);
% 
% load('CPCfilter_atpier_mode_fit_out_missing_sepandoct_n5_625.mat')
%change above for different.mat output from running the fitting model
load('EPCAPE_SMPS_APS_merged_5min_20240717.mat')
load('pier_mode_fit_out_missing_sepandoct_nofilter_0722.mat')

% indices of diameter bins

d_start_mode2 = 32;  % start diameter index for fit region of mode 2
d_end_mode2   = 49;  % end diameter index for fit region of mode 2
d_start_mode3 = 5;   % start diameter index for fit region of mode 3
d_end_mode3   = 31;  % end diameter index for fit region of mode 3

%% %% Code for graphing purposes
D=SMPS_PNSD_diam;

[size_dist2,mass_dist2] = create_PSD(D, mode2_N, mode2_diam, mode2_gsd);
[size_dist3,mass_dist3] = create_PSD(D, mode3_N, mode3_diam, mode3_gsd);


sum_dist = size_dist2 + size_dist3;
dlogdp=log10(D_merged(50))-log10(D_merged(49)); %log 10 difference of successive diameters
% %% %% Percent difference between measured and fitted mass concentration
% 
% mass_dist_fitted=mass_dist2+mass_dist3;
% mass_fitted=sum(mass_dist_fitted.*dlogdp);%value we want
% 
% mass_dist_measured=PNSD'.* ...
%             (pi / 6) .* ...
%             (1) .* ...
%             (D .* 1e-6) .^3 .* ...
%             (1e6 .* 1e6) .* 1e6;   %measured size distirbution computing to get the mass here
% 
% mass_meas_conc = nansum(mass_dist_measured.*dlogdp);%value we want
% 
% mass_concent_perc= (abs(mass_fitted-mass_meas_conc)./((mass_fitted+mass_meas_conc)/2))*100;
%% %% percent difference between measured and fitted number concentration

% tot_fit_numb_conc=mode1_N+mode2_N+mode3_N;
% meas_int_numb_concent=Nt_merged_5min;
% numb_concent_percent= (abs(tot_fit_numb_conc-meas_int_numb_concent)./((tot_fit_numb_conc+meas_int_numb_concent)/2))*100;
%% %% Plotting 2 by 2 plot of PNSD distirbution and fits
for i=67428:4:67431 % we can change these ranges
    idxplots=i:i+3;
    for j=1:4
    subplot(2,2,j)
    figure(1)
    idx=idxplots(j);
    loglog(D,movmean(PNSD(idx,:), 8), 'ok')
        hold on

       
        loglog(D, size_dist2(:,idx), 'r','LineWidth', .1) %correspond to mode 2
        loglog(D, size_dist3(:,idx),'g','LineWidth',.1)%correspond mode 3
        % loglog(D, sum_dist(:,idx), 'k', 'LineWidth', .1)% We remove
        % the percentage here if we desire to include the sum fit curve
        % to the plot
       
        % annotation on the plot for the fit bounds we set when running
        % the fitting model
        x1=xline(2e-0,'k--','sea ub');
        x1.LabelVerticalAlignment = 'bottom';
        x1.LabelHorizontalAlignment = 'center';
        x2=xline(4e-1,'k--','sea lb/acc ub');
        x2.LabelVerticalAlignment = 'bottom';
        x2.LabelHorizontalAlignment = 'center';
        x3=xline(6e-2,'k--','acc lb/ait ub');
        x3.LabelVerticalAlignment = 'bottom';
        x3.LabelHorizontalAlignment = 'center';
        x4=xline(1e-2,'k--','ait lb');
        x4.LabelVerticalAlignment = 'bottom';
        x4.LabelHorizontalAlignment = 'center';
        x5=xline(8e-1,'k--','sea x0');
        x5.LabelVerticalAlignment = 'bottom';
        x5.LabelHorizontalAlignment = 'center';
        x6=xline(1e-1,'k--','acc x0');
        x6.LabelVerticalAlignment='bottom';
        x6.LabelHorizontalAlignment = 'center';
        x7=xline(3e-2,'k--','ait x0');
        x7.LabelVerticalAlignment='bottom'
        x7.LabelHorizontalAlignment = 'center';

        hold off
        
   
        
        ylim([1e-2 1e5])%1e4
        xlim([1e-2 1e3])%1e2
        xticks([1e-2 1e-1 1e-0 1e1 1e2])
        yticks([1e-0 1e1 1e2 1e3 1e4])
        axis square
        
        % a=mass_concent_perc(idx)
        % b=numb_concent_percent(idx)
        % txt=['Mass % diff:' num2str(round(a))]
        % text(2,1000,txt)
        % tot=['Numb % diff:' num2str(round(b))]
        % text(2,1700,tot)
        title(num2str(idx))
    end
    % change the name of directory which these plot would be saved
    % below;
    %when we want to save these plots, remove percentage from two lines
    %with printname and print(printname,....)
    % **Note:everytime model is run, what was there is overwritten
    printname=strcat('/Users/atsushiosawa/Downloads/EPCAPE_PIER_APSSMPS_missingmonth/check_fitting',num2str(idxplots(1)),'_',num2str(idxplots(end)));
    print(printname,'-dpng')%png don't stretch when creating 2x2 graphs... and print from the original output
end






