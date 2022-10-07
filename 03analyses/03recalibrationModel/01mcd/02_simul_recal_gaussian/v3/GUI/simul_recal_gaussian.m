% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase. The current version assumes
% that measurement_soa shifts on every exposure trial after encountering an
% SOA, while taking MCD_corr into account. It also assumes that
% measurement_SOA is distributed as Gaussian centered on remapped_SOA.

clear all; clc; close all; %parpool(10); %50 hours per job

%% fixed parameters

% set parameters
sim_trial    = 10;
adaptor_soas = [-700, -300:100:300, 700]./1000; % s
n_soas       = length(adaptor_soas);

% set session parameters
exposure_trial = 250;

% define fixed parameters to create time series;
% see parameter details in the function below
fs          = 1e3; % hz
stim_dura   = 0.033; % in sec
ts_duration = 16; % in sec, 7 x 2 padding + 2 stimulus duration

%% free parameters to be fitted

% bias          = linspace(-0.2, 0.2, 10);%0.02; % in sec
% ta            = linspace(0.02, 0.15, 10); %0.0684;
% tv            = linspace(0.05, 0.2, 10);%0.0873;
% tav           = linspace(0.1,1.5,10); %0.7859;
% sigma_soa     = linspace(0.001, 0.2,10);%0.2;

% new range based on min max +- se
bias = linspace(-0.2, 0.2, 10);%0.02; % in sec
tv = linspace(0.0991, 0.1559, 10); %0.0873;
ta = linspace( 0.0448, 0.1152, 10); %0.0684;
tav = linspace(0.3104, 1.2086, 10);%0.7859;
sigma_soa     = linspace(0.001, 0.2,10);%0.2;

%write all the combinations
MAT_comb      = combvec(bias, ta, tv, tav, sigma_soa);
n_comb        = size(MAT_comb,2);

theta_av      = 1;
learning_rate = 0.001;

% initiate time series for the filter
% t should be the same length as the signal time series
nsample      = ts_duration * fs + 1; 
t            = linspace(0, ts_duration, nsample);

%% run simulation
%deal with sliced variables: 
%https://www.mathworks.com/help/parallel-computing/troubleshoot-variables-in-parfor-loops.html
idx_reshape = reshape(1:n_comb,[],10); %split it into 10 jobs.
job_idx     = 1;
nn_select   = idx_reshape(:, job_idx)';

n_comb_select = length(nn_select);
soa_recal   = NaN(n_comb_select, n_soas, sim_trial, exposure_trial+1);
last_recal  = NaN(n_comb_select, n_soas, sim_trial); % summarize the last recalibration effect


parfor nn = nn_select
    disp(nn)
    MAT_comb_nn  = MAT_comb(:,nn);
    bias_nn      = MAT_comb_nn(1);
    % low-pass filters (Equation 1)

    fa_nn        = fft(t.*exp(-t./MAT_comb_nn(2)));
    fv_nn        = fft(t.*exp(-t./MAT_comb_nn(3)));
    fav_nn       = fft(t.*exp(-t./MAT_comb_nn(4)));
    sigma_soa_nn = MAT_comb_nn(5);
    
    %initialize within the parfor loop
    soa_recal_temp  = NaN(n_soas, sim_trial, exposure_trial+1);
    last_recal_temp = NaN(n_soas, sim_trial);
    for i = 1:n_soas
        soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session

        for j = 1:sim_trial
            soa_recal_temp(i,j,:) = update_recal_gaussian(exposure_trial, soa, ...
                ts_duration, fs, stim_dura, fa_nn, fv_nn, fav_nn, bias_nn, ...
                sigma_soa_nn, learning_rate);

            last_recal_temp(i,j) = soa_recal_temp(i,j,end);
        end
    end
    soa_recal(nn,:,:,:) = soa_recal_temp;
    last_recal(nn,:,:)  = last_recal_temp;
    disp(last_recal_temp);
end

%% save data
SimResults = {sim_trial, adaptor_soas, exposure_trial,fs,stim_dura, MAT_comb, soa_recal, last_recal};
save(['SimResults_', num2str(job_idx), '_', datestr(datetime('now')), '.mat'],'SimResults');


