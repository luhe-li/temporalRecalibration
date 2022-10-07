% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase. The current version assumes
% that measurement_soa shifts on every exposure trial after encountering an
% SOA, whhile taking MCD_corr into account.

clear all; clc; close all;

%% fixed parameters

% set parameters
simTrial = 1;
adaptor_soas = [-700, -300:100:300, 700]./1000; % s
n_soas = length(adaptor_soas);

% set session parameters
exposure_trial = 250;
% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase. The current version assumes
% that measurement_soa shifts on every exposure trial after encountering an
% SOA, whhile taking MCD_corr into account.

clear all; clc; close all;

%% fixed parameters

% set parameters
simTrial = 1;
adaptor_soas = [-700, -300:100:300, 700]./1000; % s
n_soas = length(adaptor_soas);

% set session parameters
exposure_trial = 250;


% define fixed parameters to create time series; see parameter
% see details in the function below
duration = 16; % in s
fs = 1e3; % hz
onset = 9; % in sec
stim_dura = 0.033; % in sec

%% free parameters to be fitted later
% see parameter details in the function below
ta = 0.078;
tv = 0.068;
tav = 0.786;
learning_rate = 0.01;

% For each stimulus with an fixed soa (soa = t_a - t_v), participants formed a
% noisy sensory measurement m_soa with a bias in the head. We assume that
% m_soa is sampled from N(soa + bias, sigma).
bias = 0.02; % in sec
sigma_soa = 0.05; % in sec

% initiate recalibration effect for all soas
soa_recal = cell(1, n_soas);
last_recal = NaN(n_soas, simTrial); % summarize the last recalibration effect

%% run simulation
for i = 1:n_soas
    soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t = 1:simTrial

        % model 1: soa_m ~ N(soa', sigma)
        soa_recal{i}(t,:)= update_recal(exposure_trial, soa, bias, sigma_soa,...
            duration, fs, onset, stim_dura,...
            ta, tv, tav,...
            learning_rate);

        last_recal(i, t) = soa_recal{i}(t,end);

    end
end

%%

last_recal_iSOA_error = std(last_recal, [], 2);
last_recal_iSOA = mean(last_recal, 2);

figure; hold on
box off
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 500, 400]);
errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2)
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('Adaptor SOA (s)')
ylabel('recalibration (s)')

%
% % plot simulated recalibration effect (mu_shift)
% figure; hold on
% last_recal_effect = recal_effect(:, end);
% histogram(last_recal_effect)

% define fixed parameters to create time series; see parameter
% see details in the function below
duration = 16; % in s
fs = 1e3; % hz
onset = 9; % in sec
stim_dura = 0.033; % in sec

%% free parameters to be fitted later
% see parameter details in the function below
ta = 0.078;
tv = 0.068;
tav = 0.786;
learning_rate = 0.01;

% For each stimulus with an fixed soa (soa = t_a - t_v), participants formed a
% noisy sensory measurement m_soa with a bias in the head. We assume that
% m_soa is sampled from N(soa + bias, sigma).
bias = 0.02; % in sec
sigma_soa = 0.05; % in sec

% initiate recalibration effect for all soas
soa_recal = cell(1, n_soas);
last_recal = NaN(n_soas, simTrial); % summarize the last recalibration effect

%% run simulation
for i = 1:n_soas
    soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t = 1:simTrial

        % model 1: soa_m ~ N(soa', sigma)
        soa_recal{i}(t,:)= update_recal(exposure_trial, soa, bias, sigma_soa,...
            duration, fs, onset, stim_dura,...
            ta, tv, tav,...
            learning_rate);

        last_recal(i, t) = soa_recal{i}(t,end);

    end
end

%%

last_recal_iSOA_error = std(last_recal, [], 2);
last_recal_iSOA = mean(last_recal, 2);

figure; hold on
box off
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 500, 400]);
errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2)
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('Adaptor SOA (s)')
ylabel('recalibration (s)')

%
% % plot simulated recalibration effect (mu_shift)
% figure; hold on
% last_recal_effect = recal_effect(:, end);
% histogram(last_recal_effect)