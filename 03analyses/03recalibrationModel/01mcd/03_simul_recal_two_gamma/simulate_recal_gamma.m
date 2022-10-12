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
s_pad = 7.5; % in sec

%% free parameters to be fitted later
bias = 0.02; % in sec
ka = 2;
theta_a = 0.5;
kv = 2;
theta_v = 0.5;
kav = 2;
theta_av = 0.5;
learning_rate = 0.001;

%% initiate recalibration effect for all soas
soa_recal = cell(1, n_soas);
last_recal = NaN(n_soas, simTrial); % summarize the last recalibration effect

%% run simulation
for i = 1:n_soas
    soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t = 1:simTrial

        % model 1: soa_m ~ N(soa', sigma)
        recal_effect = update_recal_gamma(exposure_trial, soa, fs, stim_dura, s_pad,...
    ka, theta_a, kv, theta_v, kav, theta_av, bias, learning_rate);


        last_recal(i, t) = soa_recal{i}(t,end);

    end
end

%% plot

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