% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase. The current version assumes
% that measurement_soa shifts on every exposure trial after encountering an
% SOA, whhile taking MCD_corr into account. It also assumes that
% measurement_SOA is distributed as Gaussian centered on remapped_SOA.

clear all; clc; close all;

%% fixed parameters

% set parameters
sim_trial = 10;
adaptor_soas = [-700, -300:100:300, 700]./1000; % s
n_soas = length(adaptor_soas);

% set session parameters
exposure_trial = 250;

% define fixed parameters to create time series;
% see parameter details in the function below
fs = 1e3; % hz
stim_dura = 0.033; % in sec
s_pad = 7.5; % in sec

%% free parameters to be fitted
bias = 0.02; % in sec
tv = 0.0873;
ta = 0.0684;
tav = 0.7859;
theta_av = 1;
learning_rate = 0.001;
sigma_soa= 0.2;

%% initiate recalibration effect for all soas
soa_recal = cell(1, n_soas);
last_recal = NaN(n_soas, sim_trial); % summarize the last recalibration effect

%% run simulation
for i = 1:n_soas
    soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t = 1:sim_trial

        soa_recal{i}(t,:) = update_recal_gaussian(exposure_trial, soa, ...
            fs, stim_dura, s_pad, ...
            bias, sigma_soa, learning_rate, ...
            ta, tv, tav);

        last_recal(i, t) = soa_recal{i}(t,end);

    end
end

%% summarize and plot

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