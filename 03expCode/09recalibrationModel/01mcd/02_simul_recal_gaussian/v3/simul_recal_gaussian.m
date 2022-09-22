% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase. The current version assumes
% that measurement_soa shifts on every exposure trial after encountering an
% SOA, while taking MCD_corr into account. It also assumes that
% measurement_SOA is distributed as Gaussian centered on remapped_SOA.

clear all; clc; close all;

%% fixed parameters

% set parameters
sim_trial = 1000;
adaptor_soas = [-700, -300:100:300, 700]./1000; % s
n_soas = length(adaptor_soas);

% set session parameters
exposure_trial = 250;

% define fixed parameters to create time series;
% see parameter details in the function below
fs = 1e3; % hz
stim_dura = 0.033; % in sec
ts_duration = 16; % in sec, 7 x 2 padding + 2 stimulus duration

%% free parameters to be fitted

bias = 0.02; % in sec
% tv = 0.0873;
tv = 0.3;
ta = 0.3;
% ta = 0.0684;
tav = 0.7859;
learning_rate = 0.001;
sigma_soa = 0.2;

%% initiate recalibration effect for all soas

soa_recal = cell(1, n_soas);
last_recal = NaN(n_soas, sim_trial); % summarize the last recalibration effect

%% take MCD parts out from the loop

% initiate time series for the filter
% t should be the same length as the signal time series
nsample = ts_duration * fs + 1; 
t = linspace(0, ts_duration, nsample);

% temporal constants
param=[ta, tv, tav];

% low-pass filters (Equation 1)
fa=fft(t.*exp(-t./param(1)));
fv=fft(t.*exp(-t./param(2)));
fav=fft(t.*exp(-t./param(3)));

%% run simulation

for i = 1:n_soas
    soa = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t = 1:sim_trial

        soa_recal{i}(t,:) = update_recal_gaussian(exposure_trial, soa, ...
            ts_duration, fs, stim_dura, ...
            fa, fv, fav,...
            bias, sigma_soa, learning_rate);

        last_recal(i, t) = soa_recal{i}(t,end);

    end
end

%% summarize and plot recalibration effect by adaptor SOA

last_recal_iSOA_error = std(last_recal, [], 2);
last_recal_iSOA = mean(last_recal, 2);

figure; hold on; box off
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 500, 400]);
errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2)
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('Adaptor SOA (s)')
ylabel('recalibration (s)')

%% plot the histogram of recalibration effect for each adaptor SOA

figure; hold on; box off
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 1000, 800]);
for i = 1:n_soas
    subplot(3,3,i)
    histogram(last_recal(i,:),20,'FaceColor','k','FaceAlpha',0.3,'EdgeColor','w');
    title(['' num2str(adaptor_soas(i)*1e3) ' ms'])
end
sgtitle('Recal effect distribution by adaptor SOA')


%
% % plot simulated recalibration effect (mu_shift)
% figure; hold on
% last_recal_effect = recal_effect(:, end);
% histogram(last_recal_effect)