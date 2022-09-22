% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase. The current version assumes
% that measurement_soa shifts on every exposure trial after encountering an
% SOA, whhile taking MCD_corr into account.

clear all; clc; close all;

simTrial = 10;
for t = 1:simTrial

    %% fixed parameters
    % set session parameters
    soa = 0.3; % in s, adaptor_soa, fixed in session
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
    learning_rate = 0.001;

    % For each stimulus with an fixed soa (soa = t_a - t_v), participants formed a
    % noisy sensory measurement m_soa with a bias in the head. We assume that
    % m_soa is sampled from N(soa + bias, sigma).
    bias = 0.06; % in sec
    sigma_soa = 0.05; % in sec

    %%
    recal_effect(t,:) = update_recal(exposure_trial, soa, bias, sigma_soa,...
        duration, fs, onset, stim_dura,...
        ta, tv, tav,...
        learning_rate);

end

% plot simulated recalibration effect (mu_shift)
figure; hold on
last_recal_effect = recal_effect(:, end);
histogram(last_recal_effect)