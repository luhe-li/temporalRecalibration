% fig S11: simulation of recalibration effect by neural population code
% model by Roach, Heron, Whitaker, & McGraw (2011).

clear; clc; close all;

%% manage paths

restoredefaultpath;
out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Parameters

N = 29; % Number of neurons
sigma = 220.6; % Width of tuning curves
alpha = 0.41; % Maximal proportional gain reduction
sigma_a = 122.6; % Breadth of the gain field for adaptation
G0 = 1; % Unadapted response gain
adapted_soa = [-700, -300:100:300, 700];
physical_SOA = -500:50:500; % Physical SOAs

%% Model

% Preferred SOAs of neurons
SOAi = linspace(-500, 500, N);

% Tuning function
fi = @(SOA, SOAi) G0 * exp(-(SOA - SOAi).^2 / (2 * sigma^2));

% Adaptation model
Gi = @(SOAi, SOA_a) G0 * (1 - alpha * exp(-(SOAi - SOA_a).^2 / (2 * sigma_a^2)));

% Log likelihood calculation
logL = @(SOA, R) sum(R .* log(fi(SOA, SOAi))) - sum(fi(SOA, SOAi)) - sum(gammaln(R + 1));

bias = NaN(numel(adapted_soa), numel(physical_SOA));

% Precompute responses and log-likelihoods
for aa = 1:numel(adapted_soa)
    SOA_a = adapted_soa(aa);
    estimated_SOA = zeros(1, length(physical_SOA));

    for j = 1:length(physical_SOA)
        SOA = physical_SOA(j);
        R = arrayfun(@(i) fi(SOA, SOAi(i)) * Gi(SOAi(i), SOA_a), 1:N);
        SOA_range = -500:1:500;
        log_likelihoods = arrayfun(@(SOA) logL(SOA, R), SOA_range, 'UniformOutput', true);
        [~, max_idx] = max(log_likelihoods);
        estimated_SOA(j) = SOA_range(max_idx);
    end

    bias(aa, :) = estimated_SOA - physical_SOA;
end

%% plot

% Plot Bias vs. Physical SOA
figure;
hold on;
colororder(parula(length(adapted_soa)));
plot(physical_SOA, bias, 'LineWidth', 1);
xlabel('Physical SOA (ms)');
ylabel('Bias (Estimated SOA - Physical SOA; ms)');
legendEntries = arrayfun(@(x) sprintf('%d', x), adapted_soa, 'UniformOutput', false);
legend(legendEntries, 'Location', 'Best');
lgd = legend(legendEntries, 'Location', 'Best');
title(lgd, 'Adapted SOA (ms)');
hold off;
saveas(gca, fullfile(out_dir, 'bias1'),'pdf')

% Plot Mean Bias vs. Adapted SOA
figure;
hold on;
plot(adapted_soa, mean(bias, 2), 'LineWidth', 1);
xlabel('Adapted SOA (ms)');
ylabel('Mean Bias across physical SOA');
grid on;
hold off;
saveas(gca, fullfile(out_dir, 'bias2'),'pdf')