% fig S12: simulation of recalibration effect by neural population code
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
adapted_soas = [-700, -300:100:300, 700];

x_max = 1000;
adapted_SOA = -100; % example adaptor soa for subplot 1&2
stimulus_SOA = -x_max:50:x_max; % Stimulus SOAs
preferred_SOAs = linspace(-500, 500, N); % preferred SOAs of neurons (in ms)\
SOA_hyp = -1000:1000;

%% Model

% Preferred SOAs of neurons
SOAi = linspace(-500, 500, N);

% Tuning function
fi = @(SOA, SOAi) G0 * exp(-(SOA - SOAi).^2 / (2 * sigma^2));

% Adaptation model
Gi = @(SOAi, SOA_a) G0 * (1 - alpha * exp(-(SOAi - SOA_a).^2 / (2 * sigma_a^2)));

% Log likelihood calculation
logL = @(SOA, R) sum(R .* log(fi(SOA, SOAi))) - sum(fi(SOA, SOAi)) - sum(gammaln(R + 1));

bias = NaN(numel(adapted_soas), numel(stimulus_SOA));

% Precompute responses and log-likelihoods
for aa = 1:numel(adapted_soas)
    SOA_a = adapted_soas(aa);
    estimated_SOA = zeros(1, length(stimulus_SOA));

    for j = 1:length(stimulus_SOA)
        SOA = stimulus_SOA(j);
        R = arrayfun(@(i) fi(SOA, SOAi(i)) * Gi(SOAi(i), SOA_a), 1:N);
        log_likelihoods = arrayfun(@(SOA) logL(SOA, R), SOA_hyp, 'UniformOutput', true);
        [~, max_idx] = max(log_likelihoods);
        estimated_SOA(j) = SOA_hyp(max_idx);
    end

    [~, idx] = min(abs(estimated_SOA));
    PSS(aa) = stimulus_SOA(idx);
    bias(aa, :) = estimated_SOA - stimulus_SOA;
end

% Initialize response matrices
responses_unadapted = zeros(N, length(stimulus_SOA));
responses_adapted = zeros(N, length(stimulus_SOA));

% Compute the responses for each neuron to each stimulus SOA
for i = 1:N
    for j = 1:length(stimulus_SOA)
        responses_unadapted(i, j) = fi(stimulus_SOA(j), preferred_SOAs(i));
        responses_adapted(i, j) = Gi(preferred_SOAs(i), adapted_SOA) * fi(stimulus_SOA(j), preferred_SOAs(i));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 8;
titleSZ = 9;
dotSZ = 10;

figure;
set(gcf, 'Position',[0,0,420,300]);

% (a) Adapted responses
subplot(2,2,1);
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica');
hold on
plot(stimulus_SOA, responses_adapted,'k');
xlabel('Stimulus SOA (ms)');
ylabel('Neural response');
xlim([-x_max, x_max]);
hold on;
xline(adapted_SOA,'--') % Mark the adapted SOA

subplot(2,2,2)
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica');
hold on
plot(stimulus_SOA, sum(responses_unadapted, 1), 'k', 'LineWidth', 2);
plot(stimulus_SOA, sum(responses_adapted, 1), 'r', 'LineWidth', 2);
xline(adapted_SOA,'--') % Mark the adapted SOA
xlabel('Stimulus SOA (ms)');
ylabel('Population response');
xlim([-x_max, x_max]);
legend('Unadapted', 'Adapted','Adapter SOA','Location','best');
hold on;

subplot(2,2,3); 
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica');
hold on
colororder(parula(length(adapted_soas)));
plot(stimulus_SOA, bias, 'LineWidth', 1);
xlim([min(preferred_SOAs), max(preferred_SOAs)])
xlabel('Stimulus SOA (ms)');
ylabel('Bias (estimated SOA - stimulus SOA; ms)');
legendEntries = arrayfun(@(x) sprintf('%d', x), adapted_soas, 'UniformOutput', false);
legend(legendEntries, 'Location', 'Best');
lgd = legend(legendEntries, 'Location', 'Best');
title(lgd, 'Adapter SOA (ms)');
hold off;

% Plot Mean Bias vs. Adapted SOA
subplot(2,2,4); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica');
hold on
plot(adapted_soas, PSS, 'k','LineWidth', 1);
ylim([-60, 60])
xlabel('Adapter SOA (ms)');
ylabel('Mean bias across stimulus SOA (ms)');
grid on;
hold off;
saveas(gca, fullfile(out_dir, 'bias'),'pdf')