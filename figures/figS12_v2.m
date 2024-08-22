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

x_max = 1000;
adapted_SOA = [175, 700]; % example adaptor soa for subplots 1-3, 4-6
adapted_SOAs = linspace(-700, 700, 9); %  adaptor soa for subplots 7
stimulus_SOA = -x_max:1:x_max; % Stimulus SOAs
preferred_SOAs = linspace(-500, 500, N); % preferred SOAs of neurons (in ms)
SOA_hyp = stimulus_SOA;

%% Model

% Tuning function
fi = @(SOA, SOAi) G0 * exp(-(SOA - SOAi).^2 / (2 * sigma^2));

% Adaptation model
Gi = @(SOAi, SOA_a) G0 * (1 - alpha * exp(-(SOAi - SOA_a).^2 / (2 * sigma_a^2)));

% Log likelihood calculation
logL = @(SOA, R) sum(R .* log(fi(SOA, preferred_SOAs))) - sum(fi(SOA, preferred_SOAs)) - sum(gammaln(R + 1));

bias = NaN(numel(adapted_SOAs), numel(stimulus_SOA));

% Precompute responses and log-likelihoods
for aa = 1:numel(adapted_SOAs)
    SOA_a = adapted_SOAs(aa);
    estimated_SOA = zeros(1, length(stimulus_SOA));

    for j = 1:length(stimulus_SOA)
        SOA = stimulus_SOA(j);
        R = arrayfun(@(i) fi(SOA, preferred_SOAs(i)) * Gi(preferred_SOAs(i), SOA_a), 1:N);
        log_likelihoods = arrayfun(@(SOA) logL(SOA, R), SOA_hyp, 'UniformOutput', true);
        [~, max_idx] = max(log_likelihoods);
        estimated_SOA(j) = SOA_hyp(max_idx);
    end

    [~, idx] = min(abs(estimated_SOA));
    PSS(aa) = stimulus_SOA(idx);
    bias(aa, :) = estimated_SOA - stimulus_SOA;
 
end

% Compute the responses for each neuron to each stimulus SOA
for aa = 1:numel(adapted_SOA)
    i_adapted_SOA = adapted_SOA(aa);
    for i = 1:N
        for j = 1:length(stimulus_SOA)
            responses_unadapted{aa}(i, j) = fi(stimulus_SOA(j), preferred_SOAs(i));
            responses_adapted{aa}(i, j) = Gi(preferred_SOAs(i), i_adapted_SOA) * fi(stimulus_SOA(j), preferred_SOAs(i));
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure set up

lw = 0.5;
fontSZ = 8;
titleSZ = 9;
dotSZ = 10;

%%

figure;
set(gcf, 'Position',[0,0,420,120]);

subplot(1,2,1)
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
hold on
colororder(parula(numel(adapted_SOAs)));
plot(stimulus_SOA, bias, 'LineWidth', 1);
xlim([min(preferred_SOAs), max(preferred_SOAs)])
xlabel('Stimulus SOA (ms)');
ylabel('Bias (ms)');
% legendEntries = arrayfun(@(x) sprintf('%d', x), adapted_soas, 'UniformOutput', false);
% legend(legendEntries, 'Location', 'Best');
% lgd = legend(legendEntries, 'Location', 'Best');
% title(lgd, 'Adapter SOA (ms)');
xlim([-900 900])
ylim([-60, 60])

% Plot Mean Bias vs. Adapted SOA
subplot(1,2,2);
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
hold on
plot(adapted_SOAs, PSS, 'k','LineWidth', 1);
ylim([-60, 60])
xlim([-700 700])
xlabel('Adapter SOA (ms)');
ylabel('Recalibration effect (ms)');
grid on;

saveas(gcf, fullfile(out_dir, 'population_model2'), 'pdf')

%%

figure;
set(gcf, 'Position',[0,0,420,240]);

for aa = 1:2
    % (a) Adapted responses
    subplot(2,3,(aa-1)*3+1);
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    hold on
    plot(stimulus_SOA, responses_adapted{aa},'k');
    xlabel('Stimulus SOA (ms)');
    ylabel('Neural response');
    xlim([-x_max, x_max]);
    xticks([-500, 500])
    xline(adapted_SOA(aa),'-','LineWidth',lw*2) % Mark the adapted SOA

    subplot(2,3,(aa-1)*3+2)
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    hold on
    plot(stimulus_SOA, sum(responses_unadapted{aa}, 1), 'k', 'LineWidth', 1);
    plot(stimulus_SOA, sum(responses_adapted{aa}, 1), 'r', 'LineWidth', 1);
    xline(adapted_SOA(aa),'-','LineWidth',lw*2)
    xlabel('Stimulus SOA (ms)');
    ylabel('Population response');
    xlim([-x_max, x_max]);
    xticks([-500, 500])

    subplot(2,3,(aa-1)*3+3)
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    hold on
    idx = find(adapted_SOAs == adapted_SOA(aa));
    plot(stimulus_SOA, bias(idx,:), 'k','LineWidth', 1);
    yline(0)
    xlabel('Stimulus SOA (ms)');
    ylabel('Bias (ms)');
    ylim([-60, 60])
    xticks([-500, 500])
    xlim([-900 900])

end

saveas(gcf, fullfile(out_dir, 'population_model1'), 'pdf')