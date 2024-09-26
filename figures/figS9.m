% fig S9: simulation of recalibration effect by neural population code
% model by Roach, Heron, Whitaker, & McGraw (2011).

clear; clc; close all;

%% manage paths

restoredefaultpath;
out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Parameters

numNeurons = 29;               % Number of neurons
tuningWidth = 220.6;           % Width of tuning curves
maxGainReduction = 0.41;       % Maximal proportional gain reduction
adaptationWidth = 122.6;       % Breadth of the gain field for adaptation
baseGain = 1;                  % Unadapted response gain

maxSOA = 1000;
exampleAdaptorSOAs = [175, 700];          % Example adaptor SOAs for subplots
allAdaptorSOAs = linspace(-700, 700, 9);  % Adaptor SOAs for bias calculation
stimulusSOAs = -maxSOA:1:maxSOA;          % Stimulus SOAs
preferredSOAs = linspace(-500, 500, numNeurons);  % Preferred SOAs of neurons
SOAHypotheses = stimulusSOAs;

%% Color Palette

% Create a diverging color palette
cmp1 =  [88,193,238;175, 213, 128; 216, 49, 91]./255;
numColors = numel(allAdaptorSOAs);
% Interpolate to create a diverging palette
colors = zeros(numColors, 3);
for i = 1:3
    colors(:, i) = interp1([1, (numColors + 1)/2, numColors], cmp1(:, i), 1:numColors);
end

%% Model

% Tuning function
tuningFunction = @(SOA, preferredSOA) baseGain * exp(-(SOA - preferredSOA).^2 / (2 * tuningWidth^2));

% Adaptation model
adaptationGain = @(preferredSOA, adaptorSOA) baseGain * (1 - maxGainReduction * exp(-(preferredSOA - adaptorSOA).^2 / (2 * adaptationWidth^2)));

% Log likelihood calculation
logLikelihood = @(SOA, response) sum(response .* log(tuningFunction(SOA, preferredSOAs))) - sum(tuningFunction(SOA, preferredSOAs)) - sum(gammaln(response + 1));

bias = NaN(numel(allAdaptorSOAs), numel(stimulusSOAs));
PSS = zeros(1, numel(allAdaptorSOAs));

% Precompute responses and log-likelihoods
for aa = 1:numel(allAdaptorSOAs)
    adaptorSOA = allAdaptorSOAs(aa);
    estimatedSOA = zeros(1, length(stimulusSOAs));

    for j = 1:length(stimulusSOAs)
        SOA = stimulusSOAs(j);
        response = zeros(1, numNeurons);
        for i = 1:numNeurons
            response(i) = tuningFunction(SOA, preferredSOAs(i)) * adaptationGain(preferredSOAs(i), adaptorSOA);
        end
        logLikelihoods = arrayfun(@(hypSOA) logLikelihood(hypSOA, response), SOAHypotheses);
        [~, maxIdx] = max(logLikelihoods);
        estimatedSOA(j) = SOAHypotheses(maxIdx);
    end

    [~, idx] = min(abs(estimatedSOA));
    PSS(aa) = stimulusSOAs(idx);
    bias(aa, :) = estimatedSOA - stimulusSOAs;
end

% Compute the responses for each neuron to each stimulus SOA
for aa = 1:numel(exampleAdaptorSOAs)
    adaptorSOA = exampleAdaptorSOAs(aa);
    responsesAdapted{aa} = zeros(numNeurons, length(stimulusSOAs));
    for i = 1:numNeurons
        for j = 1:length(stimulusSOAs)
            responsesAdapted{aa}(i, j) = adaptationGain(preferredSOAs(i), adaptorSOA) * tuningFunction(stimulusSOAs(j), preferredSOAs(i));
        end
    end
end

%% Plotting

% Figure setup
lineWidth = 1;
fontSize = 10;

%% First Figure (1x2 subplots)
figure;
set(gcf, 'Position',[0,0,420,120]);

% (1) Bias vs. Stimulus SOA with color palette
subplot(1, 2, 1);
set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
hold on;
for aa = 1:numel(allAdaptorSOAs)
    plot(stimulusSOAs, bias(aa, :), 'Color', colors(aa, :), 'LineWidth', lineWidth);
end
xlim([-900, 900]);
ylim([-60, 60]);
xlabel('Stimulus SOA (ms)');
ylabel('Bias (ms)');
legend
% % Create a colorbar with labels corresponding to adaptor SOAs
% colormap(colors);
% caxis([allAdaptorSOAs(1), allAdaptorSOAs(end)]);
% colorbar('Ticks', allAdaptorSOAs, 'TickLabels', arrayfun(@num2str, allAdaptorSOAs, 'UniformOutput', false));
% grid on;

% (2) Mean Bias vs. Adapted SOA
subplot(1, 2, 2);
set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
hold on;
plot(allAdaptorSOAs, PSS, 'k', 'LineWidth', lineWidth);
ylim([-60, 60]);
xlim([-700, 700]);
xlabel('Adaptor SOA (ms)');
ylabel('Recalibration effect (ms)');
grid on;

set(gcf, 'Renderer', 'painters');
print(gcf, fullfile(out_dir, 'population_model2'), '-dpdf')
saveas(gcf, fullfile(out_dir, 'population_model2'), 'pdf')

%% Second Figure (2x3 subplots)
figure;
set(gcf, 'Position',[0,0,420,240]);

for aa = 1:2
    adaptorSOA = exampleAdaptorSOAs(aa);

    % Find the index of the adaptor SOA in allAdaptorSOAs to match colors
    [~, colorIdx] = min(abs(allAdaptorSOAs - adaptorSOA));
    currentColor = colors(colorIdx, :);

    % (1) Adapted responses
    subplot(2, 3, (aa - 1) * 3 + 1);
    set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
    hold on;
    plot(stimulusSOAs, responsesAdapted{aa}, 'Color', [0.5, 0.5, 0.5]);
    xlabel('Preferred SOA (ms)');
    ylabel('Neural response');
    xlim([-maxSOA, maxSOA]);
    xticks([-500, 0, 500]);
    xline(adaptorSOA, 'k-', 'LineWidth', lineWidth * 2); % Mark the adapted SOA

    % (2) Middle panel
    subplot(2, 3, (aa - 1) * 3 + 2);
    set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
    hold on;

    % Compute responses to 0 ms stimulus before and after adaptation
    responsesZeroUnadapted = tuningFunction(0, preferredSOAs);
    responsesZeroAdapted = adaptationGain(preferredSOAs, adaptorSOA) .* tuningFunction(0, preferredSOAs);

    % Compute estimated SOA before adaptation
    logLikelihoodUnadapted = arrayfun(@(SOA) logLikelihood(SOA, responsesZeroUnadapted), SOAHypotheses);
    [~, maxIdxUnadapted] = max(logLikelihoodUnadapted);
    estimatedSOAUnadapted = SOAHypotheses(maxIdxUnadapted);

    % Compute estimated SOA after adaptation
    logLikelihoodAdapted = arrayfun(@(SOA) logLikelihood(SOA, responsesZeroAdapted), SOAHypotheses);
    [~, maxIdxAdapted] = max(logLikelihoodAdapted);
    estimatedSOAAdapted = SOAHypotheses(maxIdxAdapted);

    % Plot responses
    plot(preferredSOAs, responsesZeroUnadapted, 'k', 'LineWidth', lineWidth, 'DisplayName', 'Unadapted Response');
    plot(preferredSOAs, responsesZeroAdapted, 'k:', 'LineWidth', lineWidth, 'DisplayName', 'Adapted Response');

    % Mark the physical stimulus
    xline(0, '--', 'Color',[90,174,97]./255,'LineWidth', lineWidth*2, 'DisplayName', 'Physical Stimulus');

%     % Mark the estimated SOA before adaptation
%     xline(estimatedSOAUnadapted, 'k:', 'LineWidth', lineWidth, 'DisplayName', 'Estimated SOA (unadapted)');

    % Mark the estimated SOA after adaptation
    xline(estimatedSOAAdapted, '--','Color',[153,112,171]./255, 'LineWidth', lineWidth*2, 'DisplayName', 'Estimated SOA (adapted)');

    % Mark the adaptor SOA
    xline(adaptorSOA, 'k-', 'LineWidth', lineWidth*2, 'DisplayName', 'Adaptor SOA');

    xlabel('Preferred SOA (ms)');
    ylabel('Neural response');
    xlim([-500, 500]);
    xticks([-500, 0, 500]);
    legend('Location', 'Best');

    % (3) Bias plot with matching color
    subplot(2, 3, (aa - 1) * 3 + 3);
    set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
    hold on;
    idx = find(allAdaptorSOAs == adaptorSOA);
    plot(stimulusSOAs, bias(idx, :), 'Color', currentColor, 'LineWidth', lineWidth);
    yline(0, 'k--');
    xlabel('Stimulus SOA (ms)');
    ylabel('Bias (ms)');
    ylim([-60, 60]);
    xticks([-500, 0, 500]);
    xlim([-900, 900]);
    grid on;
    title(sprintf('Bias (Adaptor SOA = %d ms)', adaptorSOA));
end

set(gcf, 'Renderer', 'painters');
print(gcf, fullfile(out_dir, 'population_model1'), '-dpdf')
saveas(gcf, fullfile(out_dir, 'population_model1'), 'pdf')