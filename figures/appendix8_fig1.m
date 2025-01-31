% Appendix 8 - fig 1
% Simulation of recalibration effect by neural population code
% model by Roach, Heron, Whitaker, & McGraw (2011).

clear; clc; close all;

%% manage paths

restoredefaultpath;
out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% parameters

numNeurons = 29;               % number of neurons
tuningWidth = 220.6;           % width of tuning curves
maxGainReduction = 0.41;       % maximal proportional gain reduction
adaptationWidth = 122.6;       % breadth of the gain field for adaptation
baseGain = 1;                  % unadapted response gain

maxSOA = 1000;
exampleAdaptorSOAs = [175, 525];          % example adaptor SOAs for subplots
allAdaptorSOAs = linspace(-700, 700, 9);  % adaptor SOAs for bias calculation
stimulusSOAs = -maxSOA:1:maxSOA;          % stimulus SOAs
preferredSOAs = linspace(-500, 500, numNeurons);  % preferred SOAs of neurons
SOAHypotheses = stimulusSOAs;

%% color Palette

cmp1 =  [88,193,238;175, 213, 128; 216, 49, 91]./255;
numColors = numel(allAdaptorSOAs);
colors = zeros(numColors, 3);
for i = 1:3
    colors(:, i) = interp1([1, (numColors + 1)/2, numColors], cmp1(:, i), 1:numColors);
end

%% model

% tuning function
tuningFunction = @(SOA, preferredSOA) baseGain * exp(-(SOA - preferredSOA).^2 / (2 * tuningWidth^2));

% adaptation model
adaptationGain = @(preferredSOA, adaptorSOA) baseGain * (1 - maxGainReduction * exp(-(preferredSOA - adaptorSOA).^2 / (2 * adaptationWidth^2)));

% log likelihood calculation
logLikelihood = @(SOA, response) sum(response .* log(tuningFunction(SOA, preferredSOAs))) - sum(tuningFunction(SOA, preferredSOAs)) - sum(gammaln(response + 1));
bias = NaN(numel(allAdaptorSOAs), numel(stimulusSOAs));
PSS = zeros(1, numel(allAdaptorSOAs));

% precompute responses and log-likelihoods
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

% compute the responses for each neuron to each stimulus SOA
for aa = 1:numel(exampleAdaptorSOAs)
    adaptorSOA = exampleAdaptorSOAs(aa);
    responsesAdapted{aa} = zeros(numNeurons, length(stimulusSOAs));
    for i = 1:numNeurons
        for j = 1:length(stimulusSOAs)
            responsesAdapted{aa}(i, j) = adaptationGain(preferredSOAs(i), adaptorSOA) * tuningFunction(stimulusSOAs(j), preferredSOAs(i));
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lineWidth = 0.5;
fontSize = 7;

%% fig C-D. simulated bias and recalibration effects

figure;
set(gcf, 'Position', [0, 0, 420, 120]);

left_margin = 0.08;
right_margin = 0.02;
bottom_margin = 0.15;
top_margin = 0.1;
horizontal_gap = 0.15;
total_width = 1 - left_margin - right_margin - horizontal_gap;
total_height = 1 - bottom_margin - top_margin;
width1 = total_width * (5 / 10);
width2 = total_width * (4 / 10);
pos1 = [left_margin, bottom_margin, width1, total_height];
pos2 = [left_margin + width1 + horizontal_gap, bottom_margin, width2, total_height];

% (C) bias vs. Stimulus SOA
ax1 = axes('Position', pos1);
set(ax1, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
hold(ax1, 'on');
for aa = 1:numel(allAdaptorSOAs)
    plot(ax1, stimulusSOAs, bias(aa, :), 'Color', colors(aa, :), 'LineWidth', lineWidth*2);
end
xlim(ax1, [-900, 900]);
ylim(ax1, [-60, 60]);
xlabel(ax1, 'Stimulus SOA (ms)');
ylabel(ax1, 'Bias (ms)');

colormap(ax1, colors);
clim(ax1, [allAdaptorSOAs(1), allAdaptorSOAs(end)]);
cb = colorbar(ax1, 'Ticks', allAdaptorSOAs, 'TickLabels', arrayfun(@num2str, allAdaptorSOAs, 'UniformOutput', false));
cb.Position(1) = pos1(1) + width1 + 0.005; 
cb.Position(3) = 0.02; 
grid(ax1, 'on');

% (D) Recalibration vs. Adapted SOA
ax2 = axes('Position', pos2);
set(ax2, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
hold(ax2, 'on');
plot(ax2, allAdaptorSOAs, PSS, 'k', 'LineWidth', lineWidth, 'DisplayName', 'PSS');
ylim(ax2, [-60, 60]);
xlim(ax2, [-700, 700]);
xlabel(ax2, 'Adapter SOA (ms)');
ylabel(ax2, 'Recalibration effect (ms)');
grid(ax2, 'on');

pos2(3) = pos2(3) - 0.05;
set(ax2, 'Position', pos2);

set(gcf, 'Renderer', 'painters');
print(gcf, fullfile(out_dir, 'population_model2'), '-dpdf');
saveas(gcf, fullfile(out_dir, 'population_model2'), 'pdf');

%% Fig A-B

figure;
set(gcf, 'Position',[0,0,420,240]);

for aa = 1:2
    adaptorSOA = exampleAdaptorSOAs(aa);
    [~, colorIdx] = min(abs(allAdaptorSOAs - adaptorSOA));
    currentColor = colors(colorIdx, :);

    % left: adapted responses
    subplot(2, 3, (aa - 1) * 3 + 1);
    set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
    hold on;
    plot(stimulusSOAs, responsesAdapted{aa}, 'Color', [0.5, 0.5, 0.5]);
    xlabel('Preferred SOA (ms)');
    ylabel('Neural response');
    xlim([-700, 700]);
    xticks([-500, 0, 500]);
    xline(adaptorSOA, 'k-', 'LineWidth', lineWidth * 2); % Mark the adapted SOA

    % middle
    subplot(2, 3, (aa - 1) * 3 + 2);
    set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
    hold on;

    % responses to 0 ms stimulus before and after adaptation
    responsesZeroUnadapted = tuningFunction(0, preferredSOAs);
    responsesZeroAdapted = adaptationGain(preferredSOAs, adaptorSOA) .* tuningFunction(0, preferredSOAs);

    % estimated SOA before adaptation
    logLikelihoodUnadapted = arrayfun(@(SOA) logLikelihood(SOA, responsesZeroUnadapted), SOAHypotheses);
    [~, maxIdxUnadapted] = max(logLikelihoodUnadapted);
    estimatedSOAUnadapted = SOAHypotheses(maxIdxUnadapted);

    % estimated SOA after adaptation
    logLikelihoodAdapted = arrayfun(@(SOA) logLikelihood(SOA, responsesZeroAdapted), SOAHypotheses);
    [~, maxIdxAdapted] = max(logLikelihoodAdapted);
    estimatedSOAAdapted = SOAHypotheses(maxIdxAdapted);

    plot(preferredSOAs, responsesZeroUnadapted, 'k', 'LineWidth', lineWidth, 'DisplayName', 'Unadapted Response');
    plot(preferredSOAs, responsesZeroAdapted, 'k:', 'LineWidth', lineWidth, 'DisplayName', 'Adapted Response');

    % mark the physical stimulus
    xline(0, '-', 'Color',[90,174,97]./255,'LineWidth', lineWidth*2, 'DisplayName', 'Physical Stimulus');

    % mark the estimated SOA after adaptation
    xline(estimatedSOAAdapted, '--','Color',[153,112,171]./255, 'LineWidth', lineWidth*2, 'DisplayName', 'Estimated SOA (adapted)');

    xlabel('Preferred SOA (ms)');
    ylabel('Neural response');
    xlim([-500, 500]);
    xticks([-500, 0, 500]);
    legend('Location', 'Best');

    % right: bias plot with matching color
    subplot(2, 3, (aa - 1) * 3 + 3);
    set(gca, 'LineWidth', lineWidth, 'FontSize', fontSize, 'TickDir', 'out');
    hold on;
    idx = find(allAdaptorSOAs == adaptorSOA);
    plot(stimulusSOAs, bias(idx, :), 'Color', currentColor, 'LineWidth', lineWidth*2);
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