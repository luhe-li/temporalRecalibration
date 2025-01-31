
clear; close all; clc;

%% model info

specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Atheoretical'}; % Column 2: specifications
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','exp_shiftMu'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% atheoretical model for baseline
athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');

%% load recal models

model_slc = 3;
n_model = numel(model_slc);
sub_slc = [1:4,6:10];
save_fig = 1;

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_models_VBMC', folders{model_slc(mm)});
    R(mm, :) = load_subject_data(result_folder, sub_slc, 'sub-*');
    
    for ss = 1:numel(sub_slc)
        bestP{mm, ss} = R{mm, ss}.diag.post_mean;
        pre_simul = R{mm, ss}.pred.pre_pmf(2,:);
        post_simul = R{mm, ss}.pred.post_pmf(:,2,:);
        psimul(mm, ss) = (sum(pre_simul,'all') + sum(post_simul,'all'))./(numel(pre_simul) + numel(post_simul));
        pcc(mm, ss) = bestP{mm, ss}(6);
        criterion(mm, ss) = bestP{mm, ss}(4);
    end
end

%% plot correlations 1

% Scatter plot with larger and thicker dots
figure;
scatter(pcc, criterion, 50, 'k', 'filled', 'LineWidth', 1.5);
xlim([0, 1]);
ylim([0, 300]);
xlabel('Prior of a common cause');
ylabel('Criterion of simultaneity (ms)');

% Calculate correlation coefficient and p-value
[R, P] = corrcoef(pcc, criterion);
r = R(1, 2);
p = P(1, 2);

% Add title with r and p values
title(sprintf('Correlation: r = %.2f, p = %.3f', r, p));

% Fit and plot linear regression line
hold on;
p_line = polyfit(pcc, criterion, 1);
xfit = linspace(0, 1, 100);
yfit = polyval(p_line, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 2);
hold off;

if save_fig
    saveas(gca, fullfile(out_dir, 'corr2'),'png');
end

%% plot correlation 2
% Scatter plot with larger and thicker dots
figure;
scatter(pcc, psimul, 50, 'k', 'filled', 'LineWidth', 1.5);
xlim([0, 1]);
ylim([0, 1]);
xlabel('Prior of a common cause');
ylabel('Probability of reporting simultaneity');

% Calculate correlation coefficient and p-value
[R, P] = corrcoef(pcc, psimul);
r = R(1, 2);
p = P(1, 2);

% Add title with r and p values
title(sprintf('Correlation: r = %.2f, p = %.3f', r, p));

% Fit and plot linear regression line
hold on;
p_line = polyfit(pcc, psimul, 1);
xfit = linspace(0, 1, 100);
yfit = polyval(p_line, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 2);
hold off;


if save_fig
    saveas(gca, fullfile(out_dir, 'corr1'),'png');
end
