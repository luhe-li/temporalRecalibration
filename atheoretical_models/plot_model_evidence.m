clear; close all; clc;

%% Model Info

specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion',  'Gaussian likelihood, shift bias'};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% Manage Paths

restoredefaultpath;
current_dir = pwd;
[project_dir, ~] = fileparts(current_dir);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(git_dir, 'vbmc')));
out_dir = fullfile(current_dir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Load Data

n_model = numel(folders);
sub_slc = [1:4, 6:10];

for mm = 1:n_model

    result_folder = fullfile(project_dir, 'fit_results','atheoretical_models', folders{mm});
    files = dir(fullfile(result_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        i_sub = sub_slc(ss);
        i_data = load(fullfile(result_folder, files(i_sub).name));
        DATA(mm, ss) = i_data;
        log_model_evi(mm, ss) = i_data.diag.bestELCBO;
        bestP{mm, ss} = i_data.diag.post_mean;
        pred{mm, ss} = i_data.pred;

    end
end

%% Plot Model Evidence

% Max subtract other log model evidence
delta_LME = max(log_model_evi, [], 1) - log_model_evi; 

figure
h = heatmap(round(delta_LME, 1), 'XLabel', 'Participant', ...
    'Colormap', flipud(bone), ...
    'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off', ...
    'FontSize', 8);
colorbar;

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% Save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position', [0 0 400 110])

flnm = 'model_evidence_atheo_models';
saveas(gca, fullfile(out_dir, flnm), 'png')