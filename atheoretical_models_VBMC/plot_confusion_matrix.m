clear; close all; clc;

%% model info

specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion',  'Gaussian likelihood, shift bias',};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load

n_model = numel(folders);
sub_slc = [1:4,6:10];

for mm = 1:n_model

    results_folder = fullfile(projectDir, 'atheoretical_models_VBMC', folders{mm});
    files = dir(fullfile(results_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        i_sub = sub_slc(ss);
        i_data = load(fullfile(results_folder, files(i_sub).name));
        DATA(mm, ss) = i_data;
        log_model_evi(mm, ss) = i_data.diag.bestELCBO;
        bestP{mm, ss} = i_data.diag.post_mean;
        pred{mm, ss} = i_data.pred;

    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. plot model evidence

% max subtract other log model evidence
delta_LME = max(log_model_evi, [], 1) - log_model_evi; 

figure
h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
%     'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'modelEvidence_atheo_models';
saveas(gca, fullfile(out_dir, flnm),'png')

%% 2. plot group log bayes factor

% max subtract other log model evidence
delta = log_model_evi - log_model_evi(4,:); 
m_delta = mean(delta, 2);
se_delta = std(delta, [], 2) ./ numel(sub_slc);

figure;
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 420 250])
hold on;
bar_handle = bar(m_delta, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none','BarWidth', 0.6);
errorbar(m_delta, se_delta, 'k', 'LineStyle', 'none', 'CapSize', 0);
yticks([0:10:40])
ylim([0, 45])

xticks(1:length(m_delta));
labels = specifications(order);
labels = cellfun(@(x) strrep(x,',','\newline'), labels,'UniformOutput',false);
xticklabels(labels);
ylabel({'\Delta log model evidence'; 'relative to heuristic-asymmetric model'});

flnm = 'group_log_model_evidence';
saveas(gca, fullfile(out_dir, flnm),'png')
