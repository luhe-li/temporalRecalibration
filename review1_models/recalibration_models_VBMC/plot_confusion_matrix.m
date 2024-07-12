% plot confusion matrix for model comparison

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
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% atheoretical model for baseline
athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');

%% load recal models

model_slc = [1:5];%[1,2,5];
n_model = numel(model_slc);
sub_slc = [3,4];%[1:4,6:10];

for mm = 1:n_model

    if mm ~= numel(folders)
        recal_folder = fullfile(projectDir, 'recalibration_models_VBMC', folders{mm});
    else
        recal_folder = fullfile(projectDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
    end
    files = dir(fullfile(recal_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        i_sub = sub_slc(ss);
        i_data = load(fullfile(recal_folder, files(i_sub).name));
        DATA(mm, ss) = i_data;
        log_model_evi(mm, ss) = i_data.diag.bestELCBO;
        bestP{mm, ss} = i_data.diag.post_mean;
        pred{mm, ss} = i_data.pred;

        % extract summary data for plot
        pred_recal(mm, ss, :) = mean(pred{mm, ss}.pss_shift,2);

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

flnm = 'ModelEvidence_recal_models';
saveas(gca, fullfile(out_dir, flnm),'png')

%% 2. plot group average log model evidence

m_LME = mean(log_model_evi,2);
figure;
bar(exp(delta_LME));

