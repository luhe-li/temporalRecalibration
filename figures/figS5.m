% fig S5: individual results of recalibration model comparison

clear; close all; clc;

%% model info

specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Fixed update, asymmetric', 'Fixed update, symmetric','Atheoretical'}; % Column 2: specifications
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym','exp_shiftMu'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), 'Google Drive','My Drive','temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load recal models

model_slc = [6,5,2,1,4,3];
n_model = numel(model_slc);
sub_slc = [1:4,6:10];

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_models_VBMC', folders{mm});
    R(mm, :) = load_subject_data(result_folder, sub_slc, 'sub-*');
    
    for ss = 1:numel(sub_slc)
        log_model_evi(mm, ss) = R{mm, ss}.diag.bestELCBO;
    end
end

%% load atheoretical model

result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

for ss = 1:numel(sub_slc)
   log_model_evi(n_model+1, ss)= atheo{ss}.diag.bestELCBO;
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. plot model evidence with atheoretical model

% max subtract other log model evidence
delta_LME = max(log_model_evi, [], 1) - log_model_evi; 

figure
set(gcf, 'Position',[0 0 400 170])
set(gca, 'FontSize', 8)

h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, log(100)], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
h.YDisplayLabels = specifications([model_slc,7]);
h.XDisplayLabels = num2cell(1:numel(sub_slc));

flnm = 'A_BF_all_models';
saveas(gca, fullfile(out_dir, flnm),'pdf')

%% 2. plot model evidence within recalibration model

% max subtract other log model evidence
delta_LME = max(log_model_evi(model_slc,:), [], 1) - log_model_evi(model_slc,:); 

figure
set(gcf, 'Position',[0 0 400 150])

subplot('Position', [0.2, 0.1, 0.6, 0.8]); 
set(gca, 'FontSize', 8)
h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, log(100)], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;

h.YDisplayLabels = specifications(model_slc);
h.XDisplayLabels = num2cell(1:numel(sub_slc));


n_win = sum(delta_LME<log(10),2);
subplot('Position', [0.9, 0.1, 0.1, 0.8]);

b = barh(1:length(n_win), n_win, 'FaceColor', [0.8, 0.8, 0.8]);
set(gca, 'YDir', 'reverse', 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Remove axis and background

yticks(1:length(n_win));
yticklabels({' ', ' ', ' ', ' '});

for i = 1:length(n_win)
    text(n_win(i) + 0.5, i, num2str(n_win(i)), 'VerticalAlignment', 'middle');
end

set(gca, 'XLim', [0 max(n_win) + 2]);
set(gca, 'YLim', [0.5 6.5]);

flnm = 'B_BF_recal_models';
saveas(gca, fullfile(out_dir, flnm),'pdf')