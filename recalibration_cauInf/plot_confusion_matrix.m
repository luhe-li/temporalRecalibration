% plot confusion matrix for model comparison

clear; close all; clc;

%% model info
specifications = {'Full model','No \sigma_{C1}', 'No \sigma_{C2}','No \sigma_{C1}, \sigma_{C2}'}; % Column 2: specifications
folders = {'cauInf_asym', 'cauInf_asym_xSigmaC1', 'cauInf_asym_xSigmaC2', 'cauInf_asym_xSigmaC1C2'}; % Column 3: folder names
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

%% load recal models

model_slc = [1,2,4];
n_model = numel(model_slc);
sub_slc = [1:4,6:10];

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_cauInf', folders{model_slc(mm)});
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
h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 2.3], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
%     'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'ModelEvidence_all_models';
saveas(gca, fullfile(out_dir, flnm),'png')

%% 2. plot model evidence within recalibration model

% max subtract other log model evidence
delta_LME = max(log_model_evi(model_slc,:), [], 1) - log_model_evi(model_slc,:); 

figure
h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 2.3], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
%     'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...

h.YDisplayLabels = specifications(model_slc);
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'ModelEvidence_recal_models';
saveas(gca, fullfile(out_dir, flnm),'png')

%% 3. plot group log bayes factor
order = [6, 5,2, 1, 4, 3];
delta = log_model_evi(order, :) - log_model_evi(6, :);
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
ylabel({'\Delta log model evidence'; 'relative to fixed-update symmetric model'});

flnm = 'group_log_model_evidence';
saveas(gca, fullfile(out_dir, flnm),'png')
