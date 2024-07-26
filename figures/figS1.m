% fig S1: model comparison of atheoretical models of recalibration

clear; close all; clc;

%% model info

specifications = {'Exponential likelihood, criterion-shift', 'Exponential likelihood, bias-shift', 'Gaussian likelihood, criterion-shift',  'Gaussian likelihood, bias-shift'};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
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
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load different atheoretical models

n_model = numel(folders);
sub_slc = [1:4,6:10];

for mm = 1:n_model

    result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', folders{mm});
    atheo{mm} = load_subject_data(result_folder, sub_slc, 'sub-*');

    for ss = 1:numel(sub_slc)

            log_model_evi(mm, ss) = atheo{mm}{ss}.diag.bestELBO.elbo; %atheo{mm}{ss}.diag.bestELCBO;
            
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot model evidence

% max model relative to others
delta_LME = max(log_model_evi, [], 1) - log_model_evi;

figure
subplot('Position', [0.2, 0.1, 0.6, 0.8]); 
h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flip(bone),...
    'ColorLimits', [0, log(100)], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
%     'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

%% bar plot of number of winning

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
set(gca, 'YLim', [0.5 4.5]);

flnm = 'modelEvidence_atheo_models';
saveas(gca, fullfile(out_dir, flnm),'pdf')