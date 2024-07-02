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
beta_lcb = 3;

maxELCBO = zeros(n_model, numel(sub_slc));
elbos = zeros(n_model, numel(sub_slc));
bestP = cell(n_model, numel(sub_slc));

diag_flnm = sprintf('VBMC_diag_results.mat');
if ~exist(fullfile(out_dir, diag_flnm),'file')

    for mm = 1:n_model

        curr_folder = fullfile(pwd, folders{mm});
        files = dir(fullfile(curr_folder, 'sub-*'));

        for ss = 1:numel(sub_slc)

            fprintf('VBMC diagnostics for model-%i, sub-%i \n',mm, ss);
            i_sub = sub_slc(ss);
            i_data = load(fullfile(curr_folder, files(i_sub).name));
            DATA(mm, i_sub) = i_data;
            [exitflag(mm, ss),bestELCBO(mm, ss),idx_best,stats{mm, ss}] = vbmc_diagnostics(i_data.model.vp);

        end
    end

    save(fullfile(out_dir, diag_flnm), 'exitflag', 'bestELCBO', 'idx_best', 'stats', 'DATA');

else
end


%% 1. plot model evidence

% subtract min across models
deltaELCBO = max(bestELCBO, [], 1) - bestELCBO;

figure
h = heatmap(round(deltaELCBO, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
%     'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));
%
% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'AIC_atheo_models';
saveas(gca, fullfile(out_dir, flnm),'pdf')

%% 2. check each model's posterior agaisnt prior to validate prior selection
% if posterior is too narrow, narrown the prior and fit again

