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
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load

n_model = numel(folders);
n_sub = 10;

for i_model = 1:n_model

    curr_folder = fullfile(pwd, folders{i_model});
    
    files = dir(fullfile(curr_folder, 'sub-*'));

    for i_sub = 1:n_sub

        i_data = load(fullfile(curr_folder, files(i_sub).name));
        DATA(i_model, i_sub) = i_data;
        NLL(i_model, i_sub) = i_data.model.minNLL;
        AIC(i_model, i_sub) = 2*i_data.model.minNLL + 2*i_data.model.initVal.num_para;
        BIC(i_model, i_sub) = 2*i_data.model.minNLL ...
            + DATA(i_sub).model.initVal.num_para ...
            * log(numel(i_data.data(1).pre_s_unique) * i_data.data(1).pre_numTrials * 2 * numel(i_data.data));
            % SOA conditions x trials per condition x pre/post-test x
            % sessions

    end

end

%% 1. plot AIC

sub_slc = [1:4, 6:10];

% subtract min across models
AIC = AIC(:,sub_slc);
deltaAIC = AIC - min(AIC, [], 1);

figure
h = heatmap(round(deltaAIC, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8); 
colorbar;

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));
% 
% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'AIC_atheo_models';
saveas(gca, fullfile(out_dir, flnm),'pdf')

%% 2. plot BIC

% subtract min across models
BIC = BIC(:,sub_slc);
deltaBIC = BIC - min(BIC, [], 1);

figure
h = heatmap(round(deltaBIC, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8); 
colorbar;

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'BIC_atheo_models';
saveas(gca, fullfile(out_dir, flnm),'pdf')
