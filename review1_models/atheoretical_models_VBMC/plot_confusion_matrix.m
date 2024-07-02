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
sub_slc = 1:10;%[1:4,6:10];
beta_lcb = 3;

maxELCBO = zeros(n_model, numel(sub_slc));
elbos = zeros(n_model, numel(sub_slc));
bestP = cell(n_model, numel(sub_slc));

for i_model = 1:n_model

    curr_folder = fullfile(pwd, folders{i_model});
    files = dir(fullfile(curr_folder, 'sub-*'));

    for idx_sub = 1:numel(sub_slc)

        i_sub = sub_slc(idx_sub);
        i_data = load(fullfile(curr_folder, files(i_sub).name));
        DATA(i_model, i_sub) = i_data;
%         vps{i_model, i_sub, :} = i_data.model.vp;

        % use runs where exitflag is true and find the largest ELCBO among
        % them, store the posterior
        success_flag = logical(i_data.model.exitflag);
        [maxELCBO, idx_best] = max(i_data.model.elbo(success_flag) - beta_lcb * i_data.model.elbo_sd(success_flag));
        temp_vp = i_data.model.vp(success_flag);
        elbos(i_model, i_sub) = maxELCBO;
        vps{i_model, i_sub} = temp_vp(idx_best);

    end
end


%% 1. plot model evidence

% subtract min across models
deltaELCBO = max(elbos, [], 1) - elbos;

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

