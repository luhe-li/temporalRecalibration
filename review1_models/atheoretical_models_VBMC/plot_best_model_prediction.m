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
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% find best-fitting model for each subject

n_model = numel(folders);
sub_slc = 1:8;%1:10;
beta_lcb = 3;

for i_model = 1:n_model

    curr_folder = fullfile(pwd, folders{i_model});
    files = dir(fullfile(curr_folder, 'sub-*'));

    for idx_sub = 1:numel(sub_slc)

        i_sub = sub_slc(idx_sub);
        i_data = load(fullfile(curr_folder, files(i_sub).name));
        DATA(i_model, i_sub) = i_data;
        vps{i_model, i_sub, :} = i_data.model.vp;

        % use runs where exitflag is true and find the largest ELCBO
        success_flag = logical(i_data.model.exitflag);
        [maxELCBO, idx_best] = max(i_data.model.elbo(success_flag) - beta_lcb * i_data.model.elbo_sd(success_flag));
        temp_vp = i_data.model.vp(success_flag);
        elbos(i_model, i_sub) = maxELCBO;
        best_vp{i_model, i_sub} = temp_vp(idx_best);

    end
end

% extract parameters for the best model
deltaELCBO = max(elbos, [], 1) - elbos;



