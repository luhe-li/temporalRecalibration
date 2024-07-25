% fig S8: table of individual parameter estimates

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

model_slc = 3;
sub_slc = [1:4,6:10];
dispLatex = 1;

result_folder = fullfile(dataDir, 'recalibration_models_VBMC', folders{model_slc});
R = load_subject_data(result_folder, sub_slc, 'sub-*');

for ss = 1:numel(sub_slc)
    bestP{ss} = R{ss}.diag.post_mean;
    est(ss,:) = bestP{ss};
end

paraID = {'\beta','\tau_A','\tau_V','Criterion','\lambda','p_{common}','\alpha','\sigma_{C=1}','\sigma_{C=2}'};

%% create latex table

est_mean                         = mean(est);
est_std                          = std(est)./sqrt(numel(sub_slc));
para.data                        = round([est; est_mean; est_std; R{ss}.model.initVal.lb;R{ss}.model.initVal.ub],3);
para.tableRowLabels              = [sprintfc('S%i', 1:numel(sub_slc)),'Mean','S.E.','Lower bound','Upper bound'];
para.tableColLabels              = paraID;

if dispLatex
    para.dataFormat = {'%.2f', 4, '%.3f',3, '%.2f',2};
    latexTable(para);
end