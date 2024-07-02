% function fit_recal_model_VBMC(i_model, useCluster)

clear; close all;
i_model = 3;

%% select models

rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric'}; % Column 2: specifications
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
currModelStr = model_info.FolderName{i_model};

%% set environment

% exclude outlier
sub_slc = [1:4, 6:10];
sub = 1;

% if ~exist('useCluster', 'var') || isempty(useCluster)
%     useCluster= false;
% end
%
% % sbatch = model, job = sub, core = run
% switch useCluster
%     case true
%         if ~exist('numCores', 'var') || isempty(numCores)
%             numCores  = maxNumCompThreads;
%         end
%         hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
%         if isnan(hpc_job_number), error('Problem with array assigment'); end
%         fprintf('Job number: %i \n', hpc_job_number);
%         sub  = sub_slc(hpc_job_number);
%     case false
%         numCores = feature('numcores');
%         sub = 1;
% %         sub = str2num(input('Enter Subject ID to analyze: ', 's'));
% end
%
% % make sure Matlab does not exceed this
% fprintf('Number of cores: %i  \n', numCores);
% maxNumCompThreads(numCores);
% if isempty(gcp('nocreate'))
%     parpool(numCores-1);
% end

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[gitDir, ~]= fileparts(fileparts(projectDir));
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'bads')));
addpath(genpath(fullfile(gitDir, 'ibs')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, currModelStr);

%% organize data

num_ses  = 9;
for  ses  = 1:num_ses
    data(ses) = organizeData(sub, ses);
    pre_resp(ses,:) = reshape(data(ses).pre_r_org, 1, []);
    post_resp(ses,:) = reshape(data(ses).post_r_org, 1, []);
end
vec_resp = [reshape(pre_resp, 1, numel(pre_resp)), reshape(post_resp, 1, numel(post_resp))]';

%% define model

% set fixed & set-up parameters
model.vec_resp = vec_resp;
model.num_ses = num_ses;
model.test_soa = data(1).pre_ms_unique;
model.n_trial = data(1).pre_numTrials; % number of trials per test soa in pre/post test phase
model.n_expo_trial = 250; % number of *real* trials in exposure phase
model.adaptor_soa =[data(1:model.num_ses).adaptor_soa]; % unsorted adaptor soa in each session
model.num_runs = 4; % fit the model multiple times, each with a different initialization
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.5*1e3; % in second, where measurements are likely to reside

% for simulation
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 1; % simulate pmf with finer axis
model.adaptor_axis_finer = 0;

% save all model information
model.model_info = model_info;
model.i_model = i_model; % current model index
model.currModelStr = currModelStr; % current model folder

% options for BADS
OPTIONS = bads('defaults');
OPTIONS.SpecifyTargetNoise = true;

% options for ibs
options_ibs = ibslike('defaults');
options_ibs.ReturnStd  = true;

%% model fitting using IBS and VMBC

currModel = str2func(['sim_' currModelStr]);

model.mode = 'initialize';
Val = currModel([], [], model);
model.initVal = Val;

model.mode  = 'optimize';
nll_ibs = @(param) ibslike(currModel, param, vec_resp, [], options_ibs, model);

NLL         = NaN(1, model.num_runs);
estP        = NaN(model.num_runs, Val.num_para);

% p = [35.693359375 28.5693359375 77.12890625 47.16796875 0.027138671875, 0.528515625 0.00128564453125 87.1825740824543 700];
% test = currModel(p, model);
% test_nll = nll_ibs(p);

for i  = 1:model.num_runs
    fprintf('[%s] Start fitting model-%s sub-%i run-%i \n', mfilename, currModelStr, sub, i);

    tempModel            = model;
    tempVal              = Val;
    tempFunc             = currModel;

    [estP(i,:),NLL(i)] = bads(nll_ibs,...
        tempVal.init(i,:), tempVal.lb,...
        tempVal.ub, tempVal.plb, tempVal.pub, [], OPTIONS);
end

model.estP            = estP;
model.NLL             = NLL;

% find the parameter with the least NLL
[model.minNLL, best_idx] = min(NLL);
bestP = estP(best_idx, :);
model.bestP = bestP;

%% model prediction by best-fitting parameters

pred =  currModel(bestP, model);

%% save the data for each participant
save(fullfile(outDir, sprintf('sub-%i_%s', sub, datestr(datetime('now')))),'data','model','pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

% end