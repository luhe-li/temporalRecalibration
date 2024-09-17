% part 1 of model recovry: simulate 100 datesets from 6 models

clear; close all; clc;

%% select models

rng('shuffle'); rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric', 'Causal inference, symmetric','Fixed updated, asymmetric', 'Fixed updated, symmetric'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% set environment
useCluster = false;

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

% job/sample = 100, core/sim_m = 6, fit once
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores  = maxNumCompThreads;
        end
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);
        i_sample = hpc_job_number;

        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores);
        end

    case false
        numCores = 6;
        i_sample = 100;
end

%% manage paths

% restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
outDir = fullfile(currentDir, mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end
if useCluster == false; projectDir = dataDir; end

%% organize data
sub_slc = [1:4,6:10];
for sess = 1:9
    data(sess) = organizeData(1, sess);
end

%% define model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 10;
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators

%% sample ground-truth from best parameter estimates

for sim_m = 1:numel(folders)

    sim_str = folders{sim_m};
    sim_func = str2func(['nll_' sim_str]);
    addpath(genpath(fullfile(pwd, sim_str)));

    model.mode = 'initialize';
    Val = sim_func([], model, []);

    % load best estimates
    clearvars bestP
    result_folder = fullfile(projectDir, 'recalibration_models_VBMC', sim_str);
    R = load_subject_data(result_folder, sub_slc, 'sub-*');
    for ss = 1:numel(sub_slc)
        bestP(ss,:) = R{ss}.diag.post_mean;
    end

    % sample 100 ground truth
    mu_GT = mean(bestP, 1);
    sd_GT = std(bestP, [], 1);
    GT_samples = generate_samples(Val, mu_GT, sd_GT, 1);
    model.mode = 'predict';

    temp_sim_func = sim_func;
    pred =  temp_sim_func(GT_samples, model, []);
    sim_data = simulateData(pred, data);
    fake_data(sim_m, i_sample).data = sim_data;
    fake_data(sim_m, i_sample).gt_p = GT_samples;
    fake_data(sim_m, i_sample).mu_GT = mu_GT;
    fake_data(sim_m, i_sample).sd_GT = sd_GT;
    fake_data(sim_m, i_sample).pred = pred;

end

save(fullfile(outDir, sprintf('sim_data_sample-%02d',i_sample)),'fake_data')

% %% optinal: check predictions
% 
% i_sample = 1;
% for sim_m = 1:numel(folders)
%     figure; 
%     sgtitle(folders{sim_m})
%     subplot 121
%     plot(fake_data(sim_m, i_sample).pred.adaptor_soa, mean(fake_data(sim_m, i_sample).pred.pss_shift,2))
%     subplot 122
%     plot(fake_data(sim_m, i_sample).pred.test_soa, fake_data(sim_m, i_sample).pred.pre_pmf')
% end
