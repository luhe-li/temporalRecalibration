% part 1 of model recovry: simulate 100 datesets from 6 models

clear; close all; clc;

%% select models

rng('shuffle'); rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric', 'Causal inference, symmetric','Fixed updated, asymmetric', 'Fixed updated, symmetric'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% set environment
useCluster = true;

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

%% utility functions
function samples = generate_samples(Val, mean_values, sd_values, num_sample)
% Extract the parameter IDs
paraIds = Val.paraID;

% Initialize the samples matrix
num_parameters = length(paraIds);
samples = zeros(num_parameters, num_sample);

% Loop through each parameter
for i = 1:num_parameters
    paraId = paraIds{i};

    switch paraId
        case {'\tau','p_{common}', '\lambda', '\alpha'}
            % Case 1: For 'tau'
            param_samples = normrnd(mean_values(i), sd_values(i), [1, num_sample]);

        case {'\sigma_{A}', '\sigma_{V}', '\sigma', 'c', '\sigma_{C=1}', '\sigma_{C=2}'}
            % Case 2: For 'sigma_a', 'sigma_v', 'sigma', 'criterion', 'sigma_C1', 'sigma_C2'
            v = sd_values(i).^2;
            log_mu = log((mean_values(i).^2) ./ sqrt(v + mean_values(i).^2));
            log_sigma = sqrt(log(v ./ (mean_values(i).^2) + 1));
            param_samples = lognrnd(log_mu, log_sigma, [1, num_sample]);

        otherwise
            error('Unknown parameter ID: %s', paraId);
    end

    % Ensure all samples are within the bounds
    for j = 1:num_sample
        while param_samples(j) < Val.lb(i) || param_samples(j) > Val.ub(i)
            switch paraId
                case {'\tau','p_{common}', '\lambda', '\alpha'}
                    param_samples(j) = normrnd(mean_values(i), sd_values(i));

                case {'\sigma_{A}', '\sigma_{V}', '\sigma', 'c', '\sigma_{C=1}', '\sigma_{C=2}'}
                    param_samples(j) = lognrnd(log_mu, log_sigma);

                otherwise
                    error('Unknown parameter ID: %s', paraId);
            end
        end
    end

    samples(i, :) = param_samples;
end
end

function fake_data = simulateData(pred, data)

fake_data = data;
nT = data.pre_numTrials;
nLevel = numel(data(1).post_ms_unique);

for ses = 1:9

    % pretest
    M = rand(nT, nLevel);
    bool_V1st = M < repmat(pred.pre_pmf(1,:),[nT, 1]);
    fake_data(ses).pre_nT_V1st = sum(bool_V1st, 1);
    fake_data(ses).pre_pResp(1,:) = sum(bool_V1st, 1)/nT;
    bool_A1st = M > repmat(1 - pred.pre_pmf(3,:),[nT, 1]);
    fake_data(ses).pre_nT_A1st = sum(bool_A1st, 1);
    fake_data(ses).pre_pResp(3,:) = sum(bool_A1st, 1)/nT;
    fake_data(ses).pre_nT_simul = repmat(nT, [1, nLevel]) - fake_data(ses).pre_nT_A1st - fake_data(ses).pre_nT_V1st;
    fake_data(ses).pre_pResp(2,:) = fake_data(ses).pre_nT_simul./nT;

    % posttest
    M = rand(nT, nLevel);
    bool_V1st = M < repmat(squeeze(pred.post_pmf(ses,1,:))',[nT, 1]);
    fake_data(ses).post_nT_V1st = sum(bool_V1st, 1);
    fake_data(ses).post_pResp(1,:) = sum(bool_V1st, 1)/nT;
    bool_A1st = M > repmat(1 - squeeze(pred.post_pmf(ses,3,:))',[nT, 1]);
    fake_data(ses).post_nT_A1st = sum(bool_A1st, 1);
    fake_data(ses).post_pResp(3,:) = sum(bool_A1st, 1)/nT;
    fake_data(ses).post_nT_simul = repmat(nT, [1, nLevel]) - fake_data(ses).post_nT_A1st - fake_data(ses).post_nT_V1st;
    fake_data(ses).post_pResp(2,:) = fake_data(ses).post_nT_simul./nT;

end
end