% check which parameter trades off with p_common by simulating NLL while
% varying other parameters

clear; close all; clc; rng('shuffle')

%% set up cluster


model_str = 'cauInf_asym';
currModel = str2func(['nll_' model_str]);
useCluster = true;

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

switch useCluster
    case true
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);

    case false

        hpc_job_number = 1;

end
%% order for job number

num_p = 10; % number of grid of each parameter
job_number = 1;
jobs = zeros(num_p*num_p, 3); % Preallocate array for job number, ii, jj

% Generate pairs of (ii, jj) for each job
for ii = 1:num_p
    for jj = 1:num_p
        jobs(job_number, :) = [job_number, ii, jj];
        job_number = job_number + 1;
    end
end

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
if useCluster
    dataDir = fullfile(projectDir);
else
    dataDir = fullfile(tempDir,'temporalRecalibrationData');
end
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load recalibration model results

sub_slc = [1:4, 6:10];
save_fig = 1;

result_folder = fullfile(dataDir, 'recalibration_models_VBMC', model_str);
addpath(genpath(result_folder));
R = load_subject_data(result_folder, sub_slc, 'sub-*');

for ss = 1:numel(sub_slc)
    p = R{ss}.diag.post_mean;
    bestP{ss} = p;
    data{ss} = R{ss}.data;
end

lb = R{1}.model.initVal.lb;
ub = R{1}.model.initVal.ub;
paraID = R{1}.model.initVal.paraID;

%% set up model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators
model.mode       = 'optimize';

%% use participant-wise best estimates of parameters, vary p_common, alpha, sigma_C1, sigma_c2

para_combi    = [6,7; 6,8; 6,9]; % index of each parameter
ii = jobs(hpc_job_number, 2);
jj = jobs(hpc_job_number, 3);
NLL = zeros(size(para_combi, 1), numel(sub_slc));

for p_comb            = 1:size(para_combi, 1)

    idx_p1           = para_combi(p_comb,1);
    idx_p2           = para_combi(p_comb,2);

    % make a grid for two selected parameters
    p1s              = linspace(lb(idx_p1), ub(idx_p1), num_p+2);
    p2s              = linspace(lb(idx_p2), ub(idx_p2), num_p+2);

    % exclude bounds
    p1s = p1s(2:end-1);
    p2s = p2s(2:end-1);
    p1ss(p_comb,:) = p1s;
    p2ss(p_comb,:) = p2s;

    for ss = 1:numel(sub_slc)

        fprintf('Start simulation for sim-%03d, parameter pair-%02d, sub-%02d \n', hpc_job_number, p_comb, ss)
        sim_p = bestP{ss};
        sim_p(idx_p1) = p1s(ii);
        sim_p(idx_p2) = p2s(jj);
        LL = currModel(sim_p, model, data{ss});
        NLL(p_comb, ss) = -LL;

    end

end

save(fullfile(out_dir, sprintf('sim_num%03d_i%02d_j%02d', hpc_job_number, ii, jj)),'NLL','para_combi',"paraID","p1ss","p2ss");

