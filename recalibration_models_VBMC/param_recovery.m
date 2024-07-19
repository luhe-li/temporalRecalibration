% parameter recovery of the causal-inference, asymmetrical likelihood model

clear; close all; rng('shuffle');

%% set environment

currModelStr = 'cauInf_asym';
useCluster = false;

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

% job/sample = 100, core/run = 3, fit once
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores  = maxNumCompThreads;
        end
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);
        recov_run = hpc_job_number;

        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores);
        end

    case false
        numCores = feature('numcores');
        recov_run = 1;
end

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[git_dir, ~] = fileparts(projectDir);
dataDir = fullfile(git_dir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end
if ~useCluster; projectDir = dataDir; end

%% organize data

for ses = 1:9
    data(ses) = organizeData(1, ses);
end

%% define model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = numCores; % fit the model multiple times, each with a different initialization
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators
model.currModelStr = currModelStr; % current model folder

% set OPTIONS
options = vbmc('defaults');
options.MaxFunEvals = 500;
options.TolStableCount = 15;

%% sample ground-truth from best parameter estimates

sub_slc = [1:4,6:10];
result_folder = fullfile(projectDir, 'recalibration_models_VBMC', currModelStr);
R = load_subject_data(result_folder, sub_slc, 'sub-*');
for ss = 1:numel(sub_slc)
    bestP(ss,:) = R{ss}.diag.post_mean;
end
mu_GT = mean(bestP, 1);
sd_GT = std(bestP, [], 1);
num_sample = 1; % for each job, sample once

currModel = str2func(['nll_' currModelStr]);
model.mode = 'initialize';
Val = currModel([], model, []);
GT_samples = generate_samples(Val, mu_GT, sd_GT, num_sample);
fprintf('[%s] GT: %.2f', mfilename, GT_samples);

%% simulation

model.mode       = 'predict';
pred =  currModel(GT_samples, model, []);

% simulate fake data using predicted PMF
fake_data = simulateData(pred, data);

%% fitting
model.mode = 'initialize';
Val = currModel([], model, fake_data);
model.initVal = Val;

% set priors
lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);

% set likelihood
model.mode = 'optimize';
llfun = @(x) currModel(x, model, fake_data);
fun = @(x) llfun(x) + lpriorfun(x);

[elbo,elbo_sd,exitflag] = deal(NaN(1,model.num_runs));

parfor i  = 1:model.num_runs
    
        fprintf('[%s] Start fitting recovery sample-%i run-%i \n', mfilename, recov_run, i);
        tempVal = Val;

        % vp: variational posterior
        % elbo: Variational Evidence Lower Bound
        [temp_vp{i}, elbo(i),elbo_sd(i),exitflag(i),temp_output{i}] = vbmc(fun, tempVal.init(i,:), tempVal.lb,...
            tempVal.ub, tempVal.plb, tempVal.pub, options);

end

% save all outputs
model.vp = temp_vp;
model.elbo = elbo;
model.elbo_sd = elbo_sd;
model.exitflag = exitflag;
model.output = temp_output;

% save in case diagnosis fails
save(fullfile(outDir, sprintf('sample-%02d', recov_run)),'fake_data','model')

% evaluate fits
[diag.exitflag, diag.bestELBO, diag.idx_best, diag.stats] = vbmc_diagnostics(temp_vp);
diag.bestELCBO = diag.bestELBO.elbo - 3*diag.bestELBO.elbo_sd;

% find best-fitting parameters
diag.Xs = vbmc_rnd(diag.bestELBO.vp,1e5);
diag.post_mean = mean(diag.Xs,1);

%% summarize for model recovery plots

summ.gt = GT_samples;
summ.est = diag.post_mean;

%% save the data for each fit
save(fullfile(outDir, sprintf('sample-%02d', recov_run)),'fake_data','model', 'diag','summ')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

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
        case '\tau'
            % Case 1: For 'tau'
            samples(i, :) = normrnd(mean_values(i), sd_values(i), [1, num_sample]);

        case {'\sigma_{A}', '\sigma_{V}', '\sigma', 'c', '\sigma_{C=1}', '\sigma_{C=2}'}
            % Case 2: For 'sigma_a', 'sigma_v', 'sigma', 'criterion', 'sigma_C1', 'sigma_C2'
            v = sd_values(i).^2;
            log_mu = log((mean_values(i).^2) ./ sqrt(v + mean_values(i).^2));
            log_sigma = sqrt(log(v ./ (mean_values(i).^2) + 1));
            samples(i, :) = lognrnd(log_mu, log_sigma, [1, num_sample]);

        case {'p_{common}', '\lambda', '\alpha'}
            % Case 3: For 'p_{common}', '\lambda', '\alpha'
            param_samples = normrnd(mean_values(i), sd_values(i), [1, num_sample]);
            param_samples(param_samples < Val.lb(i)) = Val.lb(i);
            param_samples(param_samples > Val.ub(i)) = Val.ub(i);
            samples(i, :) = param_samples;

        otherwise
            error('Unknown parameter ID: %s', paraId);
    end
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
