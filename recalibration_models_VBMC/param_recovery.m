% parameter recovery of the causal-inference, asymmetrical likelihood model

%% set environment

currModelStr = 'cauInf_asym';
useCluster = true;

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
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% organize data

for ses = 1:9
    data(ses) = organizeData(1, ses);
end

%% sample ground-truth from best parameter estimates

sub_slc = [1:4,6:10];
recal_folder = fullfile(projectDir, 'recalibration_models_VBMC', currModelStr);
files = dir(fullfile(recal_folder, 'sub-*'));
for ss = 1:numel(sub_slc)
    i_sub = sub_slc(ss);
    i_data = load(fullfile(recal_folder, files(ss).name));
    bestP(ss,:) = i_data.diag.post_mean;
end
mu_GT = mean(bestP, 1);
sd_GT = sqrt((bestP - mu_GT).^2./numel(sub_slc));
num_sample = 1; % for each job, sample once
GT_samples = sampleGT(mu_GT, sd_GT, num_sample);

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

%% simualte data and fit model

currModel = str2func(['nll_' currModelStr]);

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
llfun = @(x) currModel(x, model, data);
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
save(fullfile(outDir, sprintf('sample-%i_%s', recov_run, datestr(datetime('now')))),'fake_data','model')

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
save(fullfile(outDir, sprintf('sample-%i_%s', recov_run, datestr(datetime('now')))),'fake_data','model', 'diag','summ')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

%% utility functions

function [samples] = sampleGT(mean_values, sd_values, num_sample)

% Initialize the matrix to store the sampled values
samples = zeros(num_sample, 9);

% Parameters '\tau' are normally distributed
samples(:, 1) = normrnd(mean_values(1), sd_values(1), [num_sample, 1]);

% Parameters '\sigma_{A}', '\sigma_{V}', 'criterion', '\sigma_{C1}', '\sigma_{C2}' are log-normally distributed
m = [mean_values(2), mean_values(3), mean_values(4), mean_values(8), mean_values(9)];
v = [sd_values(2:4), sd_values(8), sd_values(9)].^2; % variance
log_mu = log((m.^2)./sqrt(v+m.^2));
log_sigma = sqrt(log(v./(m.^2)+1));
samples(:, 2) = lognrnd(log_mu(1), log_sigma(1), [num_sample, 1]);
samples(:, 3) = lognrnd(log_mu(2), log_sigma(2), [num_sample, 1]);
samples(:, 4) = lognrnd(log_mu(3), log_sigma(3), [num_sample, 1]);
samples(:, 8) = lognrnd(log_mu(4), log_sigma(4), [num_sample, 1]);
samples(:, 9) = lognrnd(log_mu(5), log_sigma(5), [num_sample, 1]);

% Parameter 'p_{common}','\lambda','\alpha' are bounded in a range, sampled from a truncated normal distribution
lambda_samples = normrnd(mean_values(5), sd_values(5), [num_sample, 1]);
p_common_samples = normrnd(mean_values(6), sd_values(6), [num_sample, 1]);
alpha_samples = normrnd(mean_values(7), sd_values(7), [num_sample, 1]);
lambda_samples(lambda_samples < 1e-4) = 1e-4;
lambda_samples(lambda_samples > 0.06) = 0.06;
p_common_samples(p_common_samples < 0) = 0;
p_common_samples(p_common_samples > 1) = 1;
alpha_samples(alpha_samples < 1e-4) = 1e-4;
alpha_samples(alpha_samples > 0.02) = 0.02;

samples(:, 5) = lambda_samples;
samples(:, 6) = p_common_samples;
samples(:, 7) = alpha_samples;

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