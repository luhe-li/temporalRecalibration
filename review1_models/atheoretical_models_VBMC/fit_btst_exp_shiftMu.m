clear all; close all; clc;

%% select best-fitting atheoretical model

currModelStr = 'exp_shiftMu';

%% set environment

% exclude outlier
sub_slc = [1:4, 6:10];

if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster= false;
end

% sbatch = model, job = sub, core = run
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores  = maxNumCompThreads;
        end
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);
        sub  = sub_slc(hpc_job_number);

        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores-1);
        end

    case false
        % for local debug
        numCores = feature('numcores');
        sub = 1;
end

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, currModelStr);

%% define model

% set fixed & set-up parameters
model.thres_R2 = 0.95;
model.n_btst = 1000;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 3; % fit the model multiple times, each with a different initialization
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.test_axis_finer = 1; % simulate pmf with finer axis
model.adaptor_axis_finer = 0;
model.currModelStr = currModelStr; % current model folder
model.num_ses = 9;

% set OPTIONS
options = vbmc('defaults');
% options.MaxFunEvals = 100; % for debug
options.TolStableCount = 30;

%% model fitting

currModel = str2func(['nll_' currModelStr]);

% initialize starting points
model.mode = 'initialize';
Val = currModel([], model, []);
model.initVal = Val;

% set priors
lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);

% Preallocate arrays
temp_vp = cell(model.n_btst, model.num_runs);
elbo = zeros(model.n_btst, model.num_runs);
elbo_sd = zeros(model.n_btst, model.num_runs);
exitflag = zeros(model.n_btst, model.num_runs);
temp_output = cell(model.n_btst, model.num_runs);
btstData = struct([]);

for i = 1:model.n_btst

    tempModel = model;
    tempCurrModel = currModel;

    % Bootstrap data
    i_btst_data = bootstrapData(sub);
    btstData(i).data = i_btst_data;

    % Set likelihood
    tempModel.mode = 'optimize';
    llfun = @(x) tempCurrModel(x, tempModel, i_btst_data);
    fun = @(x) llfun(x) + lpriorfun(x);

    % Nested for-loops are problematic in parfor, so we handle them separately
    temp_temp_vp = cell(1, tempModel.num_runs);
    temp_elbo = zeros(1, tempModel.num_runs);
    temp_elbo_sd = zeros(1, tempModel.num_runs);
    temp_exitflag = zeros(1, tempModel.num_runs);
    temp_temp_output = cell(1, tempModel.num_runs);

    for rr = 1:tempModel.num_runs

        fprintf('[%s] Start fitting sub-%i bootstrap run-%i repeat run-%i \n', mfilename, sub, i, rr);
        tempVal = Val;

        % vp: variational posterior
        % elbo: Variational Evidence Lower Bound
        [temp_temp_vp{rr}, temp_elbo(rr), temp_elbo_sd(rr), temp_exitflag(rr), temp_temp_output{rr}] = vbmc(fun, tempVal.init(i,:), tempVal.lb, ...
            tempVal.ub, tempVal.plb, tempVal.pub, options);

    end

    % Store results in the preallocated arrays
    temp_vp(i, :) = temp_temp_vp;
    elbo(i, :) = temp_elbo;
    elbo_sd(i, :) = temp_elbo_sd;
    exitflag(i, :) = temp_exitflag;
    temp_output(i, :) = temp_temp_output;

end

% Save all outputs
model.vp = temp_vp;
model.elbo = elbo;
model.elbo_sd = elbo_sd;
model.exitflag = exitflag;
model.output = temp_output;

% Save incase diagnosis fails
save(fullfile(outDir, sprintf('btst_sub-%i', sub)),'btst_data','model')

%% Evaluate fits for each bootstrap run

% Preallocate arrays for diagnostics
diagExitFlag = zeros(1, model.n_btst);
bestELBO = struct([]);
idx_best = zeros(1, model.n_btst);
stats = struct([]);
bestELCBO = zeros(1, model.n_btst);
bestP = cell(1, model.n_btst);
pred = cell(1, model.n_btst);

for jj = 1:model.n_btst

    tempCurrModel = currModel;
    tempModel = model;

    try
    % Evaluate fits
    [diagExitFlag(jj), bestELBO{jj}, idx_best(jj), stats{jj}] = vbmc_diagnostics(temp_vp(jj, :));
    
    bestELCBO(jj) = bestELBO{jj}.elbo - 3 * bestELBO{jj}.elbo_sd;

    % Find best-fitting parameters
    Xs = vbmc_rnd(bestELBO{jj}.vp, 1e5);
    post_mean = mean(Xs, 1);
    bestP{jj} = post_mean;

    % Model prediction by best-fitting parameters
    tempModel.mode = 'predict';
    pred{jj} = tempCurrModel(post_mean, tempModel, btstData(jj).data);

    catch
            sprintf('No solution has converged. Skip model prediction. \n')
            continue;
    end

end

% Save diagnosis results
diag.diagExitFlag = diagExitFlag;
diag.bestELBO = bestELBO;
diag.idx_best = idx_best;
diag.stats = stats;
diag.bestELCBO = bestELCBO;
diag.bestP = bestP;

%% save the data for each participant
save(fullfile(outDir, sprintf('btst_sub-%i', sub)),'btst_data','model','diag','pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

