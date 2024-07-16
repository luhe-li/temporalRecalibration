
clear; close all; clc;

%% select best-fitting atheoretical model

currModelStr = 'exp_shiftMu';

%% set environment

% exclude outlier
sub_slc = [1:4, 6:10];
useCluster = 1;

if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster= false;
end

% sba job = sub, core = run
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

% restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, currModelStr);

% Save incase diagnosis fails
load(fullfile(outDir, sprintf('btst_sub-%02d', sub)),'btst_data','model')

%% Evaluate fits for each bootstrap run
% Preallocate arrays for diagnostics
currModel = str2func(['nll_' currModelStr]);
diagExitFlag = zeros(1, model.n_btst);
bestELBO = cell(1, model.n_btst); % Initialize as cell array
idx_best = zeros(1, model.n_btst);
stats = cell(1, model.n_btst); % Initialize as cell array
bestELCBO = zeros(1, model.n_btst);
bestP = cell(1, model.n_btst);
pred = cell(1, model.n_btst);

parfor jj = 1:model.n_btst

    tempCurrModel = currModel;
    tempModel = model;

    % Evaluate fits
    [diagExitFlag(jj), bestELBO{jj}, idx_best(jj), stats{jj}] = vbmc_diagnostics(model.vp(jj, :));
    bestELCBO(jj) = bestELBO{jj}.elbo - 3 * bestELBO{jj}.elbo_sd;

    % Find best-fitting parameters
    Xs = vbmc_rnd(bestELBO{jj}.vp, 1e5);
    post_mean = mean(Xs, 1);
    bestP{jj} = post_mean;

    % Model prediction by best-fitting parameters
    tempModel.mode = 'predict';
    pred{jj} = tempCurrModel(post_mean, tempModel, []);

end

% Save diagnosis results
diag.diagExitFlag = diagExitFlag;
diag.bestELBO = bestELBO;
diag.idx_best = idx_best;
diag.stats = stats;
diag.bestELCBO = bestELCBO;
diag.bestP = bestP;

%% save the data for each participant
save(fullfile(outDir, sprintf('diag_btst_sub-%02d', sub)),'btst_data','model','diag','pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end