% This script fits all atheoretical models to the TOJ data, saves the model
% evidence, parameter estimates, and model predictions.

clear all; close all; clc;

%% Select models

rng('shuffle');
specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion', 'Gaussian likelihood, shift bias'};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
num_cores = feature('numcores');
if isempty(gcp('nocreate')); parpool(num_cores-1); end

%% Manage paths

restoredefaultpath;
current_dir = pwd;
[project_dir, ~] = fileparts(current_dir);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(git_dir, 'vbmc')));

%% Model setup

% Set fixed & set-up parameters
model.num_ses = 9;
model.num_runs = num_cores - 1; % Fit the model multiple times, each with a different initialization
model.bound = 10; % Bound for prior axis in seconds
model.bound_int = 1.5; % Likely bound for estimates in seconds
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5] * 1e3; % Test SOA in ms
model.sim_adaptor_soa = [-0.7, -0.3:0.1:0.3, 0.7] * 1e3; % Simulated adaptor SOA in ms
model.test_axis_finer = 1; % Simulate with finer axis
model.model_info = model_info; % Save all model information

% Set fitting options
options = vbmc('defaults');
options.TolStableCount = 30;

%% Fit model

n_sub = 10;

for i_model = 1:numel(folders)

    curr_model_str = model_info.FolderName{i_model};
    addpath(genpath(fullfile(current_dir, curr_model_str)));
    out_dir = fullfile(current_dir, curr_model_str);
    model.curr_model_str = curr_model_str; % Current model folder

    for sub = 1:n_sub

        %% Organize data

        for ses = 1:model.num_ses
            data(ses) = organize_data(sub, ses);
        end

        %% Set model

        curr_model = str2func(['nll_' curr_model_str]);

        % Initialize starting points
        model.mode = 'initialize';
        val = curr_model([], model, data);
        model.init_val = val;

        % Set priors
        lprior_fun = @(x) msplinetrapezlogpdf(x, val.lb, val.plb, val.pub, val.ub);

        % Set likelihood
        model.mode = 'optimize';
        ll_fun = @(x) curr_model(x, model, data);

        fun = @(x) ll_fun(x) + lprior_fun(x);

        [elbo, elbo_sd, exitflag] = deal(NaN(1, model.num_runs));

        parfor i = 1:model.num_runs

            fprintf('[%s] Start fitting model-%s sub-%i run-%i \n', mfilename, curr_model_str, sub, i);
            temp_val = val;

            % vp: variational posterior
            % elbo: Variational Evidence Lower Bound
            [temp_vp{i}, elbo(i), elbo_sd(i), exitflag(i), temp_output{i}] = vbmc(fun, temp_val.init(i,:), temp_val.lb, temp_val.ub, temp_val.plb, temp_val.pub, options);

        end

        % Save all outputs
        model.vp = temp_vp;
        model.elbo = elbo;
        model.elbo_sd = elbo_sd;
        model.exitflag = exitflag;
        model.output = temp_output;

        % Save in case diagnosis fails
        save(fullfile(out_dir, sprintf('sub-%i_%s', sub, curr_model_str)), 'data', 'model')

        %% Evaluate fits

        try
            [diag.exitflag, diag.best_elbo, diag.idx_best, diag.stats] = vbmc_diagnostics(temp_vp);
            diag.best_elcbo = diag.best_elbo.elbo - 3 * diag.best_elbo.elbo_sd;

            % Find best-fitting parameters
            diag.Xs = vbmc_rnd(diag.best_elbo.vp, 1e5);
            diag.post_mean = mean(diag.Xs, 1);

            % Model prediction by best-fitting parameters
            model.mode = 'predict';
            pred = curr_model(diag.post_mean, model, data);

        catch
            fprintf('[%s] No solution has converged. Skip model prediction. \n', mfilename);
            continue;
        end

        %% Save the data for each participant

        save(fullfile(out_dir, sprintf('sub-%i_%s', sub, curr_model_str)), 'data', 'model', 'diag', 'pred')

    end

    rmpath(genpath(fullfile(current_dir, curr_model_str)));

end
% Delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end