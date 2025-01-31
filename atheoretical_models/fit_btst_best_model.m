% This script fits the best-fitting atheoretical model to the bootstrapped TOJ data

clear; close all; clc;

%% Select best-fitting atheoretical model

curr_model_str = 'exp_shiftMu';
sub_slc = [1:4, 6:10]; % Exclude outlier
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
addpath(genpath(fullfile(current_dir, curr_model_str)));
out_dir = fullfile(current_dir, curr_model_str);

%% Define model

% Set fixed & set-up parameters
model.thres_R2 = 0.95;
model.n_btst = 1000;
model.expo_num_sim = 1e3; % Number of simulation for exposure phase
model.expo_num_trial = 250; % Number of *real* trials in exposure phase
model.num_runs = 3; % Fit the model multiple times, each with a different initialization
model.num_bin  = 100; % Number of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % In second, the bound for prior axis
model.bound_int = 1.4*1e3; % In second, where measurements are likely to reside
model.num_sample = 1e3; % Number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.test_axis_finer = 1; % Simulate pmf with finer axis
model.adaptor_axis_finer = 0;
model.curr_model_str = curr_model_str;
model.num_ses = 9;

% Set OPTIONS
options = vbmc('defaults');
options.TolStableCount = 30;

%% Model fitting

% Initialize starting points
curr_model = str2func(['nll_' curr_model_str]);
model.mode = 'initialize';
Val = curr_model([], model, []);
model.initVal = Val;

% Set priors
lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);

% Preallocate arrays
temp_vp = cell(model.n_btst, model.num_runs);
elbo = zeros(model.n_btst, model.num_runs);
elbo_sd = zeros(model.n_btst, model.num_runs);
exitflag = zeros(model.n_btst, model.num_runs);
temp_output = cell(model.n_btst, model.num_runs);
btst_data = struct([]);

for sub = sub_slc

    parfor i = 1:model.n_btst

        temp_model = model;
        temp_curr_model = curr_model;

        % Bootstrap data
        i_btst_data = bootstrap_data(sub);
        btst_data(i).data = i_btst_data;

        % Set likelihood
        temp_model.mode = 'optimize';
        llfun = @(x) temp_curr_model(x, temp_model, i_btst_data);
        fun = @(x) llfun(x) + lpriorfun(x);

        temp_temp_vp = cell(1, temp_model.num_runs);
        temp_elbo = zeros(1, temp_model.num_runs);
        temp_elbo_sd = zeros(1, temp_model.num_runs);
        temp_exitflag = zeros(1, temp_model.num_runs);
        temp_temp_output = cell(1, temp_model.num_runs);

        for rr = 1:temp_model.num_runs

            fprintf('[%s] Start fitting sub-%i bootstrap run-%i repeat run-%i \n', mfilename, sub, i, rr);
            temp_val = Val;

            % vp: variational posterior
            % elbo: Variational Evidence Lower Bound
            [temp_temp_vp{rr}, temp_elbo(rr), temp_elbo_sd(rr), temp_exitflag(rr), temp_temp_output{rr}] = vbmc(fun, temp_val.init(rr,:), temp_val.lb, ...
                temp_val.ub, temp_val.plb, temp_val.pub, options);

        end

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
    save(fullfile(out_dir, sprintf('btst_sub-%02d', sub)), 'btst_data', 'model')

    %% Evaluate fits for each bootstrap run

    % Preallocate arrays for diagnostics
    diag_exit_flag = zeros(1, model.n_btst);
    best_elbo = struct([]);
    idx_best = zeros(1, model.n_btst);
    stats = struct([]);
    best_elcbo = zeros(1, model.n_btst);
    best_p = cell(1, model.n_btst);
    pred = cell(1, model.n_btst);

    parfor jj = 1:model.n_btst

        temp_curr_model = curr_model;
        temp_model = model;

        % Evaluate fits
        [diag_exit_flag(jj), best_elbo{jj}, idx_best(jj), stats{jj}] = vbmc_diagnostics(temp_model.vp(jj, :));
        best_elcbo(jj) = best_elbo{jj}.elbo - 3 * best_elbo{jj}.elbo_sd;

        % Find best-fitting parameters
        Xs = vbmc_rnd(best_elbo{jj}.vp, 1e5);
        post_mean = mean(Xs, 1);
        best_p{jj} = post_mean;

        % Model prediction by best-fitting parameters
        temp_model.mode = 'predict';
        pred{jj} = temp_curr_model(post_mean, temp_model, btst_data(jj).data);

    end

    % Save diagnosis results
    diag.diag_exit_flag = diag_exit_flag;
    diag.best_elbo = best_elbo;
    diag.idx_best = idx_best;
    diag.stats = stats;
    diag.best_elcbo = best_elcbo;
    diag.best_p = best_p;

    %% save the data for each participant
    save(fullfile(out_dir, sprintf('btst_sub-%02d', sub)), 'btst_data', 'model', 'diag', 'pred')

end

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
