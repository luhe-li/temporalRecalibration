function model_recovery_s3(fit_m)

% job = sample, core = sim_m

%% select models

rng('shuffle'); rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Fixed updated, asymmetric', 'Fixed updated, symmetric'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
currModelStr = model_info.FolderName{fit_m};

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
        i_sample = 1;
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
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end
if useCluster == false; projectDir = dataDir; end

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
model.fit_m = fit_m; % current model folder
model.fit_m_str = currModelStr; % current model folder

% set OPTIONS
options = vbmc('defaults');
if ismember(fit_m, [3,4]); options.MaxFunEvals = 500; end % set max iter for cau-inf models
options.TolStableCount = 15;

%% load similated data

result_folder = fullfile(projectDir, 'recalibration_models_VBMC', 'model_recovery_s2');
d = load(fullfile(result_folder, 'sim_data.mat')); % load struct `sim`

parfor sim_m = 1:numCores

    temp_d = d;
    temp_model = model;
    sim_data = temp_d.fake_data(sim_m, i_sample).data;
    currModel = str2func(['nll_' currModelStr]);

    temp_model.mode = 'initialize';
    Val = currModel([], temp_model, []);
    temp_model.initVal = Val;

    % set priors
    lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);

    % set likelihood
    temp_model.mode = 'optimize';
    llfun = @(x) currModel(x, temp_model, sim_data);
    fun = @(x) llfun(x) + lpriorfun(x);

    fprintf('[%s] Start sim model-%s, fit model-%s, recovery sample-%i \n', mfilename, folders{sim_m}, currModelStr, i_sample);

    % vp: variational posterior
    % elbo: Variational Evidence Lower Bound
    [temp_vp, temp_elbo, temp_elbo_sd, temp_exitflag,temp_output] = vbmc(fun, Val.init(randi(temp_model.num_runs,1),:), Val.lb,...
        Val.ub, Val.plb, Val.pub, options);

    % save all outputs
    vp{sim_m} = temp_vp;
    elbo(sim_m) = temp_elbo;
    elbo_sd(sim_m) = temp_elbo_sd;
    exitflag(sim_m) = temp_exitflag;
    output{sim_m} = temp_output;
    fake_data{sim_m} = sim_data;
    fake_pred{sim_m} = temp_d.fake_data(sim_m, i_sample).pred;

    % make predictions
    Xs = vbmc_rnd(temp_vp,1e5);
    post_mean = mean(Xs,1);

    %% model prediction by best-fitting parameters

    temp_model.mode       = 'predict';
    pred =  currModel(post_mean, temp_model, sim_data);
    fit_pred{sim_m} = pred;

end

model.vp = vp;
model.elbo = elbo;
model.elbo_sd = elbo_sd;
model.exitflag = exitflag;
model.output = output;

for sim_m = 1:numCores
    summ.bestELCBO(sim_m) = model.elbo(sim_m) - 3*model.elbo_sd(sim_m);
    summ.bestELBO(sim_m) = model.elbo(sim_m);
end

save(fullfile(outDir, sprintf('fitM%02d_sample-%02d_%s', fit_m, i_sample, datestr(datetime('now')))),'fake_data','fake_pred','model','summ','fit_pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

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