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
    fake_gt{sim_m} = temp_d.fake_data(sim_m, i_sample).gt_p;
    fake_mu_gt{sim_m} = temp_d.fake_data(sim_m, i_sample).mu_GT;
    fake_SD_gt{sim_m} = temp_d.fake_data(sim_m, i_sample).sd_GT;

    % make predictions
    Xs = vbmc_rnd(temp_vp,1e5);
    post_mean = mean(Xs,1);
    fit_est{sim_m} = post_mean;

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

%% summariz some metrics

for sim_m = 1:numCores
    
    summ.bestELCBO(sim_m) = model.elbo(sim_m) - 3*model.elbo_sd(sim_m);
    summ.bestELBO(sim_m) = model.elbo(sim_m);
    summ.RMSE(sim_m) = sqrt(sum((mean(fit_pred{sim_m}.pss_shift,2) - mean(fake_pred{sim_m}.pss_shift,2)).^2)/numel(mean(fit_pred{sim_m}.pss_shift,2)));
    
end

save(fullfile(outDir, sprintf('fitM%02d_sample-%02d_%s', fit_m, i_sample, datestr(datetime('now')))),...
    'fake_data','fake_pred','model','summ','fit_est','fit_pred','fake_gt','fit_est','fake_mu_gt','fake_SD_gt')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end
