function fit_recal_model_VBMC(i_model, useCluster, sub)

%% select models

rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Fixed updated, asymmetric', 'Fixed updated, symmetric'}; 
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
currModelStr = model_info.FolderName{i_model};

%% set environment

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

%% organize data

model.num_ses  = 9;
for  ses  = 1:model.num_ses
    data(ses) = organizeData(sub, ses);
end

%% define model

% set fixed & set-up parameters
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = numCores-1; % fit the model multiple times, each with a different initialization
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3; 
model.toj_axis_finer = 1; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; 

% save all model information
model.model_info = model_info; 
model.i_model = i_model; % current model index
model.currModelStr = currModelStr; % current model folder

% set OPTIONS
options = vbmc('defaults');
if ismember(i_model, [3,5,7]); options.MaxFunEvals = 500; end % set max iter for cau-inf models
options.TolStableCount = 15;

%% model fitting

currModel = str2func(['nll_' currModelStr]);

% initialize starting points
model.mode = 'initialize';
Val = currModel([], model, data);
model.initVal = Val;

% set priors
lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);

% set likelihood
model.mode = 'optimize';
llfun = @(x) currModel(x, model, data);
fun = @(x) llfun(x) + lpriorfun(x);

[elbo,elbo_sd,exitflag] = deal(NaN(1,model.num_runs));

parfor i  = 1:model.num_runs
    
        fprintf('[%s] Start fitting model-%s sub-%i run-%i \n', mfilename, currModelStr, sub, i);
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

% save incase diagnosis fails
save(fullfile(outDir, sprintf('sub-%02d_%s', sub)),'data','model')

% evaluate fits
[diag.exitflag, diag.bestELBO, diag.idx_best, diag.stats] = vbmc_diagnostics(temp_vp);
diag.bestELCBO = diag.bestELBO.elbo - 3*diag.bestELBO.elbo_sd;

% find best-fitting parameters
diag.Xs = vbmc_rnd(diag.bestELBO.vp,1e5);
diag.post_mean = mean(diag.Xs,1);

%% model prediction by best-fitting parameters

model.mode       = 'predict';
pred =  currModel(diag.post_mean, model, data);

%% save the data for each participant
save(fullfile(outDir, sprintf('sub-%02d_%s', sub)),'data','model','diag','pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end