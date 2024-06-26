% function fit_recal_model(i_model, useCluster)

i_model = 3;

%% select models

rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric', 'Atheoretical'}; % Column 2: specifications
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym', 'atheo'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
currModelStr = model_info.FolderName{i_model};

%% set environment

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
        numJob    = numel(hpc_job_number);
        fprintf('job number: %i \n', hpc_job_number);
        sub  = hpc_job_number;
    case false
        numCores = feature('numcores');
        sub = 1;
%         sub = str2num(input('Enter Subject ID to analyze: ', 's'));
end

% % make sure Matlab does not exceed this
% fprintf('Number of cores: %i  \n', numCores);
% maxNumCompThreads(numCores);
% if isempty(gcp('nocreate'))
%     parpool(numCores-1);
% end

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'bads')));
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
model.num_runs = numCores-1; % fit the model 100 times, each with a different initialization
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

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 1;
OPTIONS.TolMesh = 1e-3;
% if useCluster; OPTIONS.Display = 'off'; end

%% model fitting

currModel = str2func(['nll_' currModelStr]);

model.mode = 'initialize';
Val = currModel([], model, data);
model.initVal = Val;

model.mode                  = 'optimize';
NLL                         = NaN(1, model.num_runs);
estP                        = NaN(model.num_runs, Val.num_para);

p = [53.857421875 86.15234375 39.23828125 54.66796875 0.00470703125 0.923046875 0.00622314453125 184.96512220759 376.953125];
test = currModel(p, model, data);

for i  = 1%:model.num_runs
    fprintf('[%s] Start fitting model-%s sub-%i run-%i \n', mfilename, currModelStr, sub, i);
    try
        tempModel            = model;
        tempVal              = Val;
        tempFunc             = currModel;

        [estP(i,:),NLL(i)] = bads(@(p) tempFunc(p, tempModel, data),...
            tempVal.init(i,:), tempVal.lb,...
            tempVal.ub, tempVal.plb, tempVal.pub, [], OPTIONS);
    catch
        sprintf('Skipped invalid NLL\n')
        continue;
    end
end

model.estP            = estP;
model.NLL             = NLL;

% find the parameter with the least NLL
[model.minNLL, best_idx] = min(NLL);
bestP = estP(best_idx, :);
model.bestP = bestP;

%% model prediction by best-fitting parameters

model.mode       = 'predict';
pred =  currModel(bestP, model, data);

%% save the data for each participant
save(fullfile(outDir, sprintf('sub-%i_%s', i_sub, datestr(datetime('now')))),'data','model','pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

% end