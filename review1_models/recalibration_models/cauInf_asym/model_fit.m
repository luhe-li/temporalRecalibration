
clear; clc; close all; rng('Shuffle');

%% set environment

useCluster= 0;
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster= false;
end

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
        sub  = str2num(input('Enter Subject ID to analyze: ', 's'));
end

% make sure Matlab does not exceed this
fprintf('Number of cores: %i  \n', numCores);
maxNumCompThreads(numCores);
if isempty(gcp('nocreate'))
    parpool(numCores-1);
end

%% manage paths

currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'bads')));
addpath(genpath(fullfile(projectDir, 'utils')));

%% organize data

model.num_ses  = 9;
for  ses  = 1:model.num_ses
    data(ses) = organizeData(sub, ses);
end

%% define model

% set mode
model.mode= 'optimize';

% set fixed & set-up parameters
model.thres_R2 = 0.95;
model.expo_num_sim= 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = numCores-1; % fit the model 100 times, each with a different initialization
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound    = 10; % in second, the bound for prior axis
model.bound_int= 1.5; % in second, where estimates are likely to reside
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.paraID   = {'\tau','\sigma_{A}','\sigma_{V}','c','\lambda','p_{common}','\alpha','\sigma_{C=1}','\sigma_{C=2}'};
model.num_para = length(model.paraID);

% hard bounds, the range for LB, UB, larger than soft bounds
paraH.tau      = [-100,   100]; % ms
paraH.sigma_a  = [  10,   120]; % ms
paraH.sigma_v  = [  10,   200]; % ms
paraH.criterion= [   1,   350]; % criterion, s
paraH.lambda   = [1e-4,  0.06]; % percentage
paraH.p_common = [1e-4, 1-1e-4]; % weight
paraH.alpha    = [1e-4,  0.02]; % percentage
paraH.sigma_C1 = [   1,   300]; % ms
paraH.sigma_C2 = [ 100,   1e3]; % ms

% soft bounds, the range for PLB, PUB
paraS.tau = [ -50,    50]; % ms
paraS.sigma_a  = [  20,    50]; % ms
paraS.sigma_v  = [  20,   120]; % ms
paraS.criterion= [  30,   150]; % criterion, s
paraS.lambda   = [0.01,  0.03]; % percentage
paraS.p_common = [ 0.3,   0.7]; % weight
paraS.alpha    = [1e-3,  2e-3]; % percentage
paraS.sigma_C1 = [  10,   100]; % ms
paraS.sigma_C2 = [ 500,   700]; % ms

% reorganize parameter bounds to feed to bads
fn= fieldnames(paraH);
for k= 1:numel(fn)
    model.lb(:,k)  = paraH.(fn{k})(1);
    model.ub(:,k)  = paraH.(fn{k})(2);
    model.plb(:,k) = paraS.(fn{k})(1);
    model.pub(:,k) = paraS.(fn{k})(2);
end
model.paraS    = paraS; model.paraH = paraH;

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 1;

%% model fitting

% define function calculating nll
funcNLL   = @(p) nll_CI_asym(p(1), p(2), p(3), p(4), p(5), p(6),p(7),p(8),p(9),...
    model, data);

% get grid initializations
numSections = model.num_runs * 2;
model.init= getInit(model.lb, model.ub, numSections, model.num_runs);

%initialize matrices that store negative log likelihood and best-fit paramters
minNLL    = NaN(1, model.num_runs);
estimatedP= NaN(model.num_runs, length(model.lb));

% test with a set of parameters if needed
% %    tau, sigma_a, sigma_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2
% p = [0.7340 82.9118 111.6321 121.8056 0.0600 0.0983 0.0065 1.1886 290.9797];
% testnll= BCImodel(p(1), p(2), p(3), p(4), p(5), p(6),p(7),p(8),p(9), model, data);

parfor i  = 1:model.num_runs
    fprintf('[%s] Start fitting sub-%i run-%i \n', mfilename, i_sub, i);
    try
        tempModel = model;
        [estimatedP(i,:),minNLL(i)] = bads(funcNLL, tempModel.init(i,:), tempModel.lb,...
            tempModel.ub, tempModel.plb, tempModel.pub, [], OPTIONS);
        disp(estimatedP(i,:))
    catch
        sprintf('Skipped invalid NLL\n')
        continue;
    end
end

% store the estimated p with min nll
[min_val, min_idx]   = min(minNLL);
model.minNLL = min_val;
model.bestP = estimatedP(min_idx,:);

% store all fits
model.NLL = minNLL;
model.estP = estimatedP;

%% model prediction by best-fitting parameters


% save the data for each participant
save(sprintf('sub-%i_%s', i_sub, datestr(datetime('now'))),'data','model','pred')

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
