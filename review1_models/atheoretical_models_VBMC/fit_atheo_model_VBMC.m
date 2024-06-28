% function fit_atheo_model_VBMC(i_model, useCluster)

i_model = 2;
%% select models

rng('Shuffle');
specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion',  'Gaussian likelihood, shift bias',};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
currModelStr = model_info.FolderName{i_model};

%% set environment

if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster = false;
end

% job = model, core = run, each job runs all subjects
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores         = maxNumCompThreads;
        end
    case false
        numCores = feature('numcores');
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
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(currentDir, currModelStr)));
outDir = fullfile(currentDir, currModelStr);

%% model set up

% subject and session
n_sub          = 10;

% set fixed & set-up parameters
model.num_ses    = 9;
model.num_runs   = 4; % fit the model multiple times, each with a different initialization
model.bound      = 10; % in second, the bound for prior axis
model.bound_int  = 1.5; % in second, where estimates are likely to reside
model.test_soa   = [-0.5, -0.3:0.05:0.3, 0.5]*1e3; % in ms
model.sim_adaptor_soa = [-0.7, -0.3:0.1:0.3, 0.7]*1e3; % in ms
model.test_axis_finer = 1; % simulate with finer axis
model.model_info = model_info; % save all model information
model.i_model = i_model; % current model index
model.currModelStr = currModelStr; % current model folder

% set OPTIONS to tell bads that my objective function is noisy
options = vbmc('defaults');
options.SpecifyTargetNoise = true;
% if useCluster; OPTIONS.Display = 'off'; end

%% fit model

for i_sub = 1:n_sub

    %% organize data

    for  ses         = 1:model.num_ses
        data(ses)        = organizeData(i_sub, ses);
    end

    %% set model

    currModel = str2func(['nll_' currModelStr]);

    % initialize starting points
    model.mode = 'initialize';
    Val = currModel([], model, data);
    model.initVal = Val;

    % set priors
    lpriorfun = @(x) msplinetrapezlogpdf(x,Val.lb,Val.plb,Val.pub,Val.ub);

    % set likelihood
    model.mode                  = 'optimize';
    llfun = @(x) currModel(x, model, data);

    fun = @(x) llfun(x) + lpriorfun(x);

%     testp = [-0.206790284646954 0.340071929623188 -0.0651280006225672 0.0651280006225673 -0.340071929623188 0.198121792259214 -0.340071929623188 0.340071929623188 -0.340071929623188 0.198121792259214 0.200370483380783 -0.257255749371589 -0.5 0.0022695536032044];
%     test = fun(testp);

    for i  = 1:model.num_runs

        fprintf('[%s] Start fitting model-%s sub-%i run-%i \n', mfilename, currModelStr, i_sub, i);
        tempVal              = Val;

        [vp(i),elbo(i),elbo_sd(i)] = vbmc(fun, tempVal.init(i,:), tempVal.lb,...
            tempVal.ub, tempVal.plb, tempVal.pub, options);
    end

%     model.estP            = estP;
%     model.NLL             = NLL;
% 
%     % find the parameter with the least NLL
%     [model.minNLL, best_idx] = min(NLL);
%     bestP = estP(best_idx, :);
%     model.bestP = bestP;
% 
%     %% model prediction by best-fitting parameters
% 
%     model.mode       = 'predict';
%     pred =  currModel(bestP, model, data);

    %% save the data for each participant

    save(fullfile(outDir, sprintf('sub-%i_%s', i_sub, datestr(datetime('now')))),'data','model','pred')

end

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

% end
