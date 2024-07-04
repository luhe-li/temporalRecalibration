function fit_atheo_model_VBMC(i_model, useCluster)

%% select models

rng('Shuffle');
specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion',  'Gaussian likelihood, shift bias',};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
currModelStr = model_info.FolderName{i_model};

% set environment

if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster = false;
end

% job = model, core = run, each job runs all subjects
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores         = maxNumCompThreads;
        end

        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores-1);
        end

    case false
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

%% model set up

% subject and session
n_sub          = 10;

% set fixed & set-up parameters
model.num_ses    = 9;
model.num_runs   = numCores-1; % fit the model multiple times, each with a different initialization
model.bound      = 10; % in second, the bound for prior axis
model.bound_int  = 1.5; % in second, where estimates are likely to reside
model.test_soa   = [-0.5, -0.3:0.05:0.3, 0.5]*1e3; % in ms
model.sim_adaptor_soa = [-0.7, -0.3:0.1:0.3, 0.7]*1e3; % in ms
model.test_axis_finer = 1; % simulate with finer axis
model.model_info = model_info; % save all model information
model.i_model = i_model; % current model index
model.currModelStr = currModelStr; % current model folder

% set OPTIONS
options = vbmc('defaults');
% options.MaxFunEvals = 50; % for debug

%% fit model

for i_sub = 9:10%1:n_sub

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
    lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);

    % set likelihood
    model.mode = 'optimize';
    llfun = @(x) currModel(x, model, data);

    fun = @(x) llfun(x) + lpriorfun(x);

    %     p = [-0.269469843924707 -0.373644351395858 0.0506209680419048 -0.0506209680419047 0.153135114221747 -0.373644351395858 -0.259666155363739 0.259666155363739 0.373644351395858 0.153135114221747 0.187225591143454 0.089219381856417 0.336266931803151 0.21974053784385] ;
    %     test = fun(p);

    [elbo,elbo_sd,exitflag] = deal(NaN(1,model.num_runs));

    parfor i  = 1:model.num_runs

        fprintf('[%s] Start fitting model-%s sub-%i run-%i \n', mfilename, currModelStr, i_sub, i);
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

    % pick the variational solution with highest ELCBO (lower confidence bound on the ELBO)
    model.beta_lcb = 3;       % Standard confidence parameter
    elcbo = elbo - model.beta_lcb*elbo_sd;
    [model.maxELCBO,idx_best] = max(elcbo);

    % find the posterior with higherst model evidence and take its mode as
    % the best-fitting parameter
    model.best_vp = temp_vp{idx_best};
    model.bestP = vbmc_mode(model.best_vp);

    %% model prediction by best-fitting parameters

    model.mode       = 'predict';
    pred =  currModel(model.bestP, model, data);

    %% save the data for each participant

    save(fullfile(outDir, sprintf('sub-%i_%s', i_sub, datestr(datetime('now')))),'data','model','pred')

end

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end
