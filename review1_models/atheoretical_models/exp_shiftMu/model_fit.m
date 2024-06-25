
clear; clc; close all; rng('Shuffle');

%% set up

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster       = false;
end

% see how many cores we have
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores         = maxNumCompThreads;
        end
    case false
        numCores = feature('numcores');
end

% make sure Matlab does not exceed this
fprintf('Number of cores: %i  \n', numCores);
maxNumCompThreads(numCores);
if isempty(gcp('nocreate'))
    parpool(numCores-1);
end

%% manage paths

currentDir       = pwd;
[projectDir, ~]  = fileparts(fileparts(currentDir));
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'bads')));
addpath(genpath(fullfile(projectDir, 'utils')));

%% model set up

% subject and session
n_sub          = 10;
n_ses          = 9;

% set fixed & set-up parameters
model.num_ses    = 9;
model.num_runs   = 100; % fit the model 100 times, each with a different initialization
model.bound      = 10; % in second, the bound for prior axis
model.bound_int  = 1.5; % in second, where estimates are likely to reside
model.test_soa   = [-0.5, -0.3:0.05:0.3, 0.5]*1e3; % in ms
model.sim_adaptor_soa = [-0.7, -0.3:0.1:0.3, 0.7]*1e3; % in ms
model.paraID     = {'\tau_{pre}','\tau_{post,1}', '\tau_{post,2}', '\tau_{post,3}', '\tau_{post,4}', '\tau_{post,5}', '\tau_{post,6}', '\tau_{post,7}', '\tau_{post,8}', '\tau__{post,9}','\sigma_{A}','\sigma_{V}','c','\lambda'};
model.num_para   = length(model.paraID);

% hard bounds, the range for LB, UB, larger than soft bounds
paraH.tau_pre               = [-200,   200]; % ms
paraH.tau1                  = [-200,   200]; % ms
paraH.tau2                  = [-200,   200]; % ms
paraH.tau3                  = [-200,   200]; % ms
paraH.tau4                  = [-200,   200]; % ms
paraH.tau5                  = [-200,   200]; % ms
paraH.tau6                  = [-200,   200]; % ms
paraH.tau7                  = [-200,   200]; % ms
paraH.tau8                  = [-200,   200]; % ms
paraH.tau9                  = [-200,   200]; % ms
paraH.sigma_a               = [  10,   200]; % ms
paraH.sigma_v               = [  10,   200]; % ms
paraH.criterion             = [   1,   350]; % criterion, s
paraH.lambda                = [1e-4,  0.06]; % percentage

% soft bounds, the range for PLB, PUB
paraS.tau_pre               = [-200,   200]; % ms
paraS.tau1                  = [-200,   200]; % ms
paraS.tau2                  = [-200,   200]; % ms
paraS.tau3                  = [-200,   200]; % ms
paraS.tau4                  = [-200,   200]; % ms
paraS.tau5                  = [-200,   200]; % ms
paraS.tau6                  = [-200,   200]; % ms
paraS.tau7                  = [-200,   200]; % ms
paraS.tau8                  = [-200,   200]; % ms
paraS.tau9                  = [-200,   200]; % ms
paraS.sigma_a               = [  20,    50]; % ms
paraS.sigma_v               = [  20,   120]; % ms
paraS.criterion             = [  30,   150]; % criterion, s
paraS.lambda                = [0.01,  0.03]; % percentage

% reorganize parameter bounds to feed to bads
fn    = fieldnames(paraH);
for k = 1:numel(fn)
    model.lb(:,k)    = paraH.(fn{k})(1);
    model.ub(:,k)    = paraH.(fn{k})(2);
    model.plb(:,k)   = paraS.(fn{k})(1);
    model.pub(:,k)   = paraS.(fn{k})(2);
end
model.paraS      = paraS; model.paraH = paraH;

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 1;
OPTIONS.TolMesh = 1e-5;
OPTIONS.Display = 'off';

%% fit model

for i_sub = 1:10

    %% organize data
    for  ses         = 1:model.num_ses
        data(ses)        = organizeData(i_sub, ses);
    end

    %% model fitting

    % define function calculating nll
    model.mode       = 'optimize';
    funcNLL          = @(p) nll_exp_shiftMu(p(1), p(2), p(3), p(4), p(5), p(6),p(7),p(8),p(9),...
        p(10),p(11),p(12),p(13),p(14),...
        model, data);

    % get grid initializations
    model.num_sections   = model.num_runs*2;
    model.init      = getInit(model.lb, model.ub, model.num_sections, model.num_runs);

    %initialize matrices that store negative log likelihood and best-fit paramters
    minNLL = NaN(1, model.num_runs);
    estimatedP = NaN(model.num_runs, length(model.lb));

%     params = num2cell([model.init(1,:)]);  
%     test = nll_exp_shiftMu(params{:}, model, data);

    parfor i         = 1:model.num_runs
        sprintf('[%s] Start fitting run-%i', mfilename, i);
        try
            tempModel = model;
            [estimatedP(i,:), minNLL(i)] = bads(funcNLL, tempModel.init(i,:), tempModel.lb,...
                tempModel.ub, tempModel.plb, tempModel.pub, [], OPTIONS);
        catch
            disp('Skipped invalid NLL')
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
    model.OPTIMSTATE = OPTIMSTATE;

   %% model prediction by best-fitting parameters

   model.mode       = 'predict';
   model.test_axis_finer = 1;
   params = num2cell([model.bestP]);  
   pred =  nll_gauss_shiftMu(params{:}, model, data);

   %% save by subject

   save(sprintf('sub-%i_%s', i_sub, datestr(datetime('now'))),'data','model','pred')

end

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
