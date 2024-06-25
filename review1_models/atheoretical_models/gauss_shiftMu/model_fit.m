
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
model.paraID     = {'\mu_{pre}','\sigma','c','\mu_{post,1}','\mu_{post,2}','\mu_{post,3}','\mu_{post,4}','\mu_{post,5}','\mu_{post,6}','\mu_{post,7}','\mu_{post,8}','\mu_{post,9}','\lambda'};
model.num_para   = length(model.paraID);

% hard bounds, the range for LB, UB, larger than soft bounds
paraH.mu         = [-300,   300]; % ms
paraH.sigma      = [0.01,   300]; % ms
paraH.c          = [0.01,   350]; % ms
paraH.mu_post1   = [-300,   300]; % ms
paraH.mu_post2   = [-300,   300]; % ms
paraH.mu_post3   = [-300,   300]; % ms
paraH.mu_post4   = [-300,   300]; % ms
paraH.mu_post5   = [-300,   300]; % ms
paraH.mu_post6   = [-300,   300]; % ms
paraH.mu_post7   = [-300,   300]; % ms
paraH.mu_post8   = [-300,   300]; % ms
paraH.mu_post9   = [-300,   300]; % ms
paraH.lambda     = [1e-4,  0.06]; % percentage

% soft bounds, the range for PLB, PUB
paraS.mu         = [-100,   100]; % ms
paraS.sigma      = [  10,   100]; % ms
paraS.c          = [  10,   100]; % ms
paraS.mu_post1   = [-100,   100]; % ms
paraS.mu_post2   = [-100,   100]; % ms
paraS.mu_post3   = [-100,   100]; % ms
paraS.mu_post4   = [-100,   100]; % ms
paraS.mu_post5   = [-100,   100]; % ms
paraS.mu_post6   = [-100,   100]; % ms
paraS.mu_post7   = [-100,   100]; % ms
paraS.mu_post8   = [-100,   100]; % ms
paraS.mu_post9   = [-100,   100]; % ms
paraS.lambda     = [0.01,  0.03]; % percentage

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

%% fit model

for i_sub = 1:10

    %% organize data
    for  ses         = 1:model.num_ses
        data(ses)        = organizeData(i_sub, ses);
    end

    %% model fitting

    % define function calculating nll
    model.mode       = 'optimize';
    funcNLL          = @(p) nll_gauss_shiftMu(p(1), p(2), p(3), p(4), p(5), p(6),p(7),p(8), p(9),...
        p(10),p(11),p(12),p(13),...
        model, data);

    % get grid initializations
    model.num_sections   = model.num_runs*2;
    model.init      = getInit(model.lb, model.ub, model.num_sections, model.num_runs);

    %initialize matrices that store negative log likelihood and best-fit paramters
    NLL = NaN(1, model.num_runs);
    estimatedP = NaN(model.num_runs, length(model.lb));

%     params = num2cell([model.init(1,:)]);  
%     test = nll_gauss_shiftMu(params{:}, model, data);

    parfor i         = 1:model.num_runs
        disp(i);
        try
            tempModel = model;
            [estimatedP(i,:), NLL(i)] = bads(funcNLL, tempModel.init(i,:), tempModel.lb,...
                tempModel.ub, tempModel.plb, tempModel.pub, [], OPTIONS);
        catch
            disp('Skipped invalid NLL')
            continue;
        end
    end

    % store the estimated p with min nll
    [min_val, min_idx]   = min(NLL);
    model.minNLL = min_val;
    model.bestP = estimatedP(min_idx,:);

    % store all fits
    model.NLL = NLL;
    model.estP = estimatedP;

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
