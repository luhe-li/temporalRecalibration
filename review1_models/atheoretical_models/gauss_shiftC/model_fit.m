
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
model.num_runs   = 80; % fit the model 100 times, each with a different initialization
model.bound      = 10; % in second, the bound for prior axis
model.bound_int  = 1.5; % in second, where estimates are likely to reside
model.test_soa   = [-0.5, -0.3:0.05:0.3, 0.5]*1e3; % in ms
model.sim_adaptor_soa = [-0.7, -0.3:0.1:0.3, 0.7]*1e3; % in ms
model.paraID     = {'\mu','\sigma','c1_{pre}','dc_{post,1}','dc_{post,2}','dc_{post,3}','dc_{post,4}','dc_{post,5}','dc_{post,6}','dc_{post,7}','dc_{post,8}','dc_{post,9}','\lambda'};
model.num_para   = length(model.paraID);

% hard bounds, the range for LB, UB, larger than soft bounds
paraH.mu         = [-300,   300]; % ms
paraH.sigma      = [0.01,   300]; % ms
paraH.c1         = [0.01,   350]; % ms
paraH.dc_post1   = [0.01,   300]; % ms
paraH.dc_post2   = [0.01,   300]; % ms
paraH.dc_post3   = [0.01,   300]; % ms
paraH.dc_post4   = [0.01,   300]; % ms
paraH.dc_post5   = [0.01,   300]; % ms
paraH.dc_post6   = [0.01,   300]; % ms
paraH.dc_post7   = [0.01,   300]; % ms
paraH.dc_post8   = [0.01,   300]; % ms
paraH.dc_post9   = [0.01,   300]; % ms
paraH.lambda     = [1e-4,  0.06]; % percentage

% soft bounds, the range for PLB, PUB
paraS.mu         = [-100,   100]; % ms
paraS.sigma      = [  10,   100]; % ms
paraS.c1         = [   1,   100]; % ms
paraS.dc_post1   = [   1,   100]; % ms
paraS.dc_post2   = [   1,   100]; % ms
paraS.dc_post3   = [   1,   100]; % ms
paraS.dc_post4   = [   1,   100]; % ms
paraS.dc_post5   = [   1,   100]; % ms
paraS.dc_post6   = [   1,   100]; % ms
paraS.dc_post7   = [   1,   100]; % ms
paraS.dc_post8   = [   1,   100]; % ms
paraS.dc_post9   = [   1,   100]; % ms
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
    funcNLL          = @(p) nll_gauss_shiftC(p(1), p(2), p(3), p(4), p(5), p(6),p(7),p(8), p(9),...
        p(10),p(11),p(12),p(13),...
        model, data);

    % get grid initializations
    model.num_sections   = model.num_runs*2;
    model.init      = getInit(model.lb, model.ub, model.num_sections, model.num_runs);

    %initialize matrices that store negative log likelihood and best-fit paramters
    minNLL = NaN(1, model.num_runs);
    estimatedP = NaN(model.num_runs, length(model.lb));

%     params = num2cell([model.init(1,:)]);  
%     test = nll_gauss_shiftC(params{:}, model, data);

    parfor i         = 1:model.num_runs
        fprintf('[%s] Start fitting run-%i', mfilename, i);
        try
            tempModel = model;
            [estimatedP(i,:),minNLL(i), ~, ~, OPTIMSTATE(i)] = bads(funcNLL, tempModel.init(i,:), tempModel.lb,...
                tempModel.ub, tempModel.plb, tempModel.pub, [], OPTIONS);
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
    model.OPTIMSTATE = OPTIMSTATE;

   %% model prediction by best-fitting parameters

   model.mode       = 'predict';
   model.test_axis_finer = 1;
   params = num2cell([model.bestP]);  
   pred =  nll_gauss_shiftC(params{:}, model, data);

   %% save by subject

   save(sprintf('sub-%i_%s', i_sub, datestr(datetime('now'))),'data','model','pred')

end

% delete current pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
