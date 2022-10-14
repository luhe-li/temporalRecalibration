% 2210

% This script fits Bayesian causal inference model to the pre- and post-
% TOJ data, for each session, each participant

% make sure to install bads for optimization

clear all; close all; clc; rng('Shuffle');
subjID = 1;
sess = 1;

%% manage paths
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03analyses')-1); % project folder
addpath(genpath([exptDir '02data'])); % data folder, necessary to run organize_data
addpath(genpath([exptDir '/03analyses/01analysesFunctions'])); % function folder
addpath('/Users/luhe/Documents/GitHub/bads') % add bads to path, to be changed later

%% organize data
data = organize_data(subjID, sess);

%% define model parametersxcv

% set fixed & set-up parameters
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 100; %fit the model 100 times, each with a different initialization
model.num_bin = 100; % num bin to approximate mu_shift
model.thre_r2 = 0.95; %if R2<0.95, we use KDE

% define the grid for free parameters
% hard bounds, the range for LB, UB, larger than soft bounds
paraH.mu1 = [-0.3, 0.3]; % s
paraH.sigma1 = [0.01, 0.35]; %s
paraH.c1 = [0.01, 0.35]; % s
paraH.lambda = [1e-4, 0.06]; % percentage
paraH.sigma2 = [0.01, 0.35]; %s
paraH.c2 = [0.01, 0.35]; % s
paraH.p_common = [1e-4, 1-1e-4]; % weight
paraH.sigma_soa  = [0.01, 0.2]; % s
paraH.sigma_c1 = [1e-4, 0.05]; % s
paraH.sigma_c2  = [0.2, 3]; % s
paraH.alpha = [1e-3, 0.1]; % percentage

% soft bounds, the range for PLB, PUB
paraS.mu1 = [-0.25, 0.25]; % s
paraS.sigma1 = [0.01, 0.3]; %s
paraS.c1 = [0.01, 0.3]; % s
paraS.lambda = [1e-4, 0.06]; % percentage
paraS.sigma2 = [0.01, 0.3]; %s
paraS.c2 = [0.01, 0.3]; % s
paraS.p_common = [0.2, 0.5]; % weight
paraS.sigma_soa  = [0.01, 0.18]; % s
paraS.sigma_c1 = [1e-4, 0.03]; % s
paraS.sigma_c2  = [0.5, 2]; % s
paraS.alpha = [1e-3, 0.01]; % percentage

% reorganize parameter bounds to feed to bads
fn = fieldnames(paraH);
for k = 1:numel(fn)
    model.lb(:,k) = paraH.(fn{k})(1);
    model.ub(:,k) = paraH.(fn{k})(2);
    model.plb(:,k) = paraS.(fn{k})(1);
    model.pub(:,k) = paraS.(fn{k})(2);
end

% set barrier function if needeed
% model.nonbcon

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 1;

%% run bads for optimization

% define function calculating nll
funcNLL = @(p) cal_nLL_CI(p(1), p(2), p(3), p(4), p(5), p(6), p(7), ...
    p(8), p(9), p(10), p(11), model, data);

% get grid initializations
numREachInit    = 1;
numSections     = 3; 
model.init      = gridInitializations_11d(model.plb, model.pub, numSections, ...
                    model.num_runs, numREachInit);     

% %% plug in parameters to test NLL function
% for i = 1:10
%     p = model.init(i,:);
%     testNLL(i) = cal_nLL_CI(p(1), p(2), p(3), p(4), p(5), p(6), p(7), ...
%         p(8), p(9), p(10), p(11), model, data);
% end

%initialize matrices that store negative log likelihood and best-fit paramters
minNLL          = NaN(1, model.num_runs);
estimatedP      = NaN(model.num_runs, length(model.lb));

parfor i = 1:model.num_runs
    disp(i);
    try
        [estimatedP(i,:),minNLL(i)] = bads(funcNLL, model.init(i,:), model.lb,...
            model.ub, model.plb, model.pub, [], OPTIONS); 
        disp(estimatedP(i,:));
        disp(round(minNLL(i),4));
    catch e
        disp('Error!')
    end
end

%% save the data

model.paraID = {'\mu_{pre}','\sigma_{pre}','c_{pre}','\lambda','\sigma_{post}','c_{post}','p_{common}','\sigma_{soa}','\sigma_{c1}','\sigma_{c2}','alpha'};
model.numPara = length(model.paraID);
figure; hold on
set(gcf, 'Position', get(0, 'Screensize'));
set(gca, 'LineWidth', 1, 'FontSize', 15)
for i = 1:model.numPara
    subplot(3,4,i)
    histogram(estimatedP(:,i), 10)
    title(model.paraID{i})
end