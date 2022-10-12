% 2210

% This script fits Bayesian causal inference model to the pre- and post-
% TOJ data, for each session, each participant

clear all; close all; clc; rng('Shuffle');
subjID = 1;
sess = 1;

%% organize data
data = organize_data(subjID, sess);

%% define model parameters

% fixed & set-up parameters
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 100; %fit the model 100 times, each with a different initialization
model.num_bin = 100; % num bin to approximate mu_shift
model.thre_r2 = 0.95; %if R2<0.95, we use KDE

% free parameters
%% debugging: test cal_nLL_CI
mu1 = 0.03; % s
sigma1 = 0.07; %s
c1 = 0.09;
lambda = 0.03;
sigma2 = 0.09;
c2 = 0.12;
p_common = 0.3;
sigma_soa  = 0.16;
sigma_c1 = 0.01;
sigma_c2  = 1;
alpha  = 0.01;

% define function calculating nll
funcNLL = @(p) cal_nLL_CI(p(1), p(2), p(3), p(4), p(5), p(6), p(7), ...
    p(8), p(9), p(10), p(11), model, data);

% set a lower bound and an upper bound

% get grid initializations

%initialize matrices that store negative log likelihood and best-fit paramters
minNLL          = NaN(1, model.numRuns);
estimatedP      = NaN(model.numRuns, length(model.lb));

for i = 1:model.num_runs
    disp(i);
    try
        [estimatedP(i,:),minNLL(i)] = bads(funcNLL,model.init(i,:),model.lb,...
            model.ub,model.plb, model.pub, model.nonbcon); 
        disp(estimatedP(i,:));
        disp(round(minNLL(i),4));
    catch
        disp('Error!')
    end
end

%% save the data