% 2210

% This script fits Bayesian causal inference model to the pre- and post-
% TOJ data, for each session, each participant

clear all; close all; clc; rng('Shuffle');
subjID = 1;
sess = 1;

%% organize data
data = organize_data(subjID, sess);

%% define model parameters

% fixed parameters

model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase

model.num_runs                 = 100; %fit the model 100 times, each with a different initialization
model.num_bin              = 100; %#horizontal/vertical bins (used for KDE)
% model.sigma_AV_V_ub           = [5, 12, 20]; %upper bound for sigma_AV_V'
model.thres_R2                = 0.925; %if R2<0.925, we use KDE

% free parameters

% define function calculating nll


nLL = cal_nLL_CI(mu1, sigma1, c1, lambda,... % pretest free para
    sigma2, c2, ... % posttest free para
    p_common, sigma_soa, sigma_c1, sigma_c2, alpha, ... % recal simul para
    model, data); % fixed para, data

% funcNLL = @(p) cal_nLL_CI();

%set a lower bound and an upper bound

%get grid initializations

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